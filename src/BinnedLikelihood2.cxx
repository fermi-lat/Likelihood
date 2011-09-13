/**
 * @file BinnedLikelihood2.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/BinnedLikelihood2.cxx,v 1.1 2011/06/27 05:05:38 jchiang Exp $
 */

#include <cmath>

#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedLikelihood2.h"
#include "Likelihood/Drm.h"
#include "Likelihood/CountsMap.h"
#define ST_DLL_EXPORTS
#include "Likelihood/SourceMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/SourceModel.h"

//namespace Likelihood {

BinnedLikelihood2::BinnedLikelihood2(const CountsMap & cmap,
                                     const Observation & observation,
                                     const std::string & smaps_file,
                                     bool computePointSources,
                                     bool applyPsfCorrections,
                                     bool performConvolution, 
                                     bool resample,
                                     double resamp_factor,
                                     double minbinsz)
   : LogLike(observation), 
     m_cmap(cmap), 
     m_smaps_file(smaps_file),
     m_computePointSources(computePointSources),
     m_applyPsfCorrections(applyPsfCorrections),
     m_performConvolution(performConvolution), 
     m_resample(resample), 
     m_resamp_factor(resamp_factor), 
     m_minbinsz(minbinsz),
     m_drm(0), 
     m_modelIsCurrent(false),
     m_npix(cmap.imageDimension(0)*cmap.imageDimension(1)),
     m_nee(cmap.imageDimension(2)) {
   cmap.getAxisVector(2, m_energies);
   m_model.resize(m_npix);
   m_fixedModelWts.resize(m_npix);
   for (size_t ipix(0); ipix < m_npix; ipix++) {
      m_model[ipix].resize(m_nee);
      m_fixedModelWts[ipix].resize(m_nee);
   }
}

BinnedLikelihood2::~BinnedLikelihood2() throw() {
   std::map<std::string, SourceMap *>::iterator srcMap(m_srcMaps.begin());
   for ( ; srcMap != m_srcMaps.end(); ++srcMap) {
      delete srcMap->second;
   }
   delete m_drm;
}

BinnedLikelihood2::BinnedLikelihood2(const BinnedLikelihood2 & other) 
   : LogLike(other), 
     m_cmap(other.m_cmap), 
     m_smaps_file(other.m_smaps_file),
     m_computePointSources(other.m_computePointSources),
     m_applyPsfCorrections(other.m_applyPsfCorrections),
     m_performConvolution(other.m_performConvolution),
     m_resample(other.m_resample),
     m_resamp_factor(other.m_resamp_factor),
     m_minbinsz(other.m_minbinsz),
     m_drm(new Drm(*other.m_drm)),
     m_modelIsCurrent(other.m_modelIsCurrent),
     m_npix(other.m_npix),
     m_nee(other.m_nee),
     m_energies(other.m_energies),
     m_model(other.m_model),
     m_fixedSources(other.m_fixedSources),
     m_fixedModelWts(other.m_fixedModelWts),
     m_posDerivs(other.m_posDerivs),
     m_negDerivs(other.m_negDerivs) {
   typedef std::map<std::string, SourceMap *> SourceMaps_t;
   SourceMaps_t::const_iterator it(other.m_srcMaps.begin());
   for ( ; it != other.m_srcMaps.end(); ++it) {
      m_srcMaps[it->first] = new SourceMap(*(it->second));
   }
}

double BinnedLikelihood2::value(optimizers::Arg & dummy) const {
   (void)(dummy);

   double npred(computeModelMap());

   const std::vector<float> & data(m_cmap.data());

   size_t indx(0);
   for (size_t k(0); k < m_nee; k++) {
      for (size_t ipix(0); ipix < m_npix; ipix++, indx++) {
         if (data[indx] > 0 && m_model[ipix][k] > 0) {
            double addend(data[indx]*std::log(m_model[ipix][k]));
            m_accumulator.add(addend);
         }
      }
   }
   m_accumulator.add(-npred);

   st_stream::StreamFormatter formatter("BinnedLikelihood2", "value", 4);
   formatter.info() << m_nevals << "  "
                    << npred << std::endl;
   m_nevals++;

   double my_total(m_accumulator.total());

// Contribution from priors.
   std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
   for ( ; par != m_parameter.end(); ++par) {
      my_total += par->log_prior_value();
   }

   saveBestFit(my_total);
   return my_total;
}

void BinnedLikelihood2::getFreeDerivs(std::vector<double> & derivs) const {
   int nparams(getNumFreeParams());
   derivs.resize(nparams, 0);
   double npred;
   if (!m_modelIsCurrent) {
      computeModelMap(npred);
   }
   const std::vector<float> & data = m_cmap.data();
   // First j value corresponding to the minimum allowed k-index.
   size_t jentry(0);
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t k(m_filledPixels[j]/m_pixels.size());
      if (k == m_kmin) {
         jentry = j;
         break;
      }
   }
   for (size_t j(jentry); j < m_filledPixels.size(); j++) {
      size_t jmin(m_filledPixels.at(j));
      size_t jmax(jmin + m_pixels.size());
      size_t k(jmin/m_pixels.size());
      if (k < m_kmin || k > m_kmax-1) {
         continue;
      }
      double emin(m_energies.at(k));
      double emax(m_energies.at(k+1));
      if (m_model.at(j) > 0) {
         long iparam(0);
         if (m_posDerivs.find(iparam) == m_posDerivs.end()) {
            m_posDerivs[iparam] = Accumulator();
            m_negDerivs[iparam] = Accumulator();
         }
         const std::map<std::string, Source *> & my_srcs = sources();
         std::map<std::string, Source *>::const_iterator src;
         for (src = my_srcs.begin(); src != my_srcs.end(); ++src) {
            std::string srcName = src->second->getName();
            if (std::count(m_fixedSources.begin(), m_fixedSources.end(),
                           srcName)) {
               continue;
            }
            const SourceMap & srcMap = sourceMap(srcName);
            const std::vector<float> & model = srcMap.model();
            std::vector<std::string> paramNames;
            src->second->spectrum().getFreeParamNames(paramNames);
            for (size_t i(0); i < paramNames.size(); i++, iparam++) {
               double my_deriv = 
                  src->second->pixelCountsDeriv(emin, emax,
                                                model.at(jmin), 
                                                model.at(jmax),
                                                paramNames.at(i));
               double addend(data.at(jmin)/m_model.at(j)*my_deriv);
               if (addend > 0) {
                  m_posDerivs[iparam].add(addend);
               } else {
                  m_negDerivs[iparam].add(addend);
               }
               if (j == jentry) {
                  const std::vector<double> & npreds = srcMap.npreds();
                  for (size_t kk(0); kk < m_energies.size()-1; kk++) {
                     if (kk >= m_kmin && kk <= m_kmax-1) {
                        addend = 
                           src->second->pixelCountsDeriv(m_energies.at(kk), 
                                                         m_energies.at(kk+1),
                                                         npreds.at(kk), 
                                                         npreds.at(kk+1),
                                                         paramNames.at(i));
                        if (-addend > 0) {
                           m_posDerivs[iparam].add(-addend);
                        } else {
                           m_negDerivs[iparam].add(-addend);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   for (size_t i(0); i < derivs.size(); i++) {
      derivs[i] = m_posDerivs[i].total() + m_negDerivs[i].total(); 
   }

   /// Derivatives from priors.
   size_t i(0);
   std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
   for ( ; par != m_parameter.end(); ++par) {
      if (par->isFree()) {
         derivs[i] += par->log_prior_deriv();
         i++;
      }
   }
}

void BinnedLikelihood2::readXml(std::string xmlFile, 
                               optimizers::FunctionFactory & funcFactory,
                               bool requireExposure, 
                               bool addPointSources, 
                               bool loadMaps,
                               bool createAllMaps) {
   SourceModel::readXml(xmlFile, funcFactory, requireExposure=false,
                        addPointSources, loadMaps);
}

std::vector<double>::const_iterator 
BinnedLikelihood2::setParamValues_(std::vector<double>::const_iterator it) {
   m_modelIsCurrent = false;
   return SourceModel::setParamValues_(it);
}

std::vector<double>::const_iterator 
BinnedLikelihood2::setFreeParamValues_(std::vector<double>::const_iterator it) {
   m_modelIsCurrent = false;
   return SourceModel::setFreeParamValues_(it);
}

void BinnedLikelihood2::addSource(Source * src) {
   m_bestValueSoFar = -1e38;
   SourceModel::addSource(src);
   if (m_srcMaps.find(src->getName()) == m_srcMaps.end()) {
      SourceMap * srcMap(getSourceMap(src->getName(), true));
      if (srcMap) {
         m_srcMaps[src->getName()] = srcMap;
      }
   }
   std::vector<double> pars;
   src->spectrum().getParamValues(pars);
   m_modelPars[src->getName()] = pars;
}

Source * BinnedLikelihood2::deleteSource(const std::string & srcName) {
   m_bestValueSoFar = -1e38;
// Check if this is a fixed source, and if so, remove it from the fixed model.
   std::vector<std::string>::iterator srcIt = 
      std::find(m_fixedSources.begin(), m_fixedSources.end(), srcName);
   if (source(srcName).fixedSpectrum() && srcIt != m_fixedSources.end()) {
      SourceMap * srcMap(getSourceMap(srcName, false));
      bool subtract;
      addSourceWts(m_fixedModelWts, srcName, srcMap, subtract=true);
      m_fixedModelNpreds.erase(srcName);
      m_fixedSources.erase(srcIt);
   }
   Source * src(SourceModel::deleteSource(srcName));
   return src;
}


void BinnedLikelihood2::syncParams() {
   SourceModel::syncParams();
}

void BinnedLikelihood2::syncSrcParams(const std::string & srcName) {
   (void)(srcName);
   syncParams();
   m_modelIsCurrent = false;
}

double BinnedLikelihood2::NpredValue(const std::string & srcName) const {
   std::map<std::string, double>::const_iterator npredIt 
      = m_fixedModelNpreds.find(srcName);
   if (npredIt != m_fixedModelNpreds.end()) {
      return npredIt->second;
   }
   const std::vector<double> & npreds(sourceMap(srcName).npreds());
   const Source * src(const_cast<BinnedLikelihood2 *>(this)->getSource(srcName));
   double value(0);
   for (size_t k(0); k < energies().size()-1; k++) {
      if (k < m_kmin || k > m_kmax-1) {
         continue;
      }
      value += src->pixelCounts(energies().at(k), energies().at(k+1),
                                npreds.at(k), npreds.at(k+1));
   }
   return value;
}

double BinnedLikelihood2::computeModelMap() const {
   double npred(0);

   ModelWeights_t modelWts;
   compute_model_wts(modelWts);

   size_t ipix(0);
   for (size_t j(0); j < m_ny; j++) {
      for (size_t i(0); i < m_nx; i++, ipix++) {
         std::vector<double> true_counts(nz);
         for (size_t k(0); k < m_nz; k++) {
            true_counts[k] = pixelCounts(m_energies[k], m_energies[k+1],
                                         modelWts[ipix][k].first, 
                                         modelWts[ipix][k].second);
         }
         m_drm->convolve(true_counts, m_model[ipix]);
         for (size_t k(0); k < m_nz; k++) {
            npred += m_model[ipix][k];
         }
      }
   }
   m_modelIsCurrent = true;
}

void BinnedLikelihood2::compute_model_wts(ModelWeights_t & modelWts) {
   if (fixedModelUpdated()) {
      buildFixedModelWts();
   }
   modelWts = m_fixedModelWts;
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (size_t i(0); i < srcNames.size(); i++) {
      if (std::count(m_fixedSources.begin(), m_fixedSource.end(),
                     srcNames[i]) == 0) {
         addSourceWts(modelWts, srcNames[i]);
      }
   }
}

bool BinnedLikelihood::fixedModelUpdated() const {
// Check if the fixed list has changed.
   std::map<std::string, Source *>::const_iterator it(m_sources.begin());
   std::vector<std::string> fixedSources;
   for ( ; it != m_sources.end(); ++it) {
      if (it->second->fixedSpectrum()) {
         fixedSources.push_back(it->first);
         if (std::count(m_fixedSources.begin(), m_fixedSources.end(), 
                        fixedSources.back()) != 1) {
            return true;
         }
      }
   }
   if (fixedSources.size() != m_fixedSources.size()) {
      return true;
   }
   
// Compare current parameter values for fixed sources with saved
   for (it = m_sources.begin(); it != m_sources.end(); ++it) {
      if (it->second->fixedSpectrum()) {
         const std::string & name(it->first);
         std::map<std::string, std::vector<double> >::const_iterator
            savedPars = m_modelPars.find(name);
         if (savedPars == m_modelPars.end()) {
            throw std::runtime_error("BinnedLikelihood::fixedModelUpdated: "
                                     "inconsistent m_modelPars.");
         }
         std::vector<double> parValues;
         it->second->spectrum().getParamValues(parValues);
         for (size_t j(0); j < parValues.size(); j++) {
            if (parValues.at(j) != savedPars->second.at(j)) {
               return true;
            }
         }
      }
   }
   return false;
}

void BinnedLikelihood::buildFixedModelWts() {
   m_fixedSources.clear();
   m_fixedModelWts.clear();
   m_fixedModelWts.resize(m_filledPixels.size(), std::make_pair(0, 0));
   m_fixedNpreds.clear();
   m_fixedNpreds.resize(m_energies.size(), 0);
   std::map<std::string, Source *>::const_iterator srcIt(m_sources.begin());
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      const std::string & srcName(srcIt->first);
      if (srcIt->second->fixedSpectrum()) {
         m_fixedSources.push_back(srcName);
         SourceMap * srcMap(0);
         std::map<std::string, SourceMap *>::const_iterator srcMapIt
            = m_srcMaps.find(srcName);
         if (srcMapIt == m_srcMaps.end()) {
            srcMap = getSourceMap(srcName, false);
         } else {
            srcMap = srcMapIt->second;
         }
         addSourceWts(m_fixedModelWts, srcName, srcMap);
         m_fixedModelNpreds[srcName] = NpredValue(srcName, *srcMap);
         for (size_t k(0); k < m_energies.size(); k++) {
            optimizers::dArg ee(m_energies[k]);
            m_fixedNpreds[k] += (srcIt->second->spectrum()(ee)*
                                 srcMap->npreds()[k]);
         }
         if (srcMapIt != m_srcMaps.end()) {
            m_srcMaps.erase(srcName);
         }
         delete srcMap;
         srcMap = 0;
         std::vector<double> pars;
         srcIt->second->spectrum().getParamValues(pars);
         m_modelPars[srcName] = pars;
      } else { 
// Process non-fixed sources.
//
// Ensure model map is available.
         std::map<std::string, SourceMap *>::const_iterator srcMapIt
            = m_srcMaps.find(srcName);
         if (srcMapIt == m_srcMaps.end()) {
            SourceMap * srcMap(getSourceMap(srcName, false));
            if (srcMap) {
               m_srcMaps[srcName] = srcMap;
            }
         }
// Delete any lingering Npred values from fixed map since they must be
// computed on-the-fly for non-fixed sources.
         if (m_fixedModelNpreds.find(srcName) != m_fixedModelNpreds.end()) {
            m_fixedModelNpreds.erase(srcName);
         }
      }
   }
}

void BinnedLikelihood2::
addSourceWts(ModelWeights_t & modelWts,
             const std::string & srcName,
             const SourceMap * srcMap,
             bool subtract) const {
   const Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
   if (src == 0) {
      return;
   }
   if (srcMap == 0) {
      std::map<std::string, SourceMap *>::const_iterator it 
         = m_srcMaps.find(srcName);
      if (it == m_srcMaps.end()) {
         throw std::runtime_error("BinnedLikelihood2::addSourceWts: "
                                  "source not found in m_srcMaps.");
      }
      srcMap = it->second;
   }
   double my_sign(1.);
   if (subtract) {
      my_sign = -1.;
   }
   const std::vector<float> & model(srcMap->model());
   size_t indx(0);
   for (size_t k(0); k < m_nee; k++) {
      double emin(m_energies[k]);
      double emax(m_energies[k+1]);
      for (size_t ipix(0); ipix < m_npix; ipix++, indx++) {
         m_modelWts[ipix][k].first += (my_sign*model[indx]*
                                       spectrum(src, emin));
         m_modelWts[ipix][k].second += (my_sign*model[indx+npix]*
                                        spectrum(src, emax));
      }
   }
}

double BinnedLikelihood2::spectrum(const Source * src, double energy) const {
   optimizers::dArg eArg(energy);
   return src->spectrum().operator()(eArg);
}

double BinnedLikelihood2::pixelCounts(double emin, double emax,
                                      double y1, double y2) const {
   return (y1*emin + y2*emax)/2.*std::log(emax/emin);
}

void BinnedLikelihood2::create_drm() {
   const astro::SkyDir & ref_dir(m_cmap.refDir());
   m_drm = new Drm(ref_dir.ra(), ref_dir.dec(), observation(), 
                   m_cmap.energies());
}

} // namespace Likelihood
