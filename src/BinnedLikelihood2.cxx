/**
 * @file BinnedLikelihood2.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/BinnedLikelihood2.cxx,v 1.3 2011/09/14 00:17:47 jchiang Exp $
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

namespace Likelihood {

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
     m_verbose(true),
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
   create_drm();
   computeCountsSpectrum();
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
     m_verbose(other.m_verbose),
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

   m_nevals++;

   double my_total(m_accumulator.total());

// Contribution from priors.
   std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
   for ( ; par != m_parameter.end(); ++par) {
      my_total += par->log_prior_value();
   }

   saveBestFit(my_total);

   st_stream::StreamFormatter formatter("BinnedLikelihood2", "value", 4);
   formatter.info() << m_nevals << "  "
                    << my_total << "  "
                    << npred << std::endl;

   return my_total;
}

void BinnedLikelihood2::getFreeDerivs(std::vector<double> & derivs) const {
   int nparams(getNumFreeParams());
   derivs.resize(nparams, 0);
   if (nparams > m_posDerivs.size()) {
      for (size_t i(m_posDerivs.size()); i < nparams; i++) {
         m_posDerivs[i] = Accumulator();
         m_negDerivs[i] = Accumulator();
      }
   }
   if (!m_modelIsCurrent) {
      computeModelMap();
   }
   const std::vector<float> & data = m_cmap.data();

   size_t iparam(0);
   std::map<std::string, Source *>::const_iterator src;
   for (src = sources().begin(); src != sources().end(); ++src) {
      std::string srcName = src->second->getName();
      // Skip fixed sources.
      if (std::count(m_fixedSources.begin(), m_fixedSources.end(), srcName)) {
         continue;
      }
      const SourceMap & srcMap(sourceMap(srcName));
      const std::vector<float> & model(srcMap.model());
      std::vector<std::string> paramNames;
      src->second->spectrum().getFreeParamNames(paramNames);
      for (size_t ipar(0); ipar < paramNames.size(); ipar++, iparam++) {
         for (size_t ipix(0); ipix < m_npix; ipix++) {
            // Compute unconvolved counts spectrum at each spatial position.
            std::vector<double> my_derivs(m_nee, 0);
            for (size_t k(0); k < m_nee; k++) {
               double emin(m_energies.at(k));
               double emax(m_energies.at(k+1));
               size_t indx(k*m_npix + ipix);
               if (m_model[ipix][k] != 0) {
                  my_derivs[k] = 
                     src->second->pixelCountsDeriv(emin, emax,
                                                   model.at(indx),
                                                   model.at(indx + m_npix),
                                                   paramNames[ipar]);
               } else {
                  my_derivs[k] = 0;
               }
            }
            // Convolve using the DRM.
            std::vector<double> convolved_derivs(m_nee);
            m_drm->convolve(my_derivs, convolved_derivs);
//            convolved_derivs = my_derivs;
            // Add to Accumulator objects.
            for (size_t k(0); k < m_nee; k++) {
               size_t indx(k*m_npix + ipix);
               double addend(convolved_derivs[k]*
                             (data[indx]/m_model[ipix][k] - 1.));
               if (addend > 0) {
                  m_posDerivs[iparam].add(addend);
               } else {
                  m_negDerivs[iparam].add(addend);
               }
            }
         }
      }
   }
   for (size_t i(0); i < nparams; i++) {
      derivs[i] = m_posDerivs[i].total() + m_negDerivs[i].total();
   }
   // Derivatives from priors.
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

void BinnedLikelihood2::addSource(Source * src, bool fromClone) {
   m_bestValueSoFar = -1e38;
   SourceModel::addSource(src, fromClone);
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

void BinnedLikelihood2::computeCountsSpectrum() {
   m_countsSpectrum.clear();
   size_t indx(0);
   for (size_t k = 0; k < m_nee; k++) {
      double ntot(0);
      for (size_t ipix(0); ipix < m_npix; ipix++, indx++) {
         ntot += m_cmap.data().at(indx);
      }
      m_countsSpectrum.push_back(ntot);
   }
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
      value += src->pixelCounts(energies().at(k), energies().at(k+1),
                                npreds.at(k), npreds.at(k+1));
   }
   return value;
}

double BinnedLikelihood2::computeModelMap() const {
   double npred(0);

   ModelWeights_t modelWts;
   const_cast<BinnedLikelihood2 *>(this)->compute_model_wts(modelWts);

   size_t ipix(0);
   for (size_t ipix(0); ipix < m_npix; ipix++) {
      std::vector<double> true_counts(m_nee);
      for (size_t k(0); k < m_nee; k++) {
         true_counts[k] = pixelCounts(m_energies[k], m_energies[k+1],
                                      modelWts[ipix][k].first, 
                                      modelWts[ipix][k].second);
      }
      m_drm->convolve(true_counts, m_model[ipix]);
//      m_model[ipix] = true_counts;
      for (size_t k(0); k < m_nee; k++) {
         npred += m_model[ipix][k];
      }
   }
   m_modelIsCurrent = true;

   return npred;
}

void BinnedLikelihood2::compute_model_wts(ModelWeights_t & modelWts) {
   if (fixedModelUpdated()) {
      buildFixedModelWts();
   }
   modelWts = m_fixedModelWts;
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (size_t i(0); i < srcNames.size(); i++) {
      if (std::count(m_fixedSources.begin(), m_fixedSources.end(),
                     srcNames[i]) == 0) {
         addSourceWts(modelWts, srcNames[i]);
      }
   }
}

bool BinnedLikelihood2::fixedModelUpdated() const {
// Check if the fixed list has changed.
   std::map<std::string, Source *>::const_iterator it(m_sources.begin());
   std::vector<std::string> fixedSources;
   for ( ; it != m_sources.end(); ++it) {
      if (it->second->fixedSpectrum()) {
         fixedSources.push_back(it->first);
         // Check for additions.
         if (std::count(m_fixedSources.begin(), m_fixedSources.end(), 
                        fixedSources.back()) != 1) {
            return true;
         }
      }
   }
   // Check for deletions.
   if (fixedSources.size() != m_fixedSources.size()) {
      return true;
   }
   
// Compare current parameter values for fixed sources with saved values.
   for (it = m_sources.begin(); it != m_sources.end(); ++it) {
      if (it->second->fixedSpectrum()) {
         const std::string & name(it->first);
         std::map<std::string, std::vector<double> >::const_iterator
            savedPars = m_modelPars.find(name);
         if (savedPars == m_modelPars.end()) {
            throw std::runtime_error("BinnedLikelihood2::fixedModelUpdated: "
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

void BinnedLikelihood2::buildFixedModelWts() {
   m_fixedSources.clear();

   m_fixedModelWts.clear();
   std::vector<std::pair<double, double> > zeros(m_nee, std::make_pair(0, 0));
   m_fixedModelWts.resize(m_npix, zeros);

   m_fixedNpreds.clear();
   m_fixedNpreds.resize(m_energies.size(), 0);

   std::map<std::string, Source *>::const_iterator it(m_sources.begin());
   for ( ; it != m_sources.end(); ++it) {
      const std::string & srcName(it->first);
      if (it->second->fixedSpectrum()) {
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
            m_fixedNpreds[k] += (it->second->spectrum()(ee)*
                                 srcMap->npreds()[k]);
         }
         // Recover memory by deleting the fixed SourceMap.
         if (srcMapIt != m_srcMaps.end()) {
            m_srcMaps.erase(srcName);
         }
         delete srcMap;
         srcMap = 0;
         // Save the parameter values for later comparison to see if
         // values have been changed by the user by hand.
         std::vector<double> pars;
         it->second->spectrum().getParamValues(pars);
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
         // Delete any lingering Npred values from fixed map since
         // they must be computed on-the-fly for non-fixed sources.
         if (m_fixedModelNpreds.find(srcName) != m_fixedModelNpreds.end()) {
            m_fixedModelNpreds.erase(srcName);
         }
      }
   }
}

void BinnedLikelihood2::
addSourceWts(ModelWeights_t & modelWts, const std::string & srcName,
             const SourceMap * srcMap, bool subtract) const {
   const Source * src(const_cast<BinnedLikelihood2 *>(this)->getSource(srcName));
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
         modelWts[ipix][k].first += (my_sign*model[indx]*
                                     spectrum(src, emin));
         modelWts[ipix][k].second += (my_sign*model[indx + m_npix]*
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

double BinnedLikelihood2::NpredValue(const std::string & srcName,
                                     const SourceMap & sourceMap) const {
   const std::vector<double> & npreds(sourceMap.npreds());
   const Source * src(const_cast<BinnedLikelihood2 *>(this)->getSource(srcName));
   double value(0);
   for (size_t k(0); k < energies().size()-1; k++) {
      value += src->pixelCounts(energies().at(k), energies().at(k+1),
                                npreds.at(k), npreds.at(k+1));
   }
   return value;
}

SourceMap * BinnedLikelihood2::getSourceMap(const std::string & srcName, 
                                            bool verbose) const {
   if (fileHasSourceMap(srcName, m_smaps_file)) {
      return new SourceMap(m_smaps_file, srcName);
   }
// Generate the map if it is not in the file.   
   Source * src(const_cast<BinnedLikelihood2 *>(this)->getSource(srcName));
   if (src->getType() == "Diffuse" || m_computePointSources) {
      return new SourceMap(src, &m_cmap, m_observation,
                           m_applyPsfCorrections,
                           m_performConvolution,
                           m_resample, m_resamp_factor,
                           m_minbinsz, verbose && m_verbose);
   }
   return 0;
}

bool BinnedLikelihood2::fileHasSourceMap(const std::string & srcName,
                                        const std::string & fitsFile) const {
   try {
      std::auto_ptr<const tip::Image> 
         image(tip::IFileSvc::instance().readImage(fitsFile, srcName));
   } catch (tip::TipException &) {
      return false;
   }
   return true;
}

} // namespace Likelihood
