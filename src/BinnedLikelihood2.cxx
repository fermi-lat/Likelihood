/**
 * @file BinnedLikelihood2.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/BinnedLikelihood2.cxx,v 1.79 2011/04/13 18:10:32 jchiang Exp $
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

BinnedLikelihood2::BinnedLikelihood2(const CountsMap & dataMap,
                                   const Observation & observation,
                                   const std::string & srcMapsFile,
                                   bool computePointSources,
                                   bool applyPsfCorrections,
                                   bool performConvolution, 
                                   bool resample,
                                   double resamp_factor,
                                   double minbinsz)
   : LogLike(observation), m_dataMap(dataMap), m_pixels(dataMap.pixels()),
     m_modelIsCurrent(false), m_srcMapsFile(srcMapsFile),
     m_computePointSources(computePointSources),
     m_applyPsfCorrections(applyPsfCorrections),
     m_performConvolution(performConvolution), 
     m_resample(resample), 
     m_resamp_factor(resamp_factor), 
     m_minbinsz(minbinsz),
     m_kmin(0), m_drm(0), m_apply_drm(false) {
   dataMap.getAxisVector(2, m_energies);
   m_kmax = m_energies.size() - 1;
   identifyFilledPixels();
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
     m_dataMap(other.m_dataMap), 
     m_pixels(other.m_pixels),
     m_energies(other.m_energies),
     m_countsSpectrum(other.m_countsSpectrum),
     m_filledPixels(other.m_filledPixels),
     m_model(other.m_model),
     m_modelIsCurrent(other.m_modelIsCurrent), 
     m_srcMapsFile(other.m_srcMapsFile),
     m_computePointSources(other.m_computePointSources),
     m_applyPsfCorrections(other.m_applyPsfCorrections),
     m_performConvolution(other.m_performConvolution),
     m_resample(other.m_resample),
     m_minbinsz(other.m_minbinsz),
     m_fixedSources(other.m_fixedSources),
     m_fixedModelWts(other.m_fixedModelWts),
     m_fixedModelNpreds(other.m_fixedModelNpreds),
     m_modelPars(other.m_modelPars),
     m_posDerivs(other.m_posDerivs),
     m_negDerivs(other.m_negDerivs),
     m_kmin(other.m_kmin),
     m_kmax(other.m_kmax),
     m_fixedNpreds(other.m_fixedNpreds),
     m_apply_drm(other.m_apply_drm) {
   typedef std::map<std::string, SourceMap *> SourceMaps_t;
   SourceMaps_t::const_iterator it(other.m_srcMaps.begin());
   for ( ; it != other.m_srcMaps.end(); ++it) {
      m_srcMaps[it->first] = new SourceMap(*(it->second));
   }
   m_drm = new Drm(*other.m_drm);
}

double BinnedLikelihood2::value(optimizers::Arg & dummy) const {
   (void)(dummy);

   double npred;
   computeModelMap(npred);

   const std::vector<float> & data(m_dataMap.data());
   double my_value(0);

   for (size_t i(0); i < m_filledPixels.size(); i++) {
      if (m_model.at(i) > 0) {
         size_t j(m_filledPixels.at(i));
         size_t k(j/m_pixels.size());
         if (k >= m_kmin && k <= m_kmax) {
            double addend = data.at(j)*std::log(m_model[i]);
            my_value += addend;
            m_accumulator.add(addend);
         }
      }
   }
   my_value -= npred;
   m_accumulator.add(-npred);

   st_stream::StreamFormatter formatter("BinnedLikelihood2", "value", 4);
   formatter.info() << m_nevals << "  "
                    << my_value << "  "
                    << npred << std::endl;
   m_nevals++;

   double my_total(m_accumulator.total());

/// Add in contribution from priors.
   std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
   for ( ; par != m_parameter.end(); ++par) {
      my_total += par->log_prior_value();
   }

   saveBestFit(my_total);
   return my_total;
}

double BinnedLikelihood2::npred() {
   double npred;
   computeModelMap(npred);
   return npred;
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

void BinnedLikelihood2::getFreeDerivs(std::vector<double> & derivs) const {
   int nparams(getNumFreeParams());
   derivs.resize(nparams, 0);
   double npred;
   if (!m_modelIsCurrent) {
      computeModelMap(npred);
   }
   const std::vector<float> & data = m_dataMap.data();
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

CountsMap * BinnedLikelihood2::createCountsMap() const {
   std::vector<float> map;
   computeModelMap(map);

   CountsMap * modelMap = new CountsMap(m_dataMap);
         
   modelMap->setImage(map);
   return modelMap;
}

void BinnedLikelihood2::readXml(std::string xmlFile, 
                               optimizers::FunctionFactory & funcFactory,
                               bool requireExposure, 
                               bool addPointSources, 
                               bool loadMaps,
                               bool createAllMaps) {
   SourceModel::readXml(xmlFile, funcFactory, requireExposure=false,
                        addPointSources, loadMaps);
   if (m_srcMapsFile == "" || createAllMaps) {
// No need to do anything, since the maps are created in the the call to 
// SourceModel::readXml via the overloaded addSource(...) call therein.
   } else {
      buildFixedModelWts();
   }
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

void BinnedLikelihood2::computeModelMap(std::vector<float> & modelMap) const {
   modelMap.clear();
   modelMap.resize(m_pixels.size()*(m_energies.size()-1), 0);
   for (size_t k(0); k < m_energies.size()-1; k++) {
      double emin(m_energies.at(k));
      double emax(m_energies.at(k+1));
      for (size_t j(0); j < m_pixels.size(); j++) {
         size_t imin(k*m_pixels.size() + j);
         size_t imax(imin + m_pixels.size());
         std::map<std::string, SourceMap *>::const_iterator srcIt
            = m_srcMaps.begin();
         for ( ; srcIt != m_srcMaps.end(); ++srcIt) {
            const std::string & name(srcIt->first);
            const SourceMap * srcMap(srcIt->second);
            const Source * src 
               = const_cast<BinnedLikelihood2 *>(this)->getSource(name);
            const std::vector<float> & model(srcMap->model());
            modelMap.at(imin) += src->pixelCounts(emin, emax,
                                                  model.at(imin), 
                                                  model.at(imax));
         }
      }
   }
}

bool BinnedLikelihood2::fixedModelUpdated() const {
// Check if the fixed list has changed.
   std::map<std::string, Source *>::const_iterator srcIt(m_sources.begin());
   std::vector<std::string> fixedSources;
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      if (srcIt->second->fixedSpectrum()) {
         fixedSources.push_back(srcIt->first);
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
   for (srcIt = m_sources.begin(); srcIt != m_sources.end(); ++srcIt) {
      if (srcIt->second->fixedSpectrum()) {
         const std::string & name(srcIt->first);
         std::map<std::string, std::vector<double> >::const_iterator
            savedPars = m_modelPars.find(name);
         if (savedPars == m_modelPars.end()) {
            throw std::runtime_error("BinnedLikelihood2::fixedModelUpdated: "
                                     "inconsistent m_modelPars.");
         }
         std::vector<double> parValues;
         srcIt->second->spectrum().getParamValues(parValues);
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

void BinnedLikelihood2::computeModelMap(double & npred) const {
   npred = 0;

   std::vector<std::pair<double, double> > modelWts;
   modelWts.resize(m_filledPixels.size());

   if (fixedModelUpdated()) {
      const_cast<BinnedLikelihood2 *>(this)->buildFixedModelWts();
   }

   for (size_t j(0); j < m_fixedModelWts.size(); j++) {
      modelWts.at(j).first = m_fixedModelWts.at(j).first;
      modelWts.at(j).second = m_fixedModelWts.at(j).second;
   }

   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (size_t i(0); i < srcNames.size(); i++) {
      if (std::count(m_fixedSources.begin(), m_fixedSources.end(),
                     srcNames.at(i)) == 0) {
         addSourceWts(modelWts, srcNames[i]);
      }
      npred += NpredValue(srcNames[i]);
   }

   m_model.clear();
   m_model.resize(m_filledPixels.size(), 0);
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t k(m_filledPixels.at(j)/m_pixels.size());
      if (k < m_kmin || k > m_kmax-1) {
         continue;
      }
      double emin(m_energies[k]);
      double emax(m_energies[k+1]);
      m_model.at(j) = pixelCounts(emin, emax, modelWts.at(j).first,
                                  modelWts.at(j).second);
   }

   m_modelIsCurrent = true;
}

void BinnedLikelihood2::
addSourceWts(std::vector<std::pair<double, double> > & modelWts,
             const std::string & srcName,
             const SourceMap * srcMap,
             bool subtract) const {
   const Source * src(const_cast<BinnedLikelihood2 *>(this)->getSource(srcName));
   if (src == 0) {
      return;
   }
   if (modelWts.size() != m_filledPixels.size()) {
      throw std::runtime_error("BinnedLikelihood2::addSourceWts: "
                               "modelWts size does not match "
                               "number of filled pixels.");
   }
   const SourceMap * sourceMap;
   if (srcMap != 0) {
      sourceMap = srcMap;
   } else {
      std::map<std::string, SourceMap *>::const_iterator srcIt 
         = m_srcMaps.find(srcName);
      if (srcIt == m_srcMaps.end()) {
         throw std::runtime_error("BinnedLikelihood2::addSourceWts: "
                                  "source not found in m_srcMaps.");
      }
      srcMap = srcIt->second;
   }

   double my_sign(1.);
   if (subtract) {
      my_sign = -1.;
   }
      
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t jmin(m_filledPixels.at(j));
      size_t jmax(jmin + m_pixels.size());
      size_t k(jmin/m_pixels.size());
      double emin(m_energies.at(k));
      double emax(m_energies.at(k+1));
      const std::vector<float> & model(srcMap->model());
      modelWts.at(j).first += my_sign*model.at(jmin)*spectrum(src, emin);
      modelWts.at(j).second += my_sign*model.at(jmax)*spectrum(src, emax);
   }
}

double BinnedLikelihood2::spectrum(const Source * src, double energy) const {
   const optimizers::Function & spectrum(src->spectrum());
   optimizers::dArg eArg(energy);
   return spectrum(eArg);
}

double BinnedLikelihood2::pixelCounts(double emin, double emax,
                                     double y1, double y2) const {
   if (::getenv("USE_LINEAR_QUADRATURE")) {
      return (y1 + y2)*(emax - emin)/2.;
   }
   return (y1*emin + y2*emax)/2.*std::log(emax/emin);
}

void BinnedLikelihood2::createSourceMaps() {
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   std::vector<std::string>::const_iterator name = srcNames.begin();
   for ( ; name != srcNames.end(); ++name) {
      Source * src = getSource(*name);
      if (src->getType() == "Diffuse" || m_computePointSources) {
         m_srcMaps[*name] = createSourceMap(*name);
      }
   }
}

SourceMap * BinnedLikelihood2::createSourceMap(const std::string & srcName) {
   Source * src = getSource(srcName);
   return new SourceMap(src, &m_dataMap, m_observation, m_applyPsfCorrections,
                        m_performConvolution, m_resample, m_resamp_factor,
                        m_minbinsz);
}

SourceMap * BinnedLikelihood2::getSourceMap(const std::string & srcName, 
                                           bool verbose) const {
   if (fileHasSourceMap(srcName, m_srcMapsFile)) {
      return new SourceMap(m_srcMapsFile, srcName);
   }
// Generate the map if it is not in the file.   
   Source * src(const_cast<BinnedLikelihood2 *>(this)->getSource(srcName));
   if (src->getType() == "Diffuse" || m_computePointSources) {
      return new SourceMap(src, &m_dataMap, m_observation,
                           m_applyPsfCorrections,
                           m_performConvolution,
                           m_resample, m_resamp_factor,
                           m_minbinsz, verbose);
   }
   return 0;
}

void BinnedLikelihood2::saveSourceMaps(const std::string & filename) {
   if (filename != "") {
      m_srcMapsFile = filename;
   }
   if (!st_facilities::Util::fileExists(filename)) {
      m_dataMap.writeOutput("BinnedLikelihood2", m_srcMapsFile);
   }
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      st_stream::StreamFormatter formatter("BinnedLikelihood2",
                                           "saveSourceMaps", 4);
      formatter.info() << srcNames.at(i) << std::endl;
      if (m_srcMaps.count(srcNames.at(i))) {
         if (fileHasSourceMap(srcNames.at(i), m_srcMapsFile)) {
//            replaceSourceMap(srcNames.at(i), m_srcMapsFile);
         } else {
            formatter.info() << "appending map for " 
                             << srcNames.at(i) << std::endl;
            appendSourceMap(srcNames.at(i), m_srcMapsFile);
         }
      }
   }
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

void BinnedLikelihood2::replaceSourceMap(const std::string & srcName,
                                        const std::string & fitsFile) const {
   const SourceMap & srcMap = sourceMap(srcName);

   std::auto_ptr<tip::Image> 
      image(tip::IFileSvc::instance().editImage(fitsFile, srcName));

   image->set(srcMap.model());
}

void BinnedLikelihood2::appendSourceMap(const std::string & srcName,
                                       const std::string & fitsFile) const {
   if (!m_srcMaps.count(srcName)) {
      throw std::runtime_error("BinnedLikelihood2::appendSourceMap: " +
                               std::string("Source ") + srcName 
                               + " is not available.");
   }

   std::vector<long> naxes;
   naxes.push_back(m_dataMap.imageDimension(0));
   naxes.push_back(m_dataMap.imageDimension(1));
   naxes.push_back(m_energies.size());

   tip::IFileSvc::instance().appendImage(fitsFile, srcName, naxes);
   tip::Image * image = tip::IFileSvc::instance().editImage(fitsFile, srcName);
      
   image->set(m_srcMaps.find(srcName)->second->model());

   m_dataMap.setKeywords(image->getHeader());

   delete image;
}

void BinnedLikelihood2::setImageDimensions(tip::Image * image, 
                                          long * dimensions) const {
   typedef std::vector<tip::PixOrd_t> DimCont_t;
   DimCont_t dims = image->getImageDimensions();
   DimCont_t::size_type num_dims = dims.size();
   if (3 != num_dims) {
      throw std::runtime_error("BinnedLikelihood2::setImageDimensions: "
                               + std::string("Cannot write a SourceMap file ")
                               + "to an image which is not 3D");
   }
   for (DimCont_t::size_type index = 0; index != num_dims; ++index) {
      dims[index] = dimensions[index];
   }
   image->setImageDimensions(dims);
}

void BinnedLikelihood2::identifyFilledPixels() {
   const std::vector<float> & data = m_dataMap.data();
   m_filledPixels.clear();
   for (unsigned int i = 0; i < data.size(); i++) {
      if (data.at(i) > 0) {
         m_filledPixels.push_back(i);
      }
   }
}

void BinnedLikelihood2::computeCountsSpectrum() {
   m_countsSpectrum.clear();
   size_t nx(m_dataMap.imageDimension(0));
   size_t ny(m_dataMap.imageDimension(1));
   size_t nz(m_dataMap.imageDimension(2));
   size_t indx(0);
   for (size_t k = 0; k < nz; k++) {
      double ntot(0);
      for (size_t j = 0; j < ny; j++) {
         for (size_t i = 0; i < nx; i++) {
            ntot += m_dataMap.data().at(indx);
            indx++;
         }
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
      if (k < m_kmin || k > m_kmax-1) {
         continue;
      }
      value += src->pixelCounts(energies().at(k), energies().at(k+1),
                                npreds.at(k), npreds.at(k+1));
   }
   return value;
}

double BinnedLikelihood2::NpredValue(const std::string & srcName,
                                    const SourceMap & sourceMap) const {
   const std::vector<double> & npreds(sourceMap.npreds());
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

void BinnedLikelihood2::getNpreds(const std::string & srcName,
                                 std::vector<double> & npreds) const {
   try {
      npreds = sourceMap(srcName).npreds();
      return;
   } catch (std::runtime_error & eObj) {
      SourceMap * srcMap(getSourceMap(srcName, false));
      npreds = srcMap->npreds();
      delete srcMap;
   }
}

std::vector<double> 
BinnedLikelihood2::countsSpectrum(const std::string & srcName) const {
   std::vector<double> counts_spectrum(m_kmax - m_kmin, 0);
// Compute ratios of individual source model predictions to that of
// the total model.
   double npred;
   if (!m_modelIsCurrent) {
      computeModelMap(npred);
   }
   std::vector<std::pair<double, double> > modelWts;
   std::pair<double, double> zeros(0, 0);
   modelWts.resize(m_filledPixels.size(), zeros);
   addSourceWts(modelWts, srcName);
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t k(m_filledPixels[j]/m_pixels.size());
      if (k < m_kmin || k > m_kmax-1) {
         continue;
      }
      double emin(m_energies[k]);
      double emax(m_energies[k+1]);
      double srcProb = pixelCounts(emin, emax, modelWts[j].first,
                                   modelWts[j].second)/m_model[j];
      size_t indx = m_filledPixels[j];
      counts_spectrum[k - m_kmin] += srcProb*m_dataMap.data()[indx];
   }
   return counts_spectrum;
}

std::vector<double> BinnedLikelihood2::fixedModelSpectrum() const {
   std::vector<double> model_spectrum;
   for (size_t k(0); k < m_energies.size()-1; k++) {
      model_spectrum.push_back(pixelCounts(m_energies[k], m_energies[k+1],
                                           m_fixedNpreds[k],
                                           m_fixedNpreds[k+1]));
   }
   return model_spectrum;
}

void BinnedLikelihood2::create_drm() {
   const astro::SkyDir & ref_dir(m_dataMap.refDir());
   m_drm = new Drm(ref_dir.ra(), ref_dir.dec(), observation(), 
                   m_dataMap.energies());
}

} // namespace Likelihood
