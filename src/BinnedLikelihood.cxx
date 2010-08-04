/**
 * @file BinnedLikelihood.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/BinnedLikelihood.cxx,v 1.65 2010/07/07 01:05:30 jchiang Exp $
 */

#include <cmath>

#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#define ST_DLL_EXPORTS
#include "Likelihood/SourceMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/SourceModel.h"

namespace Likelihood {

BinnedLikelihood::BinnedLikelihood(const CountsMap & dataMap,
                                   const Observation & observation,
                                   const std::string & srcMapsFile,
                                   bool computePointSources,
                                   bool applyPsfCorrections,
                                   bool performConvolution, 
                                   bool resample,
                                   double resamp_factor)
   : LogLike(observation), m_dataMap(dataMap), m_pixels(dataMap.pixels()),
     m_modelIsCurrent(false), m_srcMapsFile(srcMapsFile),
     m_computePointSources(computePointSources),
     m_applyPsfCorrections(applyPsfCorrections),
     m_performConvolution(performConvolution), 
     m_resample(resample), 
     m_resamp_factor(resamp_factor) {
   dataMap.getAxisVector(2, m_energies);
   identifyFilledPixels();
   computeCountsSpectrum();
}

BinnedLikelihood::~BinnedLikelihood() throw() {
   std::map<std::string, SourceMap *>::iterator srcMap(m_srcMaps.begin());
   for ( ; srcMap != m_srcMaps.end(); ++srcMap) {
      delete srcMap->second;
   }
}

double BinnedLikelihood::value(optimizers::Arg & dummy) const {
   (void)(dummy);

   double npred;
   computeModelMap(npred);

   const std::vector<float> & data(m_dataMap.data());
   double my_value(0);

   for (size_t i(0); i < m_filledPixels.size(); i++) {
      if (m_model.at(i) > 0) {
         size_t j(m_filledPixels.at(i));
         my_value += data.at(j)*std::log(m_model.at(i));
      }
   }
   my_value -= npred;

   st_stream::StreamFormatter formatter("BinnedLikelihood", "value", 4);
   formatter.info() << m_nevals << "  "
                    << my_value << "  "
                    << npred << std::endl;
   m_nevals++;
   
   saveBestFit(my_value);
   return my_value;
}

double BinnedLikelihood::npred() {
   double npred;
   computeModelMap(npred);
   return npred;
}

std::vector<double>::const_iterator 
BinnedLikelihood::setParamValues_(std::vector<double>::const_iterator it) {
   m_modelIsCurrent = false;
   return SourceModel::setParamValues_(it);
}

std::vector<double>::const_iterator 
BinnedLikelihood::setFreeParamValues_(std::vector<double>::const_iterator it) {
   m_modelIsCurrent = false;
   return SourceModel::setFreeParamValues_(it);
}

void BinnedLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
   int nparams(getNumFreeParams());
   derivs.resize(nparams, 0);
   double npred;
   if (!m_modelIsCurrent) {
      computeModelMap(npred);
   }
   const std::vector<float> & data = m_dataMap.data();
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t jmin(m_filledPixels.at(j));
      size_t jmax(jmin + m_pixels.size());
      size_t k(jmin/m_pixels.size());
      double emin(m_energies.at(k));
      double emax(m_energies.at(k+1));
      if (m_model.at(j) > 0) {
         long iparam(0);
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
            for (size_t i(0); i < paramNames.size(); i++) {
               double my_deriv = 
                  src->second->pixelCountsDeriv(emin, emax,
                                                model.at(jmin), 
                                                model.at(jmax),
                                                paramNames.at(i));
               derivs.at(iparam) += data.at(jmin)/m_model.at(j)*my_deriv;
               if (j == 0) {
                  const std::vector<double> & npreds = srcMap.npreds();
                  for (size_t kk(0); kk < m_energies.size()-1; kk++) {
                     derivs.at(iparam) -= 
                        src->second->pixelCountsDeriv(m_energies.at(kk), 
                                                      m_energies.at(kk+1),
                                                      npreds.at(kk), 
                                                      npreds.at(kk+1),
                                                      paramNames.at(i));
                  }
               }
               iparam++;
            }
         }
      }
   }
}

CountsMap * BinnedLikelihood::createCountsMap() const {
   std::vector<float> map;
   computeModelMap(map);

   CountsMap * modelMap = new CountsMap(m_dataMap);
         
   modelMap->setImage(map);
   return modelMap;
}

void BinnedLikelihood::readXml(std::string xmlFile, 
                               optimizers::FunctionFactory & funcFactory,
                               bool requireExposure, 
                               bool addPointSources, 
                               bool loadMaps) {
   SourceModel::readXml(xmlFile, funcFactory, requireExposure=false,
                        addPointSources, loadMaps);
   if (m_srcMapsFile == "") {
      createSourceMaps();
   } else {
      buildFixedModelWts();
   }
}

void BinnedLikelihood::addSource(Source * src) {
   m_bestValueSoFar = -1e38;
   SourceModel::addSource(src);
   if (m_srcMaps.find(src->getName()) == m_srcMaps.end()) {
      SourceMap * srcMap(getSourceMap(src->getName()));
      if (srcMap) {
         m_srcMaps[src->getName()] = srcMap;
      }
   }
   std::vector<double> pars;
   src->spectrum().getParamValues(pars);
   m_modelPars[src->getName()] = pars;
}

Source * BinnedLikelihood::deleteSource(const std::string & srcName) {
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

void BinnedLikelihood::computeModelMap(std::vector<float> & modelMap) const {
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
               = const_cast<BinnedLikelihood *>(this)->getSource(name);
            const std::vector<float> & model(srcMap->model());
            modelMap.at(imin) += src->pixelCounts(emin, emax,
                                                  model.at(imin), 
                                                  model.at(imax));
         }
      }
   }
}

bool BinnedLikelihood::fixedModelUpdated() const {
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
            throw std::runtime_error("BinnedLikelihood::fixedModelUpdated: "
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

void BinnedLikelihood::buildFixedModelWts() {
   m_fixedSources.clear();
   m_fixedModelWts.clear();
   m_fixedModelWts.resize(m_filledPixels.size(), std::make_pair(0, 0));
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
         if (srcMapIt != m_srcMaps.end()) {
            m_srcMaps.erase(srcName);
         }
         delete srcMap;
         std::vector<double> pars;
         srcIt->second->spectrum().getParamValues(pars);
         m_modelPars[srcName] = pars;
      } else { 
// Process non-fixed sources.
//
// Ensure model map is avaiable.
         std::map<std::string, SourceMap *>::const_iterator srcMapIt
            = m_srcMaps.find(srcName);
         if (srcMapIt == m_srcMaps.end()) {
            SourceMap * srcMap(getSourceMap(srcName));
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

void BinnedLikelihood::computeModelMap(double & npred) const {
   npred = 0;

   std::vector<std::pair<double, double> > modelWts;
   modelWts.resize(m_filledPixels.size());

   if (fixedModelUpdated()) {
      const_cast<BinnedLikelihood *>(this)->buildFixedModelWts();
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
         addSourceWts(modelWts, srcNames.at(i));
      }
      npred += NpredValue(srcNames.at(i));
   }

   m_model.clear();
   m_model.resize(m_filledPixels.size(), 0);
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t k(m_filledPixels.at(j)/m_pixels.size());
      double emin(m_energies.at(k));
      double emax(m_energies.at(k+1));
      m_model.at(j) = pixelCounts(emin, emax, modelWts.at(j).first,
                                  modelWts.at(j).second);
   }

   m_modelIsCurrent = true;
}

void BinnedLikelihood::
addSourceWts(std::vector<std::pair<double, double> > & modelWts,
             const std::string & srcName,
             const SourceMap * srcMap,
             bool subtract) const {
   const Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
   if (src == 0) {
      return;
   }
   if (modelWts.size() != m_filledPixels.size()) {
      throw std::runtime_error("BinnedLikelihood::addSourceWts: "
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
         throw std::runtime_error("BinnedLikelihood::addSourceWts: "
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

double BinnedLikelihood::spectrum(const Source * src, double energy) const {
   const optimizers::Function & spectrum(src->spectrum());
   optimizers::dArg eArg(energy);
   return spectrum(eArg);
}

double BinnedLikelihood::pixelCounts(double emin, double emax,
                                     double y1, double y2) const {
   if (::getenv("USE_OLD_PIX_EST")) {
      return (y1 + y2)*(emax - emin)/2.;
   }
   double gam(std::log(y2/y1)/std::log(emax/emin));
   if (gam == -1) {
      double y0(y2/std::pow(emax, gam));
      return y0*std::log(emax/emin);
   }

   return y2/(gam + 1.)*(emax - emin*std::pow(emin/emax, gam));
}

void BinnedLikelihood::createSourceMaps() {
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

SourceMap * BinnedLikelihood::createSourceMap(const std::string & srcName) {
   Source * src = getSource(srcName);
   return new SourceMap(src, &m_dataMap, m_observation, m_applyPsfCorrections,
                        m_performConvolution, m_resample, m_resamp_factor);
}

SourceMap * BinnedLikelihood::getSourceMap(const std::string & srcName, 
                                           bool verbose) const {
   if (fileHasSourceMap(srcName, m_srcMapsFile)) {
      return new SourceMap(m_srcMapsFile, srcName);
   }
// Generate the map if it is not in the file.   
   Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
   if (src->getType() == "Diffuse" || m_computePointSources) {
      return new SourceMap(src, &m_dataMap, m_observation,
                           m_applyPsfCorrections,
                           m_performConvolution,
                           m_resample, m_resamp_factor,
                           verbose);
   }
   return 0;
}

void BinnedLikelihood::saveSourceMaps(const std::string & filename) {
   if (filename != "") {
      m_srcMapsFile = filename;
   }
   if (!st_facilities::Util::fileExists(filename)) {
      m_dataMap.writeOutput("BinnedLikelihood", m_srcMapsFile);
   }
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      if (m_srcMaps.count(srcNames.at(i))) {
         if (fileHasSourceMap(srcNames.at(i), m_srcMapsFile)) {
//             replaceSourceMap(srcNames.at(i), m_srcMapsFile);
         } else {
            appendSourceMap(srcNames.at(i), m_srcMapsFile);
         }
      }
   }
}

bool BinnedLikelihood::fileHasSourceMap(const std::string & srcName,
                                        const std::string & fitsFile) const {
   try {
      std::auto_ptr<const tip::Image> 
         image(tip::IFileSvc::instance().readImage(fitsFile, srcName));
   } catch (tip::TipException &) {
      return false;
   }
   return true;
}

void BinnedLikelihood::replaceSourceMap(const std::string & srcName,
                                        const std::string & fitsFile) const {
   const SourceMap & srcMap = sourceMap(srcName);

   std::auto_ptr<tip::Image> 
      image(tip::IFileSvc::instance().editImage(fitsFile, srcName));

   image->set(srcMap.model());
}

void BinnedLikelihood::appendSourceMap(const std::string & srcName,
                                       const std::string & fitsFile) const {
   if (!m_srcMaps.count(srcName)) {
      throw std::runtime_error("BinnedLikelihood::appendSourceMap: " +
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

void BinnedLikelihood::setImageDimensions(tip::Image * image, 
                                          long * dimensions) const {
   typedef std::vector<tip::PixOrd_t> DimCont_t;
   DimCont_t dims = image->getImageDimensions();
   DimCont_t::size_type num_dims = dims.size();
   if (3 != num_dims) {
      throw std::runtime_error("BinnedLikelihood::setImageDimensions: "
                               + std::string("Cannot write a SourceMap file ")
                               + "to an image which is not 3D");
   }
   for (DimCont_t::size_type index = 0; index != num_dims; ++index) {
      dims[index] = dimensions[index];
   }
   image->setImageDimensions(dims);
}

void BinnedLikelihood::identifyFilledPixels() {
   const std::vector<float> & data = m_dataMap.data();
   m_filledPixels.clear();
   for (unsigned int i = 0; i < data.size(); i++) {
      if (data.at(i) > 0) {
         m_filledPixels.push_back(i);
      }
   }
}

void BinnedLikelihood::computeCountsSpectrum() {
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

void BinnedLikelihood::syncParams() {
   SourceModel::syncParams();
}

void BinnedLikelihood::syncSrcParams(const std::string & srcName) {
   (void)(srcName);
   syncParams();
   m_modelIsCurrent = false;
}

double BinnedLikelihood::NpredValue(const std::string & srcName) const {
   std::map<std::string, double>::const_iterator npredIt 
      = m_fixedModelNpreds.find(srcName);
   if (npredIt != m_fixedModelNpreds.end()) {
      return npredIt->second;
   }
   const std::vector<double> & npreds(sourceMap(srcName).npreds());
   const Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
   double value(0);
   for (size_t k(0); k < energies().size()-1; k++) {
      value += src->pixelCounts(energies().at(k), energies().at(k+1),
                                npreds.at(k), npreds.at(k+1));
   }
   return value;
}

double BinnedLikelihood::NpredValue(const std::string & srcName,
                                    const SourceMap & sourceMap) const {
   const std::vector<double> & npreds(sourceMap.npreds());
   const Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
   double value(0);
   for (size_t k(0); k < energies().size()-1; k++) {
      value += src->pixelCounts(energies().at(k), energies().at(k+1),
                                npreds.at(k), npreds.at(k+1));
   }
   return value;
}

void BinnedLikelihood::getNpreds(const std::string & srcName,
                                 std::vector<double> & npreds) const {
   try {
      npreds = sourceMap(srcName).npreds();
      return;
   } catch (std::runtime_error & eObj) {
      SourceMap * srcMap(getSourceMap(srcName));
      npreds = srcMap->npreds();
      delete srcMap;
   }
}

} // namespace Likelihood
