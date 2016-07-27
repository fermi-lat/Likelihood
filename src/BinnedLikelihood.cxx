/**
 * @file BinnedLikelihood.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/BinnedLikelihood.cxx,v 1.119 2016/07/12 02:20:44 mdwood Exp $
 */

#include <cmath>
#include <cstdlib>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "st_stream/StreamFormatter.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"

#include "st_facilities/Util.h"
#include "st_facilities/Timer.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"

#include "Likelihood/Drm.h"
#define ST_DLL_EXPORTS
#include "Likelihood/SourceMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/SourceModel.h"

namespace Likelihood {

BinnedLikelihood::BinnedLikelihood(CountsMapBase & dataMap,
                                   const Observation & observation,
                                   const std::string & srcMapsFile,
                                   bool computePointSources,
                                   bool applyPsfCorrections,
                                   bool performConvolution, 
                                   bool resample,
                                   double resamp_factor,
                                   double minbinsz)
   : LogLike(observation), 
     m_dataMap(dataMap), 
     m_weightMap(0),
     m_weightSrcMap(0),
     m_pixels(dataMap.pixels()),
     m_energies(),
     m_countsSpectrum(),
     m_filledPixels(),
     m_model(),
     m_modelIsCurrent(false),
     m_srcMapsFile(srcMapsFile),
     m_computePointSources(computePointSources),
     m_applyPsfCorrections(applyPsfCorrections),
     m_performConvolution(performConvolution), 
     m_resample(resample), 
     m_resamp_factor(resamp_factor), 
     m_minbinsz(minbinsz),
     m_verbose(true),
     m_updateFixedWeights(true),
     m_fixedSources(),
     m_fixedModelWts(),
     m_fixedModelNpreds(),
     m_fixedModelWeightedNpreds(),
     m_modelPars(),
     m_posDerivs(),
     m_negDerivs(),
     m_kmin(0), m_kmax(0), 
     m_use_edisp(false),
     m_drm(0),
     m_use_single_fixed_map(true),
     m_external_model_map(0),
     m_true_counts(),
     m_meas_counts(),
     m_fixed_counts_spec(),
     m_krefs() {
   dataMap.getEnergies(m_energies);
   m_kmax = m_energies.size() - 1;
   identifyFilledPixels();
   m_fixedModelWts.resize(m_filledPixels.size(), std::make_pair(0, 0));
   m_fixedNpreds.resize(m_energies.size(), 0);

   computeCountsSpectrum();
}

BinnedLikelihood::BinnedLikelihood(CountsMapBase & dataMap,
				   const ProjMap& weightMap,
                                   const Observation & observation,
                                   const std::string & srcMapsFile,
                                   bool computePointSources,
                                   bool applyPsfCorrections,
                                   bool performConvolution, 
                                   bool resample,
                                   double resamp_factor,
                                   double minbinsz)
   : LogLike(observation), 
     m_dataMap(dataMap),
     m_weightMap(&weightMap),
     m_weightSrcMap(0),
     m_pixels(dataMap.pixels()),
     m_energies(),
     m_countsSpectrum(),
     m_filledPixels(),
     m_model(),
     m_modelIsCurrent(false),
     m_srcMapsFile(srcMapsFile),
     m_computePointSources(computePointSources),
     m_applyPsfCorrections(applyPsfCorrections),
     m_performConvolution(performConvolution), 
     m_resample(resample), 
     m_resamp_factor(resamp_factor), 
     m_minbinsz(minbinsz),
     m_verbose(true),
     m_updateFixedWeights(true),
     m_fixedSources(),
     m_fixedModelWts(),
     m_fixedModelNpreds(),
     m_fixedModelWeightedNpreds(),
     m_modelPars(),
     m_posDerivs(),
     m_negDerivs(),
     m_kmin(0), m_kmax(0), 
     m_use_edisp(false),
     m_drm(0),
     m_use_single_fixed_map(true),
     m_external_model_map(0),
     m_true_counts(),
     m_meas_counts(),
     m_fixed_counts_spec(),
     m_krefs() {
   dataMap.getEnergies(m_energies);
   m_kmax = m_energies.size() - 1;

   if ( fileHasSourceMap("__weights__",m_srcMapsFile) ) {
     st_stream::StreamFormatter formatter("BinnedLikelihood","", 2);
     formatter.warn() << "Reading existing weights map from file " << m_srcMapsFile << std::endl;
     m_weightSrcMap = new SourceMap(m_srcMapsFile, "__weights__", m_observation,0,true);
   } else {
     m_weightSrcMap = new SourceMap(weightMap,&dataMap,observation,true);
   }
   identifyFilledPixels();
   m_fixedModelWts.resize(m_filledPixels.size(), std::make_pair(0, 0));
   m_fixedNpreds.resize(m_energies.size(), 0);
   computeCountsSpectrum();
}

BinnedLikelihood::
BinnedLikelihood(const BinnedLikelihood & other) 
   : LogLike(other),
     m_dataMap(other.m_dataMap),
     m_weightMap(other.m_weightMap),
     m_weightSrcMap(other.m_weightSrcMap ? new SourceMap(*m_weightMap,&m_dataMap,observation(),true) : 0),
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
     m_resamp_factor(other.m_resamp_factor),
     m_minbinsz(other.m_minbinsz),
     m_verbose(other.m_verbose),
     m_updateFixedWeights(other.m_updateFixedWeights),
     m_fixedSources(other.m_fixedSources),
     m_fixedModelWts(other.m_fixedModelWts),
     m_fixedModelNpreds(other.m_fixedModelNpreds),
     m_fixedModelWeightedNpreds(other.m_fixedModelWeightedNpreds),
     m_modelPars(other.m_modelPars),
     m_posDerivs(other.m_posDerivs),
     m_negDerivs(other.m_negDerivs),
     m_kmin(other.m_kmin), m_kmax(other.m_kmax),
     m_fixedNpreds(other.m_fixedNpreds),
     m_use_edisp(other.m_use_edisp),
     m_drm(new Drm(*other.m_drm)),
     m_use_single_fixed_map(other.m_use_single_fixed_map),
     m_external_model_map(other.m_external_model_map),
     m_true_counts(other.m_true_counts),
     m_meas_counts(other.m_meas_counts),
     m_fixed_counts_spec(other.m_fixed_counts_spec),
     m_krefs(other.m_krefs) {
   std::map<std::string, SourceMap *>::const_iterator 
      it(other.m_srcMaps.begin());
   for ( ; it != other.m_srcMaps.end(); ++it) {
      m_srcMaps[it->first] = new SourceMap(*(it->second));
  }
}

BinnedLikelihood::~BinnedLikelihood() throw() {
   std::map<std::string, SourceMap *>::iterator srcMap(m_srcMaps.begin());
   for ( ; srcMap != m_srcMaps.end(); ++srcMap) {
      delete srcMap->second;
   }
   delete m_weightSrcMap;
   delete m_drm;
}

double BinnedLikelihood::value(optimizers::Arg & dummy) const {
   (void)(dummy);

   // Here we want the weighted verison of the nPred
   double npred(computeModelMap(true));

   const std::vector<float> & data = m_weightSrcMap == 0 ? m_dataMap.data() : m_weightedCounts;
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

   st_stream::StreamFormatter formatter("BinnedLikelihood", "value", 4);
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

double BinnedLikelihood::npred(bool weighted) {
   return computeModelMap(weighted);
}

void BinnedLikelihood::setSourceMapImage(const std::string & name,
					 const std::vector<float>& image) {
  // Create the source map if it doesn't exist
  if (m_srcMaps.find(name) == m_srcMaps.end()) {
    m_srcMaps[name] = getSourceMap(name);
  }
  m_srcMaps[name]->setImage(image);
}

bool BinnedLikelihood::hasSourceMap(const std::string & name) const {
  if (m_srcMaps.find(name) == m_srcMaps.end())
    return false;
  else
    return true;
}

const SourceMap & BinnedLikelihood::sourceMap(const std::string & name) const {
  std::vector<std::string> srcNames;
  getSrcNames(srcNames);

  std::vector<std::string>::iterator nameIt = 
    std::find(srcNames.begin(),srcNames.end(),name);
  if(nameIt == srcNames.end())
    throw std::runtime_error("Cannot find source named: " + name);

  if (m_srcMaps.find(name) == m_srcMaps.end()) {
    m_srcMaps[name] = getSourceMap(name);
  }

  return *m_srcMaps[name];
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
   st_stream::StreamFormatter formatter("BinnedLikelihood", "getFreeDerivs",4);

   int nparams(getNumFreeParams());
   derivs.resize(nparams, 0);
   double npred;
   if (!m_modelIsCurrent) {
     // here we want the weighted version of npred
      npred = computeModelMap(true);
   }
   const std::vector<float> & data = m_weightSrcMap == 0 ? m_dataMap.data() : m_weightedCounts;

/// First j value corresponding to the minimum allowed k-index.
   size_t jentry(0);
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t k(m_filledPixels[j]/m_pixels.size());
      if (k == m_kmin) {
         jentry = j;
         break;
      }
   }

/// Get handles of sources with free parameters and names of those
/// free parameters. The vectors free_srcs and free_params will be
/// traversed in tandem inside the main loop.  Extracting this
/// information ahead of time will save substantially on execution
/// time for analyses with a large number of sources.
   std::vector<Source *> free_srcs;
   std::vector< std::vector<std::string> > free_params;
   std::map<std::string, Source *>::const_iterator src;
   for (src=sources().begin(); src != sources().end(); ++src) {
      if (!std::count(m_fixedSources.begin(), m_fixedSources.end(),
                      src->second->getName())) {
         free_srcs.push_back(src->second);
         std::vector<std::string> parnames;
         src->second->spectrum().getFreeParamNames(parnames);
         free_params.push_back(parnames);
      }
   }

// /// For debugging/development, keep track of time spent in main loop.
//    st_facilities::Timer timer;

   /// Not clear why these should be data members.
   m_posDerivs.clear();
   m_negDerivs.clear();
//   timer.start();
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
	    m_posDerivs[iparam] = Kahan_Accumulator();
	    m_negDerivs[iparam] = Kahan_Accumulator();
         }

         size_t ipar(0);
         for (std::vector<Source *>::const_iterator it(free_srcs.begin());
              it != free_srcs.end(); ++it, ipar++) {
            Source * src(*it);

            std::string srcName = src->getName();
            const SourceMap & srcMap = sourceMap(srcName);

            const std::vector<float> & model = srcMap.model();
            const std::vector<std::string> & paramNames(free_params[ipar]);

            for (size_t i(0); i < paramNames.size(); i++, iparam++) {
               double my_deriv = 
                  src->pixelCountsDeriv(emin, emax,
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
                  const std::vector<double> & npreds =  srcMap.npreds();
                  const std::vector<std::pair<double,double> > & npred_weights =  srcMap.npred_weights();
                  for (size_t kk(0); kk < m_energies.size()-1; kk++) {
                     if (kk >= m_kmin && kk <= m_kmax-1) {
                        addend = 
                           src->pixelCountsDeriv(m_energies.at(kk), 
                                                 m_energies.at(kk+1),
                                                 npreds.at(kk)*npred_weights[kk].first,
                                                 npreds.at(kk+1)*npred_weights[kk].second,
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
//    timer.stop();
//    timer.report("main loop time");

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

CountsMapBase * BinnedLikelihood::createCountsMap() const {
   std::vector<float> map;
   computeModelMap(map);

   CountsMapBase * modelMap = m_dataMap.clone();         
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
}

void BinnedLikelihood::addSource(Source * src, bool fromClone) {
   m_bestValueSoFar = -1e38;
   SourceModel::addSource(src, fromClone);
   if (use_single_fixed_map() && src->fixedSpectrum()) {
      addFixedSource(src->getName());
   } else {
      SourceMap * srcMap(getSourceMap(src->getName(), true));
      if (srcMap) {
         m_srcMaps[src->getName()] = srcMap;
         if (m_external_model_map) {
            // This should only execute from within gtmodel, in 
            // which case the SourceMap is no longer needed and thus
            // be deleted.
            updateModelMap(*m_external_model_map, srcMap);
            delete srcMap;
            srcMap = 0;
            m_srcMaps.erase(src->getName());
         }
      }
      std::vector<double> pars;
      src->spectrum().getParamValues(pars);
      m_modelPars[src->getName()] = pars;
   }
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
      m_fixedModelWeightedNpreds.erase(srcName);
      m_fixedSources.erase(srcIt);
   }
   Source * src(SourceModel::deleteSource(srcName));
   return src;
}

void BinnedLikelihood::eraseSourceMap(const std::string & srcName) {
   delete m_srcMaps[srcName];
   m_srcMaps.erase(srcName);
}

void BinnedLikelihood::
updateModelMap(std::vector<float> & modelMap,
               const SourceMap * srcMap) const {
   size_t npix(m_pixels.size());
   std::string name(srcMap->name());
   NpredValue(name); // This computes the convolved spectrum.
   const Source * src = const_cast<BinnedLikelihood *>(this)->getSource(name);
   const std::vector<float> & model(srcMap->model());

   const SourceMap* mask = srcMap->weights();
   
   for (size_t j(0); j < npix; j++) {
      for (size_t k(0); k < m_energies.size()-1; k++) {
         double emin(m_energies.at(k));
         double emax(m_energies.at(k+1));
         size_t jmin(k*npix + j);
         size_t jmax(jmin + npix);
	 // EAC, skip masked pixels
	 if ( mask && 
	      ( mask->model().at(jmin) <= 0. ) ) continue;
         double wt1(0);
         double wt2(0);
         if (use_edisp(name)) {
            if (m_true_counts[name][k] != 0) {
               double xi(m_meas_counts[name][k]
                         /m_true_counts[name][k]);
               wt1 = model[jmin]*xi;
               wt2 = model[jmax]*xi;
            } else {
               size_t ipix(jmin % npix);
               std::map<std::string, std::map<size_t, size_t> >::
                  const_iterator kref_it(m_krefs.find(name));
               std::map<size_t, size_t>::const_iterator it =
                  kref_it->second.find(ipix);
               if (it != kref_it->second.end()) {
                  size_t kref = it->second;
                  size_t jref = kref*npix + ipix;
                  wt1 = (model[jref]/m_true_counts[name][kref]
                         *m_meas_counts[name][k]);
                  wt2 = (model[jref+npix]/m_true_counts[name][kref]
                         *m_meas_counts[name][k]);
               }
            }
         } else {
            wt1 = model[jmin];
            wt2 = model[jmax];
         }
         modelMap[jmin] += src->pixelCounts(emin, emax, wt1, wt2);
      }
   }
}

void BinnedLikelihood::
saveWeightsMap(bool replace) const {
  if ( m_weightSrcMap == 0 ) return;
  bool has_weights = fileHasSourceMap(m_weightSrcMap->name(), m_srcMapsFile);
  if ( has_weights ) {
    if ( replace ) {
      switch ( m_dataMap.projection().method() ) {
      case astro::ProjBase::WCS:
	replaceSourceMap_wcs(*m_weightSrcMap,m_srcMapsFile);
	return;
      case astro::ProjBase::HEALPIX:
	replaceSourceMap_healpix(*m_weightSrcMap,m_srcMapsFile);
	return;
      default:
	break;
      }
    } else {
      // just leave it be;
      ;
    }
  } else {
    switch ( m_dataMap.projection().method() ) {
    case astro::ProjBase::WCS:
      appendSourceMap_wcs(*m_weightSrcMap,m_srcMapsFile,true);
      return;
    case astro::ProjBase::HEALPIX:
      appendSourceMap_healpix(*m_weightSrcMap,m_srcMapsFile,true);
      return;
    default:
      break;
    }
  }
}



void BinnedLikelihood::
fillWeightedCounts() {
  m_weightedCounts.clear();
  size_t ne = m_energies.size();
  size_t npix = m_pixels.size();
  m_weightedCounts.resize(ne*npix);
  for ( size_t j(0); j < npix; j++ ) {
    for (size_t k(0); k < ne-1; k++) {
      size_t idx = k*npix +j;
      double w = m_weightSrcMap->model()[idx];
      bool is_null = w <= 0 || m_dataMap.data()[idx] <= 0;
      w = is_null ? 0 : w;
      m_weightedCounts[idx] = m_dataMap.data()[idx] * w;
    }
  }
}

void BinnedLikelihood::
computeModelMap(std::vector<float> & modelMap) const {
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   computeModelMap(srcNames,modelMap);
}

void BinnedLikelihood::
computeModelMap(const std::string & srcName, std::vector<float> & modelMap) const {
   size_t npix(m_pixels.size());
   modelMap.clear();
   modelMap.resize(npix*(m_energies.size()-1), 0);
   bool hasMap = hasSourceMap(srcName);
   const SourceMap& srcMap = sourceMap(srcName);
   updateCorrectionFactors(srcName,srcMap);
   updateModelMap(modelMap, &srcMap);
   if( !hasMap ) {
     eraseSourceMap(srcName);
   }
}

void BinnedLikelihood::
computeModelMap(const std::vector<std::string>& srcNames, std::vector<float> & modelMap) const {
   size_t npix(m_pixels.size());
   modelMap.clear();
   modelMap.resize(npix*(m_energies.size()-1), 0);

   std::vector<std::string>::const_iterator name = srcNames.begin();
   for ( ; name != srcNames.end(); ++name) {  
     bool hasMap = hasSourceMap(*name);
     const SourceMap& srcMap = sourceMap(*name);
     updateCorrectionFactors(*name,srcMap);
     updateModelMap(modelMap, &srcMap);
     if( !hasMap ) {
       eraseSourceMap(*name);
     }
   }
}

void BinnedLikelihood::setCountsMap(const std::vector<float> & counts) {

  if(counts.size() != m_dataMap.data().size())
    throw std::runtime_error("Wrong size for input counts map.");
  m_dataMap.setImage(counts);
  identifyFilledPixels();
  buildFixedModelWts();
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

const std::vector<double> & BinnedLikelihood::
modelCountsSpectrum(const std::string & srcname) const {
   if (!m_modelIsCurrent) {
      computeModelMap();
   }
   std::map<std::string, std::vector<double> >::const_iterator it;
   it = m_meas_counts.find(srcname);
   if (it == m_meas_counts.end()) {
      try {
         NpredValue(srcname);
         it = m_meas_counts.find(srcname);
      } catch(...) {
         throw std::runtime_error("Cannot find model counts spectrum for "
                                  + srcname);
      }
   }
   return it->second;
}

void BinnedLikelihood::addFixedSource(const std::string & srcName) {
// Add a source to the fixed source data, under the assumption that it
// is not already there.
   std::map<std::string, Source *>::const_iterator 
      srcIt(m_sources.find(srcName));
   if (srcIt == m_sources.end()) {
      std::ostringstream message;
      message << "BinnedLikelihood::addFixedSource: "
              << "source " << srcName << " not found.";
      throw std::runtime_error(message.str());
   }

   if (std::count(m_fixedSources.begin(), m_fixedSources.end(), srcName) != 0) {
      std::ostringstream message;
      message << "BinnedLikelihood::addFixedSource: "
              << "source " << srcName << " already in fixed model.";
      throw std::runtime_error(message.str());
   }

   m_fixedSources.push_back(srcName);

   SourceMap * srcMap(0);
   std::map<std::string, SourceMap *>::const_iterator srcMapIt
      = m_srcMaps.find(srcName);
   if (srcMapIt == m_srcMaps.end()) {
      srcMap = getSourceMap(srcName, true);
   } else {
      srcMap = srcMapIt->second;
   }
   // Store the _weighted_ npreds
   m_fixedModelNpreds[srcName] = NpredValue(srcName, *srcMap);
   m_fixedModelWeightedNpreds[srcName] = NpredValue(srcName, *srcMap, true);

   addSourceWts(m_fixedModelWts, srcName, srcMap);
   double xi(1);
   for (size_t k(0); k < m_energies.size()-1; k++) {
      optimizers::dArg ee(m_energies[k]);
      if (use_edisp(srcName)) {
         if (m_true_counts[srcName][k] != 0) {
            xi = m_meas_counts[srcName][k]/m_true_counts[srcName][k];
         } else {
            xi = 1;
         }
      }
      m_fixedNpreds[k] += xi*srcIt->second->spectrum()(ee)*srcMap->npreds()[k];
   }
   // For last bin edge, use xi value from penultimate edge
   size_t kk(m_energies.size()-1);
   optimizers::dArg ee(m_energies[kk]);
   m_fixedNpreds[kk] += xi*(srcIt->second->spectrum()(ee)*
                            srcMap->npreds()[kk]);

   // Remove this source from the stored source maps to save memory
   if (srcMapIt != m_srcMaps.end()) {
      m_srcMaps.erase(srcName);
   }
   delete srcMap;
   srcMap = 0;
   std::vector<double> pars;
   srcIt->second->spectrum().getParamValues(pars);
   m_modelPars[srcName] = pars;
}

void BinnedLikelihood::deleteFixedSource(const std::string & srcName) {
// Delete a source from the fixed source data, under the assumption 
// that it is included.
   std::map<std::string, Source *>::const_iterator 
      srcIt(m_sources.find(srcName));
   if (srcIt == m_sources.end()) {
      std::ostringstream message;
      message << "BinnedLikelihood::addFixedSource: "
              << "source " << srcName << " not found.";
      throw std::runtime_error(message.str());
   }

   // Generate the SourceMap and include it in the stored maps.
   SourceMap * srcMap(getSourceMap(srcName, false));
   if (srcMap == 0) {
      throw std::runtime_error("SourceMap cannot be created for " + srcName);
   }
   m_srcMaps[srcName] = srcMap;

   // Subtract the contribution to the summed Npred spectrum.
   double xi(1);
   for (size_t k(0); k < m_energies.size(); k++) {
      optimizers::dArg ee(m_energies[k]);
      if (use_edisp(srcName)) {
         if (m_true_counts[srcName][k] != 0) {
            xi = m_meas_counts[srcName][k]/m_true_counts[srcName][k];
         } else {
            xi = 1;
         }
      }
      m_fixedNpreds[k] -= xi*srcIt->second->spectrum()(ee)*srcMap->npreds()[k];
   }
   // For last bin edge, use xi value from penultimate edge
   size_t kk(m_energies.size()-1);
   optimizers::dArg ee(m_energies[kk]);
   m_fixedNpreds[kk] -= xi*(srcIt->second->spectrum()(ee)*
                            srcMap->npreds()[kk]);

   // Subtract the source weights from the summed fixed model weights.
   bool subtract;
   addSourceWts(m_fixedModelWts, srcName, srcMap, subtract=true);

   // Remove from the list of fixed sources
   std::vector<std::string>::iterator it 
      = std::find(m_fixedSources.begin(), m_fixedSources.end(), srcName);
   m_fixedSources.erase(it);
}

void BinnedLikelihood::setUpdateFixedWeights(bool update) {
  m_updateFixedWeights = update;
}

void BinnedLikelihood::buildFixedModelWts(bool process_all) {
   m_fixedSources.clear();
   m_fixedModelWts.clear();
   m_fixedModelWts.resize(m_filledPixels.size(), std::make_pair(0, 0));
   m_fixedNpreds.clear();
   m_fixedNpreds.resize(m_energies.size(), 0);
   std::map<std::string, Source *>::const_iterator srcIt(m_sources.begin());
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      const std::string & srcName(srcIt->first);
      if (srcIt->second->fixedSpectrum() && !process_all) {
         addFixedSource(srcName);
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
         if (m_fixedModelWeightedNpreds.find(srcName) != m_fixedModelWeightedNpreds.end()) {
            m_fixedModelWeightedNpreds.erase(srcName);
         }	 
      }
   }

   computeFixedCountsSpectrum();
}

double BinnedLikelihood::computeModelMap(bool weighted) const {
   double npred(0);

   std::vector<std::pair<double, double> > modelWts;
   modelWts.resize(m_filledPixels.size());

   if (fixedModelUpdated() && m_updateFixedWeights) {
      const_cast<BinnedLikelihood *>(this)->buildFixedModelWts();
   }

   for (size_t j(0); j < m_fixedModelWts.size(); j++) {
      modelWts.at(j).first = m_fixedModelWts.at(j).first;
      modelWts.at(j).second = m_fixedModelWts.at(j).second;
   }

   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (size_t i(0); i < srcNames.size(); i++) {
      npred += NpredValue(srcNames[i],weighted);
      if (std::count(m_fixedSources.begin(), m_fixedSources.end(),
                     srcNames.at(i)) == 0) {
         addSourceWts(modelWts, srcNames[i]);
      }
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

   return npred;
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

   std::vector<double> spec;
   for (size_t k(0); k < m_energies.size(); k++) {
      spec.push_back(spectrum(src, m_energies[k]));
   }

   size_t npix(m_pixels.size());
   const std::vector<float> & model(srcMap->model());
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t jmin(m_filledPixels.at(j));
      size_t jmax(jmin + npix);
      size_t k(jmin/npix);
      if (use_edisp(srcName)) {
         if (m_true_counts[srcName][k] != 0) {
            double xi(m_meas_counts[srcName][k]/m_true_counts[srcName][k]);
            modelWts[j].first += my_sign*model[jmin]*spec[k]*xi;
            modelWts[j].second += my_sign*model[jmax]*spec[k+1]*xi;
         } else {
            size_t ipix(jmin % npix);
            std::map<std::string, std::map<size_t, size_t> >::
               const_iterator kref_it(m_krefs.find(srcName));
            std::map<size_t, size_t>::const_iterator it =
               kref_it->second.find(ipix);
            if (it != kref_it->second.end()) {
               size_t kref = it->second;
               size_t jref = kref*npix + ipix;
               modelWts[j].first += (my_sign*model[jref]*spec[kref]
                                     /m_true_counts[srcName][kref]
                                     *m_meas_counts[srcName][k]);
               modelWts[j].second += (my_sign*model[jref+npix]*spec[kref+1]
                                      /m_true_counts[srcName][kref]
                                      *m_meas_counts[srcName][k]);
            }
         }
      } else {
         modelWts[j].first += my_sign*model[jmin]*spec[k];
         modelWts[j].second += my_sign*model[jmax]*spec[k+1];
      }
   }
}

double BinnedLikelihood::spectrum(const Source * src, double energy) const {
   const optimizers::Function & spectrum(src->spectrum());
   optimizers::dArg eArg(energy);
   return spectrum(eArg);
}

double BinnedLikelihood::pixelCounts(double emin, double emax,
                                     double y1, double y2) const {
   if (::getenv("USE_LINEAR_QUADRATURE")) {
      return (y1 + y2)*(emax - emin)/2.;
   }
   return (y1*emin + y2*emax)/2.*std::log(emax/emin);
}

void BinnedLikelihood::loadSourceMaps(bool recreate, bool saveMaps) {
  std::vector<std::string> srcNames;
  getSrcNames(srcNames);
  loadSourceMaps(srcNames,recreate,saveMaps);
}

void BinnedLikelihood::loadSourceMaps(const  std::vector<std::string>& srcNames,
				      bool recreate, bool saveMaps) {

   std::vector<std::string>::const_iterator name = srcNames.begin();
   for ( ; name != srcNames.end(); ++name) {     

     if (m_srcMaps.find(*name) == m_srcMaps.end() || recreate) 
       loadSourceMap(*name,recreate,false);

     if(saveMaps) {
       if(fileHasSourceMap(*name, m_srcMapsFile)) {
	 replaceSourceMap(*name, m_srcMapsFile);
       } else if(!fileHasSourceMap(*name, m_srcMapsFile)) {
	 appendSourceMap(*name, m_srcMapsFile);
       }
     }
   }

   if(recreate) {
     buildFixedModelWts();
   }
}

void BinnedLikelihood::loadSourceMap(const std::string & srcName, bool recreate,
				     bool buildFixedWeights) {
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);

   std::vector<std::string>::iterator nameIt = 
     std::find(srcNames.begin(),srcNames.end(),srcName);
   if(nameIt == srcNames.end())
     throw std::runtime_error("Cannot find source named: " + srcName);
     
   Source * src = getSource(srcName);
 
   if(!(src->getType() == "Diffuse" || m_computePointSources))
     return;

   std::map<std::string, SourceMap *>::iterator mapIt = m_srcMaps.find(srcName);

   if( mapIt != m_srcMaps.end() ) {
     delete m_srcMaps[srcName];
     m_srcMaps.erase(srcName);
   }

   SourceMap * srcMap = 0;

   if(recreate) {
     srcMap = createSourceMap(srcName);
   } else {
     srcMap = getSourceMap(srcName);
   }

   if(srcMap) {
     m_srcMaps[srcName] = srcMap;

     std::vector<std::string>::iterator srcIt = 
       std::find(m_fixedSources.begin(), m_fixedSources.end(), srcName);
     if (srcIt != m_fixedSources.end() && buildFixedWeights) {
       buildFixedModelWts();
     }
   }
}

SourceMap * BinnedLikelihood::createSourceMap(const std::string & srcName) {
   Source * src = getSource(srcName);
   return new SourceMap(src, &m_dataMap, m_observation, m_applyPsfCorrections,
                        m_performConvolution, m_resample, m_resamp_factor,
                        m_minbinsz, m_verbose, m_weightSrcMap);
}

SourceMap * BinnedLikelihood::getSourceMap(const std::string & srcName, 
                                           bool verbose) const {
   if (fileHasSourceMap(srcName, m_srcMapsFile)) {
      return new SourceMap(m_srcMapsFile, srcName, m_observation, m_weightSrcMap);
   }
// Generate the map if it is not in the file.   
   Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
  
   if (src->getType() == "Diffuse" || m_computePointSources) {
      return new SourceMap(src, &m_dataMap, m_observation,
                           m_applyPsfCorrections,
                           m_performConvolution,
                           m_resample, m_resamp_factor,
                           m_minbinsz, verbose && m_verbose,
			   m_weightSrcMap);
   }
   return 0;
}

void BinnedLikelihood::saveSourceMaps(const std::string & filename,
				      bool replace) {
   if (filename != "") {
      m_srcMapsFile = filename;
   }
   if (!st_facilities::Util::fileExists(filename)) {
      m_dataMap.writeOutput("BinnedLikelihood", m_srcMapsFile);
   }
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      st_stream::StreamFormatter formatter("BinnedLikelihood",
                                           "saveSourceMaps", 4);
      formatter.info() << srcNames.at(i) << std::endl;
      if (m_srcMaps.count(srcNames.at(i))) {
         if (fileHasSourceMap(srcNames.at(i), m_srcMapsFile) && replace) {
	   replaceSourceMap(srcNames.at(i), m_srcMapsFile);
         } else if(!fileHasSourceMap(srcNames.at(i), m_srcMapsFile)) {
            formatter.info() << "appending map for " 
                             << srcNames.at(i) << std::endl;
            appendSourceMap(srcNames.at(i), m_srcMapsFile);
         }
      }
   }
   saveWeightsMap(replace);
}

bool BinnedLikelihood::fileHasSourceMap(const std::string & srcName,
                                        const std::string & fitsFile) const {
   try {
      std::auto_ptr<const tip::Extension> ext(tip::IFileSvc::instance().readExtension(fitsFile,srcName));
   } catch (tip::TipException &) {
      return false;
   }
   return true;
}

void BinnedLikelihood::replaceSourceMap(const std::string & srcName,
                                        const std::string & fitsFile) const {
   const SourceMap & srcMap = sourceMap(srcName);
   switch ( m_dataMap.projection().method() ) {
   case astro::ProjBase::WCS:
     replaceSourceMap_wcs(srcMap,fitsFile);
     return;
   case astro::ProjBase::HEALPIX:
     replaceSourceMap_healpix(srcMap,fitsFile);
     return;
   default:
     break;
   }
   std::string errMsg("BinnedLikelihood did not recognize projection method used for CountsMap: ");
   errMsg += m_dataMap.filename();
   throw std::runtime_error(errMsg);
}

void BinnedLikelihood::replaceSourceMap_wcs(const SourceMap& srcMap, 
					    const std::string & fitsFile) const {
   std::auto_ptr<tip::Image> 
       image(tip::IFileSvc::instance().editImage(fitsFile, srcMap.name()));  
   image->set(srcMap.model());
}

void BinnedLikelihood::replaceSourceMap_healpix(const SourceMap& srcMap, 
						const std::string & fitsFile) const {
   std::auto_ptr<tip::Table> 
       table(tip::IFileSvc::instance().editTable(fitsFile, srcMap.name()));  
   m_dataMap.setKeywords(table->getHeader());
   
   const CountsMapHealpix& dataMap_healpix = static_cast<const CountsMapHealpix&>(m_dataMap);
   if ( !dataMap_healpix.allSky() ) {
     tip::Header & header(table->getHeader());     
     header["INDXSCHM"].set("EXPLICIT");
     header["REFDIR1"].set(dataMap_healpix.isGalactic() ? dataMap_healpix.refDir().l() :  dataMap_healpix.refDir().ra() );
     header["REFDIR2"].set(dataMap_healpix.isGalactic() ? dataMap_healpix.refDir().b() :  dataMap_healpix.refDir().dec() );
     header["MAPSIZE"].set(dataMap_healpix.mapRadius());     
     std::string pixname("PIX");
     table->appendField(pixname, std::string("J"));
     tip::IColumn* col = table->getColumn(table->getFieldIndex(pixname));
     int writeValue(-1);
     for(int iloc(0); iloc < dataMap_healpix.nPixels(); iloc++ ) {
       writeValue = dataMap_healpix.localToGlobalIndex(iloc);
       col->set(iloc,writeValue);
     }
   }

   long nPix = m_dataMap.imageDimension(0);
   long nEBins = m_dataMap.energies().size();
   long idx(0);
   double writeValue(0);
   for (long e_index = 0; e_index != nEBins; e_index++ ) {
      std::ostringstream e_channel;
      e_channel<<"CHANNEL"<<e_index+1;
      // Check to see if the column already exist
      tip::FieldIndex_t col_idx = table->getFieldIndex(e_channel.str());
      if ( col_idx < 0 ) {
	table->appendField(e_channel.str(), std::string("D"));
	col_idx = table->getFieldIndex(e_channel.str());
      }
      tip::IColumn* col = table->getColumn(col_idx);
      for(long hpx_index = 0; hpx_index != nPix; ++hpx_index, idx++) {
	writeValue = double(srcMap.model()[idx]);
	col->set(hpx_index,writeValue);
      }
   }
}

void BinnedLikelihood::appendSourceMap(const std::string & srcName,
                                       const std::string & fitsFile,
				       bool isWeights ) const {
   const SourceMap & srcMap = sourceMap(srcName);
   switch ( m_dataMap.projection().method() ) {
   case astro::ProjBase::WCS:
     appendSourceMap_wcs(srcMap,fitsFile,isWeights);
     return;
   case astro::ProjBase::HEALPIX:
     appendSourceMap_healpix(srcMap,fitsFile,isWeights);
     return;
   default:
     break;
   }
   std::string errMsg("BinnedLikelihood did not recognize projection method used for CountsMap: ");
   errMsg += m_dataMap.filename();
   throw std::runtime_error(errMsg);
   return;
}

void BinnedLikelihood::appendSourceMap_wcs(const SourceMap& srcMap,	
					   const std::string & fitsFile,
					   bool isWeights) const {
   std::vector<long> naxes;
   naxes.push_back(m_dataMap.imageDimension(0));
   naxes.push_back(m_dataMap.imageDimension(1));
   long nEBins = isWeights ? m_dataMap.energies().size()-1 : m_dataMap.energies().size();
   naxes.push_back(nEBins);

   tip::IFileSvc::instance().appendImage(fitsFile, srcMap.name(), naxes);
   tip::Image * image = tip::IFileSvc::instance().editImage(fitsFile,srcMap.name());
      
   image->set(srcMap.model());
   m_dataMap.setKeywords(image->getHeader());
   delete image;
}

void BinnedLikelihood::appendSourceMap_healpix(const SourceMap& srcMap, 
					       const std::string & fitsFile,
					       bool isWeights) const {
   tip::IFileSvc::instance().appendTable(fitsFile, srcMap.name());
   std::auto_ptr<tip::Table> 
       table(tip::IFileSvc::instance().editTable(fitsFile, srcMap.name()));  
   m_dataMap.setKeywords(table->getHeader());
   const CountsMapHealpix& dataMap_healpix = static_cast<const CountsMapHealpix&>(m_dataMap);
   if ( !dataMap_healpix.allSky() ) {
     tip::Header & header(table->getHeader());     
     header["INDXSCHM"].set("EXPLICIT");
     header["REFDIR1"].set(dataMap_healpix.isGalactic() ? dataMap_healpix.refDir().l() :  dataMap_healpix.refDir().ra() );
     header["REFDIR2"].set(dataMap_healpix.isGalactic() ? dataMap_healpix.refDir().b() :  dataMap_healpix.refDir().dec() );
     header["MAPSIZE"].set(dataMap_healpix.mapRadius());     
     std::string pixname("PIX");
     table->appendField(pixname, std::string("J"));
     tip::IColumn* col = table->getColumn(table->getFieldIndex(pixname));
     int writeValue(-1);
     for(int iloc(0); iloc < dataMap_healpix.nPixels(); iloc++ ) {
       writeValue = dataMap_healpix.localToGlobalIndex(iloc);
       col->set(iloc,writeValue);
     }
   }
   long nPix = m_dataMap.imageDimension(0);
   long nEBins = isWeights ? m_dataMap.energies().size() : m_dataMap.energies().size() -1;
   long idx(0);
   double writeValue(0.);
   for (long e_index = 0; e_index != nEBins; e_index++ ) {
      std::ostringstream e_channel;
      e_channel<<"CHANNEL"<<e_index+1;
      table->appendField(e_channel.str(), std::string("D"));
      tip::FieldIndex_t col_idx = table->getFieldIndex(e_channel.str());
      tip::IColumn* col = table->getColumn(col_idx);
      for(long hpx_index = 0; hpx_index != nPix; ++hpx_index, idx++) {
	 writeValue = double(srcMap.model()[idx]);
	 col->set(hpx_index,writeValue);
      }
   }
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
   if ( m_weightSrcMap ) {
      fillWeightedCounts();
   }
   const std::vector<float> & data = m_weightSrcMap == 0 ? m_dataMap.data() : m_weightedCounts;
   m_filledPixels.clear();
   for (unsigned int i = 0; i < data.size(); i++) {
      if (data.at(i) > 0) {
         m_filledPixels.push_back(i);
      }
   }
}

void BinnedLikelihood::computeCountsSpectrum() {
   // EAC_FIX, CountsMap should be able to do this 
   switch ( m_dataMap.projection().method() ) {
   case astro::ProjBase::WCS:
     computeCountsSpectrum_wcs();
     return;
   case astro::ProjBase::HEALPIX:
     computeCountsSpectrum_healpix();
     return;
   default:
     break;
   }
   std::string errMsg("BinnedLikelihood did not recognize projection method used for CountsMap: ");
   errMsg += m_dataMap.filename();
   throw std::runtime_error(errMsg);

}

void BinnedLikelihood::computeCountsSpectrum_wcs() {
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

void BinnedLikelihood::computeCountsSpectrum_healpix() {
   m_countsSpectrum.clear();
   size_t nx(m_dataMap.imageDimension(0));
   size_t ny(m_dataMap.imageDimension(1));
   size_t indx(0);
   for (size_t k = 0; k < ny; k++) {
      double ntot(0);
      for (size_t i = 0; i < nx; i++) {
	ntot += m_dataMap.data().at(indx);
	indx++;
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

double BinnedLikelihood::NpredValue(const std::string & srcName, bool weighted) const {
   if ( weighted ) { 
     std::map<std::string, double>::const_iterator npredIt 
       = m_fixedModelWeightedNpreds.find(srcName);
     if (npredIt != m_fixedModelWeightedNpreds.end()) {
       return npredIt->second;
     }
   } else {
     std::map<std::string, double>::const_iterator npredIt 
       = m_fixedModelNpreds.find(srcName);
     if (npredIt != m_fixedModelNpreds.end()) {
       return npredIt->second;
     }
   }
   return NpredValue(srcName, sourceMap(srcName), weighted);
}

// This version forces the recalculation of Npred, whether the source is
// fixed or not. It is also called from buildFixedModelWts.
double BinnedLikelihood::NpredValue(const std::string & srcName,
                                    const SourceMap & sourceMap, 
				    bool weighted) const {
  updateCorrectionFactors(srcName,sourceMap);

  const std::vector<double>& true_counts_spec = m_true_counts[srcName];
  const std::vector<double>& meas_counts_spec = m_meas_counts[srcName];

   double value(0);
   for (size_t k(0); k < energies().size()-1; k++) {
      if (k < m_kmin || k > m_kmax-1) {
         continue;
      }
      if (use_edisp(srcName)) {
         value += meas_counts_spec[k];
      } else {
         value += true_counts_spec[k];
      }
   }
   return value;
}

void BinnedLikelihood::updateCorrectionFactors(const std::string & srcName, 
					       const SourceMap & sourceMap) const {
   const std::vector<double> & npreds = sourceMap.npreds();
   const Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
   std::vector<double> true_counts_spec;
   const std::vector<std::pair<double,double> >& npred_weights = sourceMap.npred_weights();
   for (size_t k(0); k < energies().size()-1; k++) {
      double my_value(src->pixelCounts(energies().at(k), energies().at(k+1),
				       npreds.at(k)*npred_weights[k].first, npreds.at(k+1)*npred_weights[k].second));
      true_counts_spec.push_back(my_value);
   }
   std::vector<double> meas_counts_spec;
   if (use_edisp(srcName)) {
      std::vector<optimizers::Parameter> pars;
      src->spectrum().getFreeParams(pars);
      const_cast<BinnedLikelihood *>(this)
         ->edisp_correction_factors(srcName, true_counts_spec,
                                    meas_counts_spec);
      m_meas_counts[srcName] = meas_counts_spec;
   } else {
      m_meas_counts[srcName] = true_counts_spec;
   }
   m_true_counts[srcName] = true_counts_spec;
}

// Compute measured count spectrum scaling wrt total true counts
// summed over energy bins.
void BinnedLikelihood::
edisp_correction_factors(const std::string & srcName,
                         const std::vector<double> & true_counts_spec,
                         std::vector<double> & meas_counts_spec) {
   meas_counts_spec.resize(m_energies.size()-1);
   std::vector<double> xi;
   drm().convolve(true_counts_spec, meas_counts_spec);
   // Find the reference pixel in the energy dimension (k) for each
   // sky location (ipix) that contains at least 1 count and for which
   // the true model counts is positive.
   std::map<size_t, size_t> krefs;
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t ipix = m_filledPixels[j] % m_pixels.size();
      size_t k = m_filledPixels[j]/m_pixels.size();
      if (true_counts_spec[k] > 0) {
         krefs[ipix] = k;
      }
   }
   m_krefs[srcName] = krefs;
}

void BinnedLikelihood::getNpreds(const std::string & srcName,
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
BinnedLikelihood::countsSpectrum(const std::string & srcName,
                                 bool use_klims) const {
   size_t kmin(0);
   size_t kmax(m_energies.size() - 1);
   if (use_klims) {
      kmin = m_kmin;
      kmax = m_kmax;
   }
   std::vector<double> counts_spectrum(kmax - kmin, 0);
// Compute ratios of individual source model predictions to that of
// the total model.
   double npred;
   if (!m_modelIsCurrent) {
      npred = computeModelMap();
   }
   std::vector<std::pair<double, double> > modelWts;
   std::pair<double, double> zeros(0, 0);
   modelWts.resize(m_filledPixels.size(), zeros);
   addSourceWts(modelWts, srcName);
   for (size_t j(0); j < m_filledPixels.size(); j++) {
      size_t k(m_filledPixels[j]/m_pixels.size());
      if (k < kmin || k > kmax-1) {
         continue;
      }
      double emin(m_energies[k]);
      double emax(m_energies[k+1]);
      double srcProb = pixelCounts(emin, emax, modelWts[j].first,
                                   modelWts[j].second)/m_model[j];
      size_t indx = m_filledPixels[j];
      counts_spectrum[k - kmin] += srcProb*m_dataMap.data()[indx];
   }
   return counts_spectrum;
}

const std::vector<double> & BinnedLikelihood::fixedModelSpectrum() const {
   return m_fixed_counts_spec;
}

void BinnedLikelihood::computeFixedCountsSpectrum() {
   m_fixed_counts_spec.clear();
   for (size_t k(0); k < m_energies.size()-1; k++) {
      m_fixed_counts_spec.push_back(pixelCounts(m_energies[k], m_energies[k+1],
                                                m_fixedNpreds[k],
                                                m_fixedNpreds[k+1]));
   }
}

Drm & BinnedLikelihood::drm() {
   if (m_drm == 0) {
      m_drm = new Drm(m_dataMap.refDir().ra(), m_dataMap.refDir().dec(), 
                      observation(), m_dataMap.energies());
   }
   return *m_drm;
}

void BinnedLikelihood::set_edisp_flag(bool use_edisp) {
   m_use_edisp = use_edisp;
}

bool BinnedLikelihood::use_edisp(const std::string & srcname) const {
   bool use_src_edisp(m_use_edisp);
   if (srcname != "") {
      /// Determine from Source if it already accounts for convolution
      /// with the energy dispersion, e.g., for Galactic diffuse
      /// background models.
      const std::map<std::string, Source *>::const_iterator 
         & it(m_sources.find(srcname));
      if (it != m_sources.end()) {
         use_src_edisp = it->second->use_edisp();
      }
   }
   return ( (m_use_edisp || ::getenv("USE_BL_EDISP")) && use_src_edisp );
}

void BinnedLikelihood::set_use_single_fixed_map(bool use_sfm) {
   m_use_single_fixed_map = use_sfm;
}

bool BinnedLikelihood::use_single_fixed_map() const {
   return m_use_single_fixed_map;
}

} // namespace Likelihood
