/**
 * @file BinnedLikelihood.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedLikelihood.cxx,v 1.140 2017/10/10 17:24:15 echarles Exp $
 */

#include <cmath>
#include <cstdlib>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "st_stream/StreamFormatter.h"

#include "tip/IFileSvc.h"
#include "tip/Extension.h"

#include "st_facilities/Util.h"
#include "st_facilities/Timer.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/CompositeSource.h"
#include "Likelihood/FitUtils.h"
#include "Likelihood/FileUtils.h"
#include "Likelihood/WeightMap.h"

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
    m_dataCache(dataMap,observation,0,srcMapsFile),
    m_kmin(0), m_kmax(m_dataCache.num_ebins()),
    m_srcMapsFile(srcMapsFile),
    m_config(computePointSources,
	     applyPsfCorrections,
	     performConvolution,
	     resample,
	     resamp_factor,
	     minbinsz,
	     PsfIntegConfig::adaptive,1e-3,1e-6,
	     true,0,true,false,false,false),
    m_drm(new Drm(m_dataCache.countsMap().refDir().ra(), m_dataCache.countsMap().refDir().dec(), 
		  observation, m_dataCache.countsMap().energies())),
    m_srcMapCache(m_dataCache,observation,srcMapsFile,m_config,m_drm),
    m_modelIsCurrent(false),
    m_updateFixedWeights(true){

  std::cerr << "This version of the constructor of BinnedLikelihood is deprecated." << std::endl
	    << "It will be removed in an upcoming release." << std::endl
	    << "Please switch to using the version that takes a BinnedLikeConfigObject." << std::endl;

  m_fixedModelCounts.resize(m_dataCache.nFilled(), 0);
  m_fixed_counts_spec.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_wt.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_edisp.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_edisp_wt.resize(m_dataCache.num_ebins(), 0); 
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
                                   double minbinsz,
				   bool overwriteWeights)
  : LogLike(observation), 
    m_dataCache(dataMap,observation,&weightMap,srcMapsFile, overwriteWeights),
    m_kmin(0),m_kmax(m_dataCache.num_ebins()),
    m_srcMapsFile(srcMapsFile),
    m_config(computePointSources,
	     applyPsfCorrections,
	     performConvolution,
	     resample,
	     resamp_factor,
	     minbinsz,	     
	     PsfIntegConfig::adaptive,1e-3,1e-6,
	     true,0,true,false,false,false),
    m_drm(new Drm(m_dataCache.countsMap().refDir().ra(), m_dataCache.countsMap().refDir().dec(), 
		  observation, m_dataCache.countsMap().energies())),
    m_srcMapCache(m_dataCache,observation,srcMapsFile,m_config,m_drm),
    m_modelIsCurrent(false),
    m_updateFixedWeights(true){    

  std::cerr << "This version of the constructor of BinnedLikelihood is deprecated." << std::endl
	    << "It will be removed in an upcoming release." << std::endl
	    << "Please switch to using the version that takes a BinnedLikeConfigObject." << std::endl;

  m_fixedModelCounts.resize(m_dataCache.nFilled(), 0);
  m_fixed_counts_spec.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_wt.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_edisp.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_edisp_wt.resize(m_dataCache.num_ebins(), 0); 

}

BinnedLikelihood::BinnedLikelihood(CountsMapBase & dataMap,
				   const Observation & observation,	
				   const BinnedLikeConfig& config,
				   const std::string & srcMapsFile,
				   const ProjMap* weightMap,
				   bool overwriteWeights)
  : LogLike(observation), 
    m_dataCache(dataMap,observation,weightMap,srcMapsFile,overwriteWeights),
    m_kmin(0),m_kmax(m_dataCache.num_ebins()),
    m_srcMapsFile(srcMapsFile),
    m_config(config),
    m_drm(new Drm(m_dataCache.countsMap().refDir().ra(), m_dataCache.countsMap().refDir().dec(), 
		  observation, m_dataCache.countsMap().energies(), config.drm_bins())),
    m_srcMapCache(m_dataCache,observation,srcMapsFile,m_config,m_drm),
    m_modelIsCurrent(false),
    m_updateFixedWeights(true){    
  m_fixedModelCounts.resize(m_dataCache.nFilled(), 0);
  m_fixed_counts_spec.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_wt.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_edisp.resize(m_dataCache.num_ebins(), 0); 
  m_fixed_counts_spec_edisp_wt.resize(m_dataCache.num_ebins(), 0); 
}

BinnedLikelihood::~BinnedLikelihood() throw() {
  delete m_drm;
}


void BinnedLikelihood::setCountsMap(const std::vector<float> & counts) {
  m_dataCache.setCountsMap(counts);
  buildFixedModelWts();
}


void BinnedLikelihood::setWeightsMap(const ProjMap* wmap) {
  m_srcMapCache.setWeightsMap(wmap);
  buildFixedModelWts();
}

double BinnedLikelihood::value(const optimizers::Arg & dummy, 
			       bool include_priors) const {
  (void)(dummy);

  // Here we want the weighted verison of the nPred
  double npred = computeModelMap_internal(true);
  
  const std::vector<float> & data = m_dataCache.data( m_dataCache.has_weights() );  

  const std::vector<size_t>& pix_ranges = m_dataCache.firstPixels();
  for (size_t k(m_kmin); k < m_kmax; k++ ) {
    size_t j_start = pix_ranges[k];
    size_t j_stop = pix_ranges[k+1];    
    size_t ipix_base = k * m_dataCache.num_pixels();
    for (size_t j(j_start); j < j_stop; j++) {
      if ( m_model[j] <= 0 ) {
	continue;
      }
      size_t i = ipix_base + m_dataCache.filledPixels()[j];
      double addend = data[i]*std::log(m_model[j]);
      m_accumulator.add(addend);
    }
  }
  m_accumulator.add(-npred);
  
  double my_total(m_accumulator.total());

  if ( include_priors ) {
    /// Add in contribution from priors.
    std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
    for ( ; par != m_parameter.end(); ++par) {
      my_total += par->log_prior_value();
    }
  }
  
  saveBestFit(my_total);
  st_stream::StreamFormatter formatter("BinnedLikelihood", "value", 4);
  formatter.warn() << m_nevals << "  "
		   << my_total << "  "
		   << npred << std::endl;
  m_nevals++;

  return my_total;
}



void BinnedLikelihood::getFreeDerivs(const optimizers::Arg & dummy, 
				     std::vector<double> & derivs, 
				     bool include_priors) const {
  st_stream::StreamFormatter formatter("BinnedLikelihood", "getFreeDerivs",4);
  
  int nparams(getNumFreeParams());
  derivs.resize(nparams, 0);
  double npred(0.);
  if (!m_modelIsCurrent) {
    // here we want the weighted version of npred
    npred = computeModelMap_internal(true);
  }
  const std::vector<float> & data = m_dataCache.data( m_dataCache.has_weights() );  

  /// Update the cached vectors of spectral derivatives inside the
  /// various source maps
  std::vector<Source *> free_srcs;
  std::map<std::string, Source *>::const_iterator src;

  for (src=sources().begin(); src != sources().end(); ++src) {
    if (!std::count(m_fixedSources.begin(), m_fixedSources.end(),
		    src->second->getName())) {
      free_srcs.push_back(src->second);
      std::vector<std::string> parnames;
      src->second->spectrum().getFreeParamNames(parnames);
      SourceMap & srcMap = sourceMap(src->first);
      srcMap.setSpectralValues();
      srcMap.setSpectralDerivs(parnames);
    }
  }

  // /// For debugging/development, keep track of time spent in main loop.
  //st_facilities::Timer timer;
  
  std::vector<Kahan_Accumulator> posDerivs(nparams);
  std::vector<Kahan_Accumulator> negDerivs(nparams);

  long freeIndex(0);

  //timer.start();

  // The data/model is used for each of the deriavtive, so pre-compute it here
  std::vector<double> data_over_model(m_dataCache.nFilled());
  const std::vector<size_t>& pix_ranges = m_dataCache.firstPixels();

  for (size_t k(m_kmin); k < m_kmax; k++ ) {
    size_t j_start = pix_ranges[k];
    size_t j_stop = pix_ranges[k+1];      
    size_t ipix_base = k * m_dataCache.num_pixels();
    for (size_t j(j_start); j < j_stop; j++) {
      size_t ipix = ipix_base + m_dataCache.filledPixels()[j];
      data_over_model[j] = m_model[j] > 0. ? data[ipix] / m_model[j] : 0.;
    }
  }

  // We only need to loop on the free sources
  for (std::vector<Source *>::const_iterator it(free_srcs.begin());
       it != free_srcs.end(); ++it ) {
    

    Source * src(*it);
    SourceMap & srcMap = sourceMap(src->getName());
 
    FitUtils::addFreeDerivs(posDerivs, negDerivs, freeIndex,
			    srcMap, data_over_model, m_dataCache, m_kmin, m_kmax);
    
    // Update index of the next free parameter, based on the number of free parameters of this source
    freeIndex += src->spectrum().getNumFreeParams();

  } // Loop on free sources

  //timer.stop();
  //timer.report("main loop time");
  
   for (size_t i(0); i < derivs.size(); i++) {
    derivs[i] = posDerivs[i].total() + negDerivs[i].total(); 
  }

  /// Derivatives from priors.
  if ( include_priors ) {
    size_t i(0);
    std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
    for ( ; par != m_parameter.end(); ++par) {
      if (par->isFree()) {
	derivs[i] += par->log_prior_deriv();
	i++;
      }
    }
  }
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


  CountsMapBase * BinnedLikelihood::createCountsMap() const {
    CountsMapBase * modelMap = m_dataCache.countsMap().clone();         
    std::vector<float> map;
    computeModelMap(map);
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


  void BinnedLikelihood::addSource(Source * src, bool fromClone, SourceMap* srcMap, bool loadMap) {
    m_bestValueSoFar = -1e38;
    SourceModel::addSource(src, fromClone);
    if ( m_config.use_single_fixed_map() && src->fixedSpectrum()) {
      addFixedSource(src->getName());      
    } else {
      if ( loadMap ) {
        if ( srcMap == 0 ) {
          m_srcMapCache.loadSourceMap(*src,false);
        } else {
          m_srcMapCache.eraseSourceMap(src->getName());
          m_srcMapCache.insertSourceMap(src->getName(),*srcMap);
        }
      }
    }
  }

  void BinnedLikelihood::addSource(Source * src, BinnedLikeConfig* config, bool fromClone) {
    m_bestValueSoFar = -1e38;
    SourceModel::addSource(src, fromClone);
    if ( m_config.use_single_fixed_map() && src->fixedSpectrum()) {
      addFixedSource(src->getName());
    } else {
      m_srcMapCache.loadSourceMap(*src,false,config);
    }
  }


  


  Source * BinnedLikelihood::deleteSource(const std::string & srcName) {
    m_bestValueSoFar = -1e38;
    // Check if this is a fixed source, and if so, remove it from the fixed model.
    std::vector<std::string>::iterator srcIt = 
      std::find(m_fixedSources.begin(), m_fixedSources.end(), srcName);
    bool subtract(true);
    if (source(srcName).fixedSpectrum() && srcIt != m_fixedSources.end()) {
      SourceMap * srcMap(getSourceMap(srcName, false));
      addSourceCounts(m_fixedModelCounts, srcName, srcMap, subtract);
      addFixedNpreds(srcName, srcMap, subtract);
      m_fixedSources.erase(srcIt);
    }
    Source * src(SourceModel::deleteSource(srcName));
    return src;
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
    // call the other version of NpredValue
    return NpredValue(srcName, sourceMap(srcName), weighted);
  }

  // This version forces the recalculation of Npred, whether the source is
  // fixed or not. It is also called from buildFixedModelWts.
  double BinnedLikelihood::NpredValue(const std::string & srcName,
				      SourceMap & sourceMap, 
				      bool weighted) const {
    const Source& src = source(srcName);
    return m_srcMapCache.NpredValue(src,sourceMap,m_kmin,m_kmax,weighted);
  }

  bool BinnedLikelihood::hasSourceMap(const std::string & name) const {
    return m_srcMapCache.hasSourceMap(name);
  }

  SourceMap & BinnedLikelihood::sourceMap(const std::string & name) const {
    const Source& src = source(name);
    SourceMap* srcMap =  m_srcMapCache.getSourceMap(src);
    return *srcMap;
  }

  SourceMap * BinnedLikelihood::getSourceMap(const std::string & srcName, 
					     bool verbose) const {
    const Source& src = source(srcName);
    return m_srcMapCache.getSourceMap(src,verbose);
  }


  SourceMap * BinnedLikelihood::createSourceMap(const std::string & srcName) {
    Source * src = getSource(srcName);
    return m_srcMapCache.createSourceMap(*src);
  }

  SourceMap * BinnedLikelihood::createExternalSourceMap(Source& aSrc) {
    return m_srcMapCache.createSourceMap(aSrc);
  }

  void BinnedLikelihood::eraseSourceMap(const std::string & srcName) {
    m_srcMapCache.eraseSourceMap(srcName);
  }


  void BinnedLikelihood::loadSourceMaps(bool recreate, bool saveMaps) {
    std::vector<std::string> srcNames;
    getSrcNames(srcNames);
    loadSourceMaps(srcNames,recreate,saveMaps);
  }


  void BinnedLikelihood::loadSourceMaps(const  std::vector<std::string>& srcNames,
					bool recreate, bool saveMaps) {  
    std::vector<const Source*> srcs;
    getSources(srcNames,srcs);
    m_srcMapCache.loadSourceMaps(srcs,recreate,saveMaps);
  }

  void BinnedLikelihood::loadSourceMap(const std::string & srcName, bool recreate,
				       bool buildFixedWeights) {

    const Source& src = source(srcName);
    m_srcMapCache.loadSourceMap(src,recreate);
 
    std::vector<std::string>::iterator srcIt = 
      std::find(m_fixedSources.begin(), m_fixedSources.end(), srcName);
    if (srcIt != m_fixedSources.end() && buildFixedWeights) {
      buildFixedModelWts();
    }
  }
  

  void BinnedLikelihood::setSourceMapImage(const std::string & name,
					   const std::vector<float>& image) {
    const Source& src = source(name);
    m_srcMapCache.setSourceMapImage(src,image);
  }

  void BinnedLikelihood::saveSourceMaps(const std::string & filename,
					bool replace,
					bool includecountsmap) {
    if (filename != "") {
      m_srcMapsFile = filename;
    }

    if (!st_facilities::Util::fileExists(filename)) {
      if ( includecountsmap ) {
	m_dataCache.countsMap().writeOutput("BinnedLikelihood", m_srcMapsFile);
      } else { 
	m_dataCache.countsMap().writeEmptyOutput("BinnedLikelihood", m_srcMapsFile);	
      }
    }
    
    std::vector<std::string> srcNames;
    getSrcNames(srcNames);
    
    std::vector<const Source*> srcs;
    getSources(srcNames,srcs);

    m_srcMapCache.saveSourceMaps(m_srcMapsFile,srcs,replace);

    tip::Extension* whdu = saveWeightsMap(replace);
    if ( whdu != 0 ) {
      delete whdu;
    }

  }
  

  void BinnedLikelihood::saveSourceMap_partial(const std::string & filename,
						const Source& source,
						int kmin, int kmax,
						bool replace) {
    if (!st_facilities::Util::fileExists(filename)) {
      m_dataCache.countsMap().writeEnergies("BinnedLikelihood", filename, kmin, kmax);
    }

    m_srcMapCache.saveSourceMap_partial(filename,source,kmin,kmax,replace);
  }

  

  tip::Extension* BinnedLikelihood::saveTotalFixedSourceMap(bool replace) {
    // FIXME

    std::vector<float> srcMap(source_map_size(),0.);
    std::vector<std::string> fixedSrcNames;
    std::vector<const Source*> fixedSrcs;
    for ( std::map<std::string, Source*>::const_iterator itr = m_sources.begin();
	  itr != m_sources.end(); itr++ ) {
      if (itr->second->fixedSpectrum()) {
	fixedSrcNames.push_back(itr->first);
	fixedSrcs.push_back(itr->second);
      }
    }
    fillSummedSourceMap(fixedSrcNames,srcMap);
    bool has_fixed = FileUtils::fileHasExtension(m_srcMapsFile, "__FIXED__");    
    tip::Extension* ext(0);
    if ( has_fixed ) {
      if ( replace ) {
	ext = FileUtils::replace_image_from_float_vector(m_srcMapsFile,"__FIXED__",
							 m_dataCache.countsMap(),srcMap,true);
      } else {
	// just leave it be;
	;
      }
    } else {
      ext = FileUtils::append_image_from_float_vector(m_srcMapsFile,"__FIXED__",
						      m_dataCache.countsMap(),srcMap,true);
    }
    
    tip::Extension* pext = FileUtils::write_model_parameters_to_table(m_srcMapsFile,"__FIXED_PAR__",
								      fixedSrcs);
    delete pext;
    return ext;
  }


  double BinnedLikelihood::computeModelMap_internal(bool weighted) const {
    using namespace std;
    std::vector<double> modelCounts;
    modelCounts.resize(m_dataCache.nFilled(), 0.);

    if (fixedModelUpdated() && m_updateFixedWeights) {
      const_cast<BinnedLikelihood *>(this)->buildFixedModelWts();
    }

    // We don't apply the weights to the model we will be filling
    // But we do apply the energy dispersion
    // This was already handled in buildFixedModelWts
    std::copy(m_fixedModelCounts.begin(), m_fixedModelCounts.end(), modelCounts.begin());

    std::vector<std::string> srcNames;
    getSrcNames(srcNames);
    
    // Get the version of the fixed model counts to use
    // Always pick one of the two that uses energy dispersion
    const std::vector<double>& fixedModelCounts = weighted ? m_fixed_counts_spec_edisp_wt : m_fixed_counts_spec_edisp;

    double npred(0.);
    double npred_check(0.);
    for ( size_t kx(m_kmin); kx < m_kmax; kx++ ) {
      npred += fixedModelCounts[kx];
    }

    for (size_t i(0); i < srcNames.size(); i++) {     
      // EAC FIXME, in principle we should only need to call NpredValue
      // (which updates and caches the computation for the Npred) for the free sources
      // but because it is possible to free and modify the source parameters directly, 
      // that doesn't always work. 
      // So we have to call NpredValue on all the sources 
      double npred_src = NpredValue(srcNames[i], weighted);
      npred_check += npred_src;
      if (std::count(m_fixedSources.begin(), m_fixedSources.end(),
		     srcNames[i]) == 0) {
	addSourceCounts(modelCounts, srcNames[i]);
	npred += npred_src;
      } 
    }

    m_model.clear();
    m_model.resize(m_dataCache.nFilled(), 0);

    const std::vector<size_t>& pix_ranges = m_dataCache.firstPixels();
 
    for (size_t k(m_kmin); k < m_kmax; k++ ) {
      size_t j_start = pix_ranges[k];
      size_t j_stop = pix_ranges[k+1];      
      for (size_t j(j_start); j < j_stop; j++) {
	m_model[j] = modelCounts[j];
      }
    }
    m_modelIsCurrent = true;

    return npred;
  }

  void BinnedLikelihood::computeModelMap(std::vector<float> & modelMap,
					 bool use_mask) const {
    std::vector<std::string> srcNames;
    getSrcNames(srcNames);
    computeModelMap(srcNames,modelMap,use_mask);
  }
  
  void BinnedLikelihood::computeModelMap(const std::string & srcName, 
					 std::vector<float> & modelMap,
					 bool use_mask) const {
    const Source& src = source(srcName);
    m_srcMapCache.computeModelMap(src,modelMap,use_mask);			  
  }

  void BinnedLikelihood::computeModelMap(const std::vector<std::string>& srcNames, 
					 std::vector<float> & modelMap,
					 bool use_mask) const {
    std::vector<const Source*> srcs;
    getSources(srcNames,srcs);
    m_srcMapCache.computeModelMap(srcs,modelMap,use_mask);
  }


  void BinnedLikelihood::updateModelMap(std::vector<float> & modelMap,
					SourceMap * srcMap,
					bool use_mask) const {
    const Source& src = source(srcMap->name());
    m_srcMapCache.updateModelMap(modelMap,src,srcMap,use_mask);
  }

  void BinnedLikelihood::updateModelMap_fromSrcMap(std::vector<float> & modelMap,
						   const Source& src,
						   SourceMap * srcMap,
						   bool use_mask) const {
    m_srcMapCache.updateModelMap_fromSrcMap(modelMap,src,srcMap,use_mask);
  }


  void BinnedLikelihood::getNpreds(const std::string & srcName,
				   std::vector<double> & npreds) const {
    const Source& src = source(srcName);
    m_srcMapCache.getNpreds(src,npreds);
  }


  const std::vector<double>& 
  BinnedLikelihood::modelCountsSpectrum(const std::string & srcName,
					bool weighted) const {
    const Source& src = source(srcName);
    return m_srcMapCache.modelCountsSpectrum(src, weighted);
  }


  bool BinnedLikelihood::fixedModelUpdated() const {
    // Check if the fixed list has changed.

    // Check to see if we have build the m_fixed_counts_spec yet
    if ( m_fixed_counts_spec.size() == 0 && m_fixedSources.size() != 0 ) {
      return true;
    }

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
      const Source* src = srcIt->second;
      if (src->fixedSpectrum()) {
	const SourceMap* srcMap = m_srcMapCache.getSourceMap(*src);
	if ( srcMap == 0 ) {
	  throw std::runtime_error("BinnedLikelihood::fixedModelUpdated: "
				   "inconsistent m_srcMaps.");
	}
	if ( srcMap->spectrum_changed() ) {	  
	  return true;
	}	
      }
    }
    return false;
  }




  void BinnedLikelihood::buildFixedModelWts(bool process_all) {
//    std::cout << "Entering BinnedLikelihood::buildFixedModelWts(" << std::boolalpha << process_all << ") - m_fixedModelBuilt = " 
//          << m_fixedModelBuilt << std::noboolalpha << std::endl;
    //if (m_fixedModelBuilt) return;
    m_fixedSources.clear();

    m_fixedModelCounts.clear();
    m_fixedModelCounts.resize(m_dataCache.nFilled(), 0.);

    m_fixed_counts_spec.clear();
    m_fixed_counts_spec_wt.clear();
    m_fixed_counts_spec_edisp.clear();
    m_fixed_counts_spec_edisp_wt.clear();

    m_fixed_counts_spec.resize(m_dataCache.num_ebins(), 0); 
    m_fixed_counts_spec_wt.resize(m_dataCache.num_ebins(), 0); 
    m_fixed_counts_spec_edisp.resize(m_dataCache.num_ebins(), 0); 
    m_fixed_counts_spec_edisp_wt.resize(m_dataCache.num_ebins(), 0); 

    std::map<std::string, Source *>::const_iterator srcIt(m_sources.begin());
    for ( ; srcIt != m_sources.end(); ++srcIt) {
//      std::cout << "BinnedLikelihood::buildFixedModelWts(): processing source " << srcIt->first << std::endl;
      const std::string & srcName(srcIt->first);
      const Source* src = srcIt->second;      
      if (src->fixedSpectrum() ) {
          //if (m_fixedModelBuilt) continue;
          if ( ! process_all ) {
            addFixedSource(srcName);
          } else {
            m_srcMapCache.getSourceMap(*src, false);
          }
      } else { 
          // Process non-fixed sources.
          //
          // Ensure model map is available.
          SourceMap * srcMap = m_srcMapCache.getSourceMap(*src, false);
      }
    }

    computeFixedCountsSpectrum();
//    std::cout << "Leaving BinnedLikelihood::buildFixedModelWts()" << std::endl;  
  }


  void BinnedLikelihood::fillSummedSourceMap(const std::vector<std::string>& sources, std::vector<float>& model) {
    std::vector<const Source*> srcs;
    getSources(sources,srcs);
    m_srcMapCache.fillSummedSourceMap(srcs,model);
  }

   void BinnedLikelihood::fillSingleSourceMap(const std::string& sourceName, 
					      std::vector<float>& model,
					      FileUtils::SrcMapType& mapType,
					      int kmin, int kmax) {
     const Source& src = source(sourceName);
     m_srcMapCache.fillSingleSourceMap(src, model, mapType, kmin, kmax);
  }


  void BinnedLikelihood::addFixedSource(const std::string & srcName) {
    //std::cout << "BinnedLikelihood::addFixedSource() - processing source " << srcName << std::endl;
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
    const Source& src = *(srcIt->second);

    SourceMap * srcMap = m_srcMapCache.getSourceMap(src);
    bool has_wts = srcMap->weights() != 0;
    if ( m_config.save_all_srcmaps() ) {
      srcMap->setSaveModel(true);
    }
    addSourceCounts(m_fixedModelCounts, srcName, srcMap, false, true);
    addFixedNpreds(srcName, srcMap, false);

    // Remove this source from the stored source maps to save memory
    if ( !srcMap->save_model() ) {
      if(srcMap->clear_model( m_config.delete_local_fixed() )){
        m_fixedModelBuilt = true;
      }
      
    }
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
    SourceMap * srcMap = getSourceMap(srcName, false);
    srcMap->reloadIfCleared();
    //bool has_wts = srcMap->weights() != 0;
    if (srcMap == 0) {
      throw std::runtime_error("SourceMap cannot be created for " + srcName);
    }

    // Subtract the source weights from the summed fixed model weights.
    bool subtract;
    addFixedNpreds(srcName, srcMap, subtract=true);
    addSourceCounts(m_fixedModelCounts, srcName, srcMap, subtract=true);
    // Remove from the list of fixed sources
    std::vector<std::string>::iterator it 
      = std::find(m_fixedSources.begin(), m_fixedSources.end(), srcName);
    m_fixedSources.erase(it);
  }



  std::vector<double> 
  BinnedLikelihood::countsSpectrum(const std::string & srcName,
				   bool use_klims) const {

    size_t kmin(0);
    size_t kmax(m_dataCache.num_ebins());
    if (use_klims) {
      kmin = m_kmin;
      kmax = m_kmax;
    }
    std::vector<double> counts_spectrum(kmax - kmin, 0);
    // Compute ratios of individual source model predictions to that of
    // the total model.
    double npred;
    if (!m_modelIsCurrent) {
      npred = computeModelMap_internal();
    }
    std::vector<double> modelCounts(m_dataCache.nFilled(), 0.);
    addSourceCounts(modelCounts, srcName);

    const std::vector<size_t>& pix_ranges = m_dataCache.firstPixels(); 
    for (size_t k(kmin); k < kmax; k++ ) {
      size_t j_start = pix_ranges[k];
      size_t j_stop = pix_ranges[k+1];  
      size_t ipix_base = k * m_dataCache.num_pixels();
      for (size_t j(j_start); j < j_stop; j++) {
	double srcProb = modelCounts[j] / m_model[j];
	size_t indx = ipix_base + m_dataCache.filledPixels()[j];
	counts_spectrum[k - kmin] += srcProb*m_dataCache.data()[indx];
      }
    }
    return counts_spectrum;
  }


  int BinnedLikelihood::edisp_val(const std::string & srcname) const {
    int src_edisp_val = m_config.edisp_val();
    if ( src_edisp_val == 0 ) return src_edisp_val;    
    if (srcname != "") {
      /// Determine from Source if it already accounts for convolution
      /// with the energy dispersion, e.g., for Galactic diffuse
      /// background models.
      const std::map<std::string, Source *>::const_iterator 
	& it(m_sources.find(srcname));
      if (it != m_sources.end()) {
	if ( ! it->second->use_edisp() ) {
	  src_edisp_val = 0;
	}
      }
    }
    return src_edisp_val;
  }


  Drm & BinnedLikelihood::drm() {
    if (m_drm == 0) {
      m_drm = new Drm(m_dataCache.countsMap().refDir().ra(), m_dataCache.countsMap().refDir().dec(), 
		      observation(), m_dataCache.countsMap().energies());
    }
    return *m_drm;
  }


  void BinnedLikelihood::initialize_composite(CompositeSource& comp) const {
    comp.buildSourceMapCache(m_dataCache,m_srcMapsFile,m_drm);
  }


  tip::Extension* BinnedLikelihood::saveWeightsMap(bool replace) const {
    return m_dataCache.saveWeightsMap(m_srcMapsFile,replace);
  }


  double BinnedLikelihood::spectrum(const Source * src, double energy)  {
    const optimizers::Function & spectrum(src->spectrum());
    optimizers::dArg eArg(energy);
    return spectrum(eArg);
  }


  void BinnedLikelihood::setImageDimensions(tip::Image * image, 
					    long * dimensions) {
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

 
  void BinnedLikelihood::addSourceCounts(std::vector<double> & modelCounts,
					 const std::string & srcName,
					 SourceMap * srcMap,
					 bool subtract,
					 bool latchParams) const {

    
    const Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
    if (src == 0) {
      return;
    }
    if (modelCounts.size() != m_dataCache.nFilled()) {
      throw std::runtime_error("BinnedLikelihood::addSourceCounts: "
			       "modelCounts size does not match "
			       "number of filled pixels.");
    }
    SourceMap * sourceMap = srcMap;
    if ( sourceMap == 0 ) {
      sourceMap = getSourceMap(srcName);
    } 
    
    sourceMap->reloadIfCleared();
    sourceMap->setSpectralValues(latchParams);    

    FitUtils::addSourceCounts(modelCounts,*sourceMap,
			      m_dataCache, subtract);

    if ( !sourceMap->save_model() ) {
      if (std::count(m_fixedSources.begin(), m_fixedSources.end(), srcName) != 0) {
        sourceMap->clear_model( m_config.delete_local_fixed() );
      }  
    }

  }

  void BinnedLikelihood::addFixedNpreds(const std::string & srcName,
					SourceMap * srcMap, 
					bool subtract) {

    FitUtils::addFixedNpreds(m_fixed_counts_spec, m_fixed_counts_spec_wt, 
			     m_fixed_counts_spec_edisp, m_fixed_counts_spec_edisp_wt,
			     *srcMap, m_dataCache, subtract);
    
  }

  float BinnedLikelihood::npred_explicit(const BinnedLikelihood& like,
					 bool weighted) {
    std::vector<float> model(like.dataCache().data_map_size());
    like.computeModelMap(model);
    float retVal(0.);
    if ( weighted ) {
      const WeightMap* wts = like.weightMap();
      if ( wts == 0 ) {
	throw std::runtime_error("BinnedLikelihood::npred_explicit called with weighted option, but no weights are present");
	return 0;
      }
      const std::vector<float>& wt_vals = wts->model();
      FitUtils::innerProduct(model.begin(), model.end(), wt_vals.begin(), wt_vals.end(), retVal);
      /*
      size_t npix = like.dataCache().num_pixels();
      size_t ne = like.dataCache().num_ebins();
      size_t idx_start(0);
      for ( size_t ie(0); ie < ne; ie++, idx_start+=npix) {
	float val_k(0.);
	float val_k_nowts(0.);
	FitUtils::innerProduct(model.begin()+idx_start, model.begin()+idx_start+npix, 
			       wt_vals.begin()+idx_start, wt_vals.begin()+idx_start+npix, val_k);
	FitUtils::sumVector(model.begin()+idx_start, model.begin()+idx_start+npix, val_k_nowts);
	std::cout << "npred_explicit:" << ' ' << ie << ' ' << val_k << ' ' << val_k_nowts << std::endl;
      }
      */
    } else {
      FitUtils::sumVector(model.begin(), model.end(), retVal);
    }
    return retVal;
  }
     
  float BinnedLikelihood::npred_explicit_src(const BinnedLikelihood& like, 
					     const std::string& srcName, 
					     bool weighted){
    std::vector<float> model(like.dataCache().data_map_size());
    like.computeModelMap(srcName, model);
    float retVal(0.);
    if ( weighted ) {
      const WeightMap* wts = like.weightMap();
      if ( wts == 0 ) {
	throw std::runtime_error("BinnedLikelihood::npred_explicit called with weighted option, but no weights are present");
	return 0;
      }
      const std::vector<float>& wt_vals = wts->model();
      FitUtils::innerProduct(model.begin(), model.end(), wt_vals.begin(), wt_vals.end(), retVal);
    } else {
      FitUtils::sumVector(model.begin(), model.end(), retVal);
    }
    return retVal;
  }
  
  double BinnedLikelihood::nll_explict(const BinnedLikelihood& like, 
				       bool weighted){
    std::vector<float> model(like.dataCache().data_map_size());
    like.computeModelMap(model);
    const std::vector<float>& data = like.countsMap().data();
    double retVal(0.);
    if ( weighted ) {
      const WeightMap* wts = like.weightMap();
      if ( wts == 0 ) {
	throw std::runtime_error("BinnedLikelihood::npred_explicit called with weighted option, but no weights are present");
	return 0;
      }
      const std::vector<float>& wt_vals = wts->model();
      retVal = FitUtils::logLikePoisson(data.begin(), data.end(), model.begin(), model.end(), 
					wt_vals.begin(), wt_vals.end());
    } else {
      retVal = FitUtils::logLikePoisson(data.begin(), data.end(), model.begin(), model.end());
    }
    return retVal;
  }

  void BinnedLikelihood::spectrum_explict(const BinnedLikelihood& like, 
					  const std::string& srcName, 
					  std::vector<double>& spec,
					  bool weighted) {

    std::vector<float> model(like.dataCache().data_map_size());
    const WeightMap* wts(0);
    if ( weighted ) {
      wts = like.weightMap();
      if ( wts == 0 ) {
	throw std::runtime_error("BinnedLikelihood::spectrum_explicit called with weighted option, but no weights are present");
	return;
      }
    }
    like.computeModelMap(srcName, model);
    size_t npix = like.dataCache().num_pixels();
    size_t ne = like.dataCache().num_ebins();
    size_t i_start(0);
    size_t i_end(npix);
    spec.clear();
    spec.resize(ne, 0);
    for ( size_t ie(0); ie < ne; ie++ ) {
      float sum(0.);
      for ( size_t i(i_start); i < i_end; i++ ) { 
	if ( wts == 0 ) {
	  sum += model[i];
	} else { 
	  sum += model[i] * wts->model()[i];
	}
      }
      i_start += npix;
      i_end += npix; 
      spec[ie] = sum;
    }   
  }

  void BinnedLikelihood::check_npreds(const BinnedLikelihood& like, 
				      bool weighted) {
    std::vector<std::string> srcNames;
    like.getSrcNames(srcNames);
    for ( std::vector<std::string>::const_iterator itr = srcNames.begin(); itr != srcNames.end();
	  itr++ ) {
      double npred_func = like.NpredValue(*itr, weighted);
      float npred_explicit = BinnedLikelihood::npred_explicit_src(like, * itr, weighted);
      double npred_check = npred_func - npred_explicit;
      if ( std::fabs(npred_check) > 0.01 ) {
	std::cout << "BinnedLikelihood::check_npreds " << *itr << ' ' 
		  << npred_func << ' ' << npred_explicit << std::endl;
	std::vector<double> spec_explicit;
	int edisp_src = like.edisp_val(*itr);
	SourceMap& srcMap = like.sourceMap(*itr);
	const std::vector<double>& spec_func = srcMap.counts_spectra(edisp_src, weighted);
	BinnedLikelihood::spectrum_explict(like, *itr, spec_explicit, weighted);
	const std::vector<std::vector<std::pair<double,double> > >& w_npreds = srcMap.weighted_npreds();
	for ( size_t ie(0); ie < spec_func.size(); ie++ ) {
	  std::cout << ie << ' ' << spec_func[ie] << ' ' << spec_explicit[ie];
	  const std::vector<std::pair<double,double> >& ww_npreds = w_npreds[ie];
	  for ( size_t iw(0); iw < ww_npreds.size(); iw++ ) {
	    std::cout << ' ' << ww_npreds[iw].first << ' ' << ww_npreds[iw].second;
	  }
	  std::cout << std::endl;
	}
      }
    }
  }


} // namespace Likelihood
