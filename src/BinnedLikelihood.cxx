/**
 * @file BinnedLikelihood.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedLikelihood.cxx,v 1.139 2017/10/06 01:38:12 echarles Exp $
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
	     true,false,true,false,false,false),
    m_drm(new Drm(m_dataCache.countsMap().refDir().ra(), m_dataCache.countsMap().refDir().dec(), 
		  observation, m_dataCache.countsMap().energies())),
    m_srcMapCache(m_dataCache,observation,srcMapsFile,m_config,m_drm),
    m_modelIsCurrent(false),
    m_updateFixedWeights(true){
  m_fixedModelWts.resize(m_dataCache.nFilled(), std::make_pair(0, 0));
  m_fixedNpreds.resize(m_dataCache.num_energies(), 0); 
  m_fixedNpred_xis.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
  m_fixedNpred_wts.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
  m_fixedNpred_xiwts.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
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
	     true,false,true,false),
    m_drm(new Drm(m_dataCache.countsMap().refDir().ra(), m_dataCache.countsMap().refDir().dec(), 
		  observation, m_dataCache.countsMap().energies())),
    m_srcMapCache(m_dataCache,observation,srcMapsFile,m_config,m_drm),
    m_modelIsCurrent(false),
    m_updateFixedWeights(true){    
  m_fixedModelWts.resize(m_dataCache.nFilled(), std::make_pair(0, 0));
  m_fixedNpreds.resize(m_dataCache.num_energies(), 0);
  m_fixedNpred_xis.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
  m_fixedNpred_wts.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
  m_fixedNpred_xiwts.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
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
		  observation, m_dataCache.countsMap().energies())),
    m_srcMapCache(m_dataCache,observation,srcMapsFile,m_config,m_drm),
    m_modelIsCurrent(false),
    m_updateFixedWeights(true){    
  m_fixedModelWts.resize(m_dataCache.nFilled(), std::make_pair(0, 0));
  m_fixedNpreds.resize(m_dataCache.num_energies(), 0);
  m_fixedNpred_xis.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
  m_fixedNpred_wts.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
  m_fixedNpred_xiwts.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
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

double BinnedLikelihood::value(optimizers::Arg & dummy) const {
  (void)(dummy);

  // Here we want the weighted verison of the nPred
  double npred = computeModelMap_internal(true);
  
  const std::vector<float> & data = m_dataCache.data( m_dataCache.has_weights() );
  
  for (size_t i(0); i < m_dataCache.nFilled(); i++) {
    if (m_model.at(i) > 0) {
      size_t j(m_dataCache.filledPixels()[i]);
      size_t k(j/num_pixels());
      if (k >= m_kmin && k <= m_kmax) {
	double addend = data.at(j)*std::log(m_model[i]);
	m_accumulator.add(addend);
      }
    }
  }
  m_accumulator.add(-npred);
  
  double my_total(m_accumulator.total());



  /// Add in contribution from priors.
  std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
  for ( ; par != m_parameter.end(); ++par) {
    my_total += par->log_prior_value();
  }
  
  saveBestFit(my_total);
  st_stream::StreamFormatter formatter("BinnedLikelihood", "value", 4);
  formatter.warn() << m_nevals << "  "
		   << my_total << "  "
		   << npred << std::endl;
  m_nevals++;

  return my_total;
}



void BinnedLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
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
      srcMap.setSpectralValues(m_dataCache.energies());
      srcMap.setSpectralDerivs(m_dataCache.energies(),parnames);
    }
  }

  // /// For debugging/development, keep track of time spent in main loop.
  //st_facilities::Timer timer;
  
  std::vector<Kahan_Accumulator> posDerivs(nparams);
  std::vector<Kahan_Accumulator> negDerivs(nparams);

  long freeIndex(0);

  //timer.start();

  // We only need to loop on the free sources
  for (std::vector<Source *>::const_iterator it(free_srcs.begin());
       it != free_srcs.end(); ++it ) {

    Source * src(*it);
    bool use_edisp_val = m_srcMapCache.use_edisp(src);

    std::string srcName = src->getName();
    SourceMap & srcMap = sourceMap(srcName);
    bool has_wts = srcMap.weights() != 0;
    const std::vector< std::vector<double> > & specDerivs = srcMap.cached_specDerivs();

    // We need this stuff for the second term...
    const Drm_Cache* drm_cache = srcMap.drm_cache();
    const std::vector<double> & npreds = srcMap.npreds();
    const std::vector<std::pair<double,double> > & npred_weights =  srcMap.npred_weights();

    // First term, derivate of n_obs log n_model = ( n_obs / n_model ) * ( d model / d param ) 
    // loop over all the filled pixels
    for (size_t j(0); j < m_dataCache.nFilled(); j++) {
      size_t jmin(m_dataCache.filledPixels()[j]);
      size_t k(jmin/num_pixels());
      size_t ipix(jmin % num_pixels());
      if (k < m_kmin || k > m_kmax-1) {
	continue;
      }
      double emin(m_dataCache.energies()[k]);
      double emax(m_dataCache.energies()[k+1]);
      double log_e_ratio(m_dataCache.log_energy_ratios()[k]);
 
      // EAC, deal with energy dispersion 
      // FIXME, do we really have to deal with kref here?
      int kref(-1);	
      double xi(1.);
      if (use_edisp_val) {
	xi = drm_cache->get_correction(k,kref,has_wts);
	if ( kref < 0 ) {
	  // Still have true counts in this energy bin, just use kk as kref
	  kref = k;
	}
      } else {
	kref = k;
      }

      size_t jref(kref*num_pixels() + ipix);
      size_t jmax(jref + num_pixels());

      if (m_model.at(j) > 0) {
	long iparam(freeIndex);
     	for (size_t i(0); i < specDerivs.size(); i++, iparam++) {	  
	  double my_deriv = xi* pixelCounts(emin,emax,
					    srcMap[jref]*specDerivs[i][k],
					    srcMap[jmax]*specDerivs[i][k+1],
					    log_e_ratio);
	  double addend = (data.at(jmin)/m_model.at(j))*my_deriv;
	  if (addend > 0) {
	    posDerivs[iparam].add(addend);
	  } else {
	    negDerivs[iparam].add(addend);
	  }
	} 	
      }
    }
    
    // Second term, the derivatives of the nPreds. 
    // Loop over the energy layers
    for (size_t kk(0); kk < m_dataCache.num_ebins(); kk++) {
      
      if (kk < m_kmin || kk > m_kmax-1) {
	continue;
      }
      double emin(m_dataCache.energies()[kk]);
      double emax(m_dataCache.energies()[kk+1]);
      double log_e_ratio(m_dataCache.log_energy_ratios()[kk]);

      // EAC, deal with energy dispersion
      int kref(-1);	
      double xi(1.);
      if (use_edisp_val) {
	xi = drm_cache->get_correction(kk,kref,has_wts);
	if ( kref < 0 ) {
	  // Still have true counts in this energy bin, just use kk as kref
	  kref = kk;
	}
      } else {
	kref = kk;
      }

      long iparam2(freeIndex);
      for (size_t i(0); i < specDerivs.size(); i++, iparam2++) {

	double addend = xi*pixelCounts(emin,emax,
				       npreds[kref]*npred_weights[kref].first*specDerivs[i][kk],
				       npreds[kref+1]*npred_weights[kref].second*specDerivs[i][kk+1],
				       log_e_ratio);
	if (-addend > 0) {
	  posDerivs[iparam2].add(-addend);
	} else {
	  negDerivs[iparam2].add(-addend);
	}
      }
    } // Loop on energy bins

    // Update index of the next free parameter, based on the number of free parameters of this source
    freeIndex += src->spectrum().getNumFreeParams();

  } // Loop on free sources

  //timer.stop();
  //timer.report("main loop time");
  
 
  for (size_t i(0); i < derivs.size(); i++) {
    derivs[i] = posDerivs[i].total() + negDerivs[i].total(); 
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
    std::vector<float> map;
    computeModelMap(map);
    CountsMapBase * modelMap = m_dataCache.countsMap().clone();         
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


  void BinnedLikelihood::addSource(Source * src, bool fromClone, SourceMap* srcMap) {
    m_bestValueSoFar = -1e38;
    SourceModel::addSource(src, fromClone);
    if ( m_config.use_single_fixed_map() && src->fixedSpectrum()) {
      addFixedSource(src->getName());
    } else {
      if ( srcMap == 0 ) {
	m_srcMapCache.loadSourceMap(*src,false);
      } else {
	m_srcMapCache.insertSourceMap(src->getName(),*srcMap);
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
      addSourceWts(m_fixedModelWts, srcName, srcMap, subtract);
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

    std::vector<std::pair<double, double> > modelWts;
    modelWts.resize(m_dataCache.nFilled());
  
    if (fixedModelUpdated() && m_updateFixedWeights) {
      const_cast<BinnedLikelihood *>(this)->buildFixedModelWts();
    }

    for (size_t j(0); j < m_fixedModelWts.size(); j++) {
      modelWts.at(j).first = m_fixedModelWts.at(j).first;
      modelWts.at(j).second = m_fixedModelWts.at(j).second;
    }	   

    std::vector<std::string> srcNames;
    getSrcNames(srcNames);
    
    double npred(0.);

    for (size_t i(0); i < srcNames.size(); i++) {      
      npred += NpredValue(srcNames[i],weighted);
      if (std::count(m_fixedSources.begin(), m_fixedSources.end(),
		     srcNames.at(i)) == 0) {
	addSourceWts(modelWts, srcNames[i]);
      }
    }

    m_model.clear();
    m_model.resize(m_dataCache.nFilled(), 0);
    for (size_t j(0); j < m_dataCache.nFilled(); j++) {
      size_t k(m_dataCache.filledPixels()[j]/num_pixels());
      if (k < m_kmin || k > m_kmax-1) {
	continue;
      }
      double emin(m_dataCache.energies()[k]);
      double emax(m_dataCache.energies()[k+1]);      
      m_model.at(j) = pixelCounts(emin, emax, modelWts.at(j).first,
				  modelWts.at(j).second, m_dataCache.log_energy_ratios()[k]);
    }
  
    m_modelIsCurrent = true;

    return npred;
  }

  void BinnedLikelihood::computeModelMap(std::vector<float> & modelMap) const {
    std::vector<std::string> srcNames;
    getSrcNames(srcNames);
    computeModelMap(srcNames,modelMap);
  }
  
  void BinnedLikelihood::computeModelMap(const std::string & srcName, 
					 std::vector<float> & modelMap) const {
    const Source& src = source(srcName);
    m_srcMapCache.computeModelMap(src,modelMap);			  
  }

  void BinnedLikelihood::computeModelMap(const std::vector<std::string>& srcNames, 
					 std::vector<float> & modelMap) const {
    std::vector<const Source*> srcs;
    getSources(srcNames,srcs);
    m_srcMapCache.computeModelMap(srcs,modelMap);
  }


  void BinnedLikelihood::updateModelMap(std::vector<float> & modelMap,
					SourceMap * srcMap) const {
    const Source& src = source(srcMap->name());
    m_srcMapCache.updateModelMap(modelMap,src,srcMap);
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
    m_fixedSources.clear();
    m_fixedModelWts.clear();
    m_fixedModelWts.resize(m_dataCache.nFilled(), std::make_pair(0, 0));
    m_fixedNpreds.clear();
    m_fixedNpreds.resize(m_dataCache.energies().size(), 0);
    m_fixedNpred_xis.clear();
    m_fixedNpred_xis.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
    m_fixedNpred_wts.clear();
    m_fixedNpred_wts.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));
    m_fixedNpred_xiwts.clear();
    m_fixedNpred_xiwts.resize(m_dataCache.num_ebins(), std::make_pair(0, 0));

    std::map<std::string, Source *>::const_iterator srcIt(m_sources.begin());
    for ( ; srcIt != m_sources.end(); ++srcIt) {
      const std::string & srcName(srcIt->first);
      const Source* src = srcIt->second;      
      if (src->fixedSpectrum() && !process_all) {
	addFixedSource(srcName);
      } else { 
	// Process non-fixed sources.
	//
	// Ensure model map is available.
	SourceMap * srcMap = m_srcMapCache.getSourceMap(*src, false);
      }
    }
  
    computeFixedCountsSpectrum();
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
    
    addSourceWts(m_fixedModelWts, srcName, srcMap, false, true);

    bool use_edisp_val = use_edisp(srcName);
    const Drm_Cache* drm_cache = srcMap->drm_cache();
    const std::vector<double>& spec = srcMap->specVals();
    const std::vector<double>& npred_vals = srcMap->npreds();
    const std::vector<std::pair<double,double> >& npred_weights = srcMap->npred_weights();

    double xi(1);
    int kref(-1);
    size_t ne(m_dataCache.num_ebins());

    for (size_t k(0); k < ne; k++) {
      double npred_0 = npred_vals[k] * spec[k];
      double npred_1 = npred_vals[k+1] * spec[k+1];
      m_fixedNpreds[k] += npred_0;
      m_fixedNpred_wts[k].first += (npred_0 * npred_weights[k].first);
      m_fixedNpred_wts[k].second += (npred_1 * npred_weights[k].second);
      
      if (use_edisp_val) {
	xi = drm_cache->get_correction(k,kref,has_wts);
	if ( kref < 0 ) { 
	  // Still have true counts in this energy bin, just use k as kref
	  kref = k;
	} 
      } else {
	kref = k;
      }
      if ( kref != 0 ) {
	npred_0 = npred_vals[kref] * spec[kref];
	npred_1 = npred_vals[kref+1] * spec[kref+1];
      }
      m_fixedNpred_xis[k].first += (npred_0*xi);
      m_fixedNpred_xis[k].second += (npred_1*xi);
      m_fixedNpred_xiwts[k].first += (npred_0*xi*npred_weights[kref].first);
      m_fixedNpred_xiwts[k].second += (npred_1*xi*npred_weights[kref].first);
    }

    // Add the last bin edge
    m_fixedNpreds[ne] += (npred_vals[ne] * spec[ne]);

    // Remove this source from the stored source maps to save memory
    if ( !srcMap->save_model() ) {
      srcMap->clear_model();
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
    bool has_wts = srcMap->weights() != 0;
    if (srcMap == 0) {
      throw std::runtime_error("SourceMap cannot be created for " + srcName);
    }

    bool use_edisp_val = use_edisp(srcName);
    const Drm_Cache* drm_cache = srcMap->drm_cache();
    const std::vector<double>& spec = srcMap->specVals();
    const std::vector<double>& npred_vals = srcMap->npreds();
    const std::vector<std::pair<double,double> >& npred_weights = srcMap->npred_weights();

    // Subtract the contribution to the summed Npred spectrum.
    double xi(1);
    int kref(-1);
    size_t ne(m_dataCache.num_ebins());

    for (size_t k(0); k < ne; k++) {
      double npred_0 = npred_vals[k] * spec[k];
      double npred_1 = npred_vals[k+1] * spec[k+1];
      m_fixedNpreds[k] -= npred_0;
      m_fixedNpred_wts[k].first -= (npred_0 * npred_weights[k].first);
      m_fixedNpred_wts[k].second -= (npred_1 * npred_weights[k].second);

      if (use_edisp_val) {
	// FIXME, test this
	// xi = drm_cache->get_correction(k,kref,has_wts);
	xi = drm_cache->get_correction(k,kref,has_wts);
	if ( kref < 0 ) { 
	  // Still have true counts in this energy bin, just use k as kref
	  kref = k;
	}      
      } else {
	kref = k;
      }
      if ( kref != 0 ) {
	npred_0 = npred_vals[kref] * spec[kref];
	npred_1 = npred_vals[kref+1] * spec[kref+1];
      }
      m_fixedNpred_xis[k].first -= (npred_0*xi);
      m_fixedNpred_xis[k].second -= (npred_1*xi);
      m_fixedNpred_xiwts[k].first -= (npred_0*xi*npred_weights[kref].first);
      m_fixedNpred_xiwts[k].second -= (npred_1*xi*npred_weights[kref].first);
    }
    m_fixedNpreds[ne] -= (npred_vals[ne] * spec[ne]);
  
    // Subtract the source weights from the summed fixed model weights.
    bool subtract;
    addSourceWts(m_fixedModelWts, srcName, srcMap, subtract=true);
  
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
    std::vector<std::pair<double, double> > modelWts;
    std::pair<double, double> zeros(0, 0);
    modelWts.resize(m_dataCache.nFilled(), zeros);
    addSourceWts(modelWts, srcName);
    for (size_t j(0); j < m_dataCache.nFilled(); j++) {
      size_t k(m_dataCache.filledPixels()[j]/num_pixels());
      if (k < kmin || k > kmax-1) {
	continue;
      }
      double emin(m_dataCache.energies()[k]);
      double emax(m_dataCache.energies()[k+1]);
      double srcProb = pixelCounts(emin, emax, 
				   modelWts[j].first,
				   modelWts[j].second,
				   m_dataCache.log_energy_ratios()[k])/m_model[j];
      size_t indx = m_dataCache.filledPixels()[j];
      counts_spectrum[k - kmin] += srcProb*m_dataCache.data()[indx];
    }
    return counts_spectrum;
  }


  bool BinnedLikelihood::use_edisp(const std::string & srcname) const {
    bool use_src_edisp = m_config.use_edisp();
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
    return ( m_config.use_edisp() && use_src_edisp );
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


 



double BinnedLikelihood::pixelCounts(double emin, double emax, double y1, double y2, double log_ratio) const {
  return m_config.use_linear_quadrature() ? 
    FitUtils::pixelCounts_linearQuad(emin,emax,y1,y2) :
    FitUtils::pixelCounts_loglogQuad(emin,emax,y1,y2,log_ratio);
}


void BinnedLikelihood::computeFixedCountsSpectrum() {
  m_fixed_counts_spec.clear();
  for (size_t k(0); k < m_dataCache.num_ebins(); k++) {
    m_fixed_counts_spec.push_back(pixelCounts(m_dataCache.energies()[k], 
					      m_dataCache.energies()[k+1],
					      m_fixedNpred_xis[k].first,
					      m_fixedNpred_xis[k].second,
					      m_dataCache.log_energy_ratios()[k]));
  }
}
  

  void BinnedLikelihood::addSourceWts(std::vector<std::pair<double, double> > & modelWts,
				      const std::string & srcName,
				      SourceMap * srcMap,
				      bool subtract,
				      bool latchParams) const {
    const Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));
    if (src == 0) {
      return;
    }
    if (modelWts.size() != m_dataCache.nFilled()) {
      throw std::runtime_error("BinnedLikelihood::addSourceWts: "
			       "modelWts size does not match "
			       "number of filled pixels.");
    }
    SourceMap * sourceMap = srcMap;
    if ( sourceMap == 0 ) {
      sourceMap = getSourceMap(srcName);
    } 
    
    sourceMap->setSpectralValues(m_dataCache.energies(),latchParams);
    
    bool use_edisp_val = use_edisp(srcName);
    const Drm_Cache* drm_cache = sourceMap->drm_cache();
    
    SourceMapCache::addSourceWts_static(modelWts,*sourceMap,num_pixels(),
					m_dataCache.filledPixels(),drm_cache,use_edisp_val,subtract);
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

} // namespace Likelihood
