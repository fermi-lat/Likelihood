/**
 * @file BinnedLikelihood.cxx
 * @brief Photon events are binned in sky direction and energy.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedLikelihood.cxx,v 1.123 2016/09/13 19:26:23 echarles Exp $
 */

#include <cmath>
#include <cstdlib>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "st_stream/StreamFormatter.h"

#include "tip/IFileSvc.h"

#include "st_facilities/Util.h"
#include "st_facilities/Timer.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/FitUtils.h"
#include "Likelihood/FileUtils.h"
#include "Likelihood/WeightMap.h"

#include "Likelihood/Drm.h"
#define ST_DLL_EXPORTS
#include "Likelihood/SourceMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/SourceModel.h"

namespace Likelihood {


bool BinnedLikelihood::fileHasSourceMap(const std::string& srcName, 
					const std::string& fitsFile) {
  try {
    std::auto_ptr<const tip::Extension> ext(tip::IFileSvc::instance().readExtension(fitsFile,srcName));
  } catch (tip::TipException &) {
    return false;
  }
  return true;
}

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
    m_pixels(dataMap.pixels()),
    m_kmin(0), m_kmax(0), 
    m_modelIsCurrent(false),
    m_weightMap_orig(0),
    m_weightMap(0),
    m_weightedCounts(0),
    m_updateFixedWeights(true),
    m_drm(0),
    m_srcMapsFile(srcMapsFile),
    m_config(computePointSources,
	     applyPsfCorrections,
	     performConvolution,
	     resample,
	     resamp_factor,
	     minbinsz,
	     PsfIntegConfig::adaptive,1e-3,1e-6,
	     true,false,true,false) {
  dataMap.getEnergies(m_energies);
  log_energy_ratios(m_energies,m_log_energy_ratios);
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
    m_pixels(dataMap.pixels()),
    m_kmin(0), m_kmax(0), 
    m_modelIsCurrent(false),
    m_weightMap_orig(&weightMap),
    m_weightMap(0),
    m_weightedCounts(0),
    m_updateFixedWeights(true),
    m_drm(0),
    m_srcMapsFile(srcMapsFile),
    m_config(computePointSources,
	     applyPsfCorrections,
	     performConvolution,
	     resample,
	     resamp_factor,
	     minbinsz,	     
	     PsfIntegConfig::adaptive,1e-3,1e-6,
	     true,false,true,false){  
  
  dataMap.getEnergies(m_energies);
  log_energy_ratios(m_energies,m_log_energy_ratios);
  m_kmax = m_energies.size() - 1;
  
  if ( fileHasSourceMap("__weights__",m_srcMapsFile) ) {
    st_stream::StreamFormatter formatter("BinnedLikelihood","", 2);
    formatter.warn() << "Reading existing weights map from file " << m_srcMapsFile << std::endl;
    m_weightMap = new WeightMap(m_srcMapsFile, &dataMap, m_observation);
  } else {
    m_weightMap = new WeightMap(weightMap,&dataMap,observation,true);
  }
  identifyFilledPixels();
  m_fixedModelWts.resize(m_filledPixels.size(), std::make_pair(0, 0));
  m_fixedNpreds.resize(m_energies.size(), 0);
  computeCountsSpectrum();
}

BinnedLikelihood::~BinnedLikelihood() throw() {
  std::map<std::string, SourceMap *>::iterator srcMap(m_srcMaps.begin());
  for ( ; srcMap != m_srcMaps.end(); ++srcMap) {
    delete srcMap->second;
  }
  delete m_weightMap;
  delete m_drm;
}


void BinnedLikelihood::setCountsMap(const std::vector<float> & counts) {
  if(counts.size() != m_dataMap.data().size())
    throw std::runtime_error("Wrong size for input counts map.");
  m_dataMap.setImage(counts);
  identifyFilledPixels();
  buildFixedModelWts();
}


double BinnedLikelihood::value(optimizers::Arg & dummy) const {
  (void)(dummy);
  
  // Here we want the weighted verison of the nPred
  double npred = computeModelMap_internal(true);
  
  const std::vector<float> & data = m_weightMap == 0 ? m_dataMap.data() : m_weightedCounts->data();
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



void BinnedLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
  st_stream::StreamFormatter formatter("BinnedLikelihood", "getFreeDerivs",4);
  
  int nparams(getNumFreeParams());
  derivs.resize(nparams, 0);
  double npred;
  if (!m_modelIsCurrent) {
    // here we want the weighted version of npred
    npred = computeModelMap_internal(true);
  }
  const std::vector<float> & data = m_weightMap == 0 ? m_dataMap.data() : m_weightedCounts->data();
  
  /// First j value corresponding to the minimum allowed k-index.
  size_t jentry(0);
  for (size_t j(0); j < m_filledPixels.size(); j++) {
    size_t k(m_filledPixels[j]/m_pixels.size());
    if (k == m_kmin) {
      jentry = j;
      break;
    }
  }
  
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
      srcMap.setSpectralValues(m_energies);
      srcMap.setSpectralDerivs(m_energies,parnames);
    }
  }
  
  // /// For debugging/development, keep track of time spent in main loop.
  //    st_facilities::Timer timer;
  
  std::vector<Kahan_Accumulator> posDerivs(nparams);
  std::vector<Kahan_Accumulator> negDerivs(nparams);


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
     
      for (std::vector<Source *>::const_iterator it(free_srcs.begin());
	   it != free_srcs.end(); ++it ) {
	Source * src(*it);
	
	std::string srcName = src->getName();
	SourceMap & srcMap = sourceMap(srcName);
	
	const std::vector<float> & model = srcMap.model();
	const std::vector< std::vector<double> > & specDerivs = srcMap.cached_specDerivs();

	for (size_t i(0); i < specDerivs.size(); i++, iparam++) {
	  double my_deriv = pixelCounts(emin,emax,
					model.at(jmin)*specDerivs[i][k],
					model.at(jmax)*specDerivs[i][k+1],
					m_log_energy_ratios[k]);
	  double addend(data.at(jmin)/m_model.at(j)*my_deriv);
	  if (addend > 0) {
	    posDerivs[iparam].add(addend);
	  } else {
	    negDerivs[iparam].add(addend);
	  }
	}
      }
    }
  }

  size_t iparam2(0);
  for (std::vector<Source *>::const_iterator it2(free_srcs.begin());
       it2 != free_srcs.end(); ++it2 ) {
    Source * src(*it2);
    std::string srcName = src->getName();
    SourceMap & srcMap = sourceMap(srcName);

    const std::vector<double> & npreds =  srcMap.npreds();
    const std::vector<std::pair<double,double> > & npred_weights =  srcMap.npred_weights();
    const std::vector< std::vector<double> > & specDerivs = srcMap.cached_specDerivs();
   
    for (size_t i(0); i < specDerivs.size(); i++, iparam2++) {
      for (size_t kk(0); kk < m_energies.size()-1; kk++) {
	if (kk < m_kmin || kk > m_kmax-1) {
	  continue;
	}
	double addend = pixelCounts(m_energies.at(kk),m_energies.at(kk+1),
				    npreds.at(kk)*npred_weights[kk].first*specDerivs[i][kk],
				    npreds.at(kk+1)*npred_weights[kk].second*specDerivs[i][kk+1],
				    m_log_energy_ratios[kk]);
	if (-addend > 0) {
	  posDerivs[iparam2].add(-addend);
	} else {
	  negDerivs[iparam2].add(-addend);
	}
      }
    }
  }
  //    timer.stop();
  //    timer.report("main loop time");
  
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
    if ( m_config.use_single_fixed_map() && src->fixedSpectrum()) {
      addFixedSource(src->getName());
    } else {
      SourceMap * srcMap(getSourceMap(src->getName(), true));
      if (srcMap) {
	m_srcMaps[src->getName()] = srcMap;
      }
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

    sourceMap.setSpectralValues(m_energies);
    updateCorrectionFactors(srcName,sourceMap);
    double value = sourceMap.summed_counts(m_kmin,m_kmax,
					   use_edisp(srcName),weighted);
    return value;
  }

  bool BinnedLikelihood::hasSourceMap(const std::string & name) const {
    if (m_srcMaps.find(name) == m_srcMaps.end())
      return false;
    else
      return true;
  }

  SourceMap & BinnedLikelihood::sourceMap(const std::string & name) const {
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

  SourceMap * BinnedLikelihood::getSourceMap(const std::string & srcName, 
					     bool verbose) const {
  
    // Check to see if we already have the map
    std::map<std::string, SourceMap *>::iterator itrFind = m_srcMaps.find(srcName);
    if ( itrFind != m_srcMaps.end() ) return itrFind->second;
 
    // Check to see if the source map is in the file  
    Source * src = (const_cast<BinnedLikelihood *>(this))->getSource(srcName);

    if (fileHasSourceMap(srcName, m_srcMapsFile)) {
      return new SourceMap(m_srcMapsFile, *src, &m_dataMap, m_observation, m_weightMap);
    }
  
    // Generate the source map    
    if (src->getType() == "Diffuse" || m_config.computePointSources() ) {
      Drm* the_drm (0);
      if ( use_edisp(srcName) ) {
	BinnedLikelihood* nct = const_cast<BinnedLikelihood*>(this);
	the_drm = &(nct->drm());
      }
      return new SourceMap(*src, &m_dataMap, m_observation,
			   m_config.psf_integ_config(),
			   the_drm,m_weightMap);
    }
    return 0;
  }


  SourceMap * BinnedLikelihood::createSourceMap(const std::string & srcName) {
    Source * src = getSource(srcName);
    return new SourceMap(*src, &m_dataMap, m_observation, m_config.psf_integ_config(), m_drm, m_weightMap);
  }

  void BinnedLikelihood::eraseSourceMap(const std::string & srcName) {
    delete m_srcMaps[srcName];
    m_srcMaps.erase(srcName);
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
  
    if(!(src->getType() == "Diffuse" || m_config.computePointSources() ))
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
  

  void BinnedLikelihood::setSourceMapImage(const std::string & name,
					   const std::vector<float>& image) {
    // Create the source map if it doesn't exist
    if (m_srcMaps.find(name) == m_srcMaps.end()) {
      m_srcMaps[name] = getSourceMap(name);
    }
    m_srcMaps[name]->setImage(image);
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
  


  double BinnedLikelihood::computeModelMap_internal(bool weighted) const {
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
				  modelWts.at(j).second, m_log_energy_ratios[k]);
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
    size_t npix(m_pixels.size());
    modelMap.clear();
    modelMap.resize(npix*(m_energies.size()-1), 0);
    bool hasMap = hasSourceMap(srcName);
    SourceMap& srcMap = sourceMap(srcName);
    updateCorrectionFactors(srcName,srcMap);
    updateModelMap(modelMap, &srcMap);
    if( !hasMap ) {
      BinnedLikelihood* nct = const_cast<BinnedLikelihood*>(this);
      nct->eraseSourceMap(srcName);
    }
  }

  void BinnedLikelihood::computeModelMap(const std::vector<std::string>& srcNames, 
					 std::vector<float> & modelMap) const {
    size_t npix(m_pixels.size());
    modelMap.clear();
    modelMap.resize(npix*(m_energies.size()-1), 0);
  
    std::vector<std::string>::const_iterator name = srcNames.begin();
    for ( ; name != srcNames.end(); ++name) {  
      bool hasMap = hasSourceMap(*name);
      SourceMap& srcMap = sourceMap(*name);
      updateCorrectionFactors(*name,srcMap);
      updateModelMap(modelMap, &srcMap);
      if( !hasMap ) {
	BinnedLikelihood* nct = const_cast<BinnedLikelihood*>(this);
	nct->eraseSourceMap(*name);
      }
    }
  }


  void BinnedLikelihood::updateModelMap(std::vector<float> & modelMap,
					SourceMap * srcMap) const {
    size_t npix(m_pixels.size());
    std::string name(srcMap->name());
    NpredValue(name); // This computes the convolved spectrum.
    const Source * src = const_cast<BinnedLikelihood*>(this)->getSource(name);
    const std::vector<float> & model(srcMap->model());
    const std::vector<double> & specVals = srcMap->specVals();
  
    const WeightMap* mask = srcMap->weights();
  
    bool use_edisp_val = use_edisp(name);
    const Drm_Cache* drm_cache = srcMap->drm_cache();

    int kref(-1);
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
	if ( use_edisp_val ) {

	  double xi = drm_cache->get_correction(k,kref);
	  if ( kref < 0 ) {
	    // Correction factor is for this energy, use it
	    wt1 = specVals[k]*model[jmin]*xi;
	    wt2 = specVals[k+1]*model[jmax]*xi;
	  } else {
	    // Correction factor is for different energy bin, use kref
	    size_t ipix(jmin % npix);
	    size_t jref = kref*npix + ipix;
	    wt1 = specVals[k]*model[jref]*xi;
	    wt2 = specVals[k+1]*model[jref+npix]*xi;
	  }
	} else {
	  wt1 = specVals[k]*model[jmin];
	  wt2 = specVals[k+1]*model[jmax];
	}
	modelMap[jmin] += pixelCounts(emin, emax,wt1, wt2, m_log_energy_ratios[k]);
      }
    }
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


  const std::vector<double>& 
  BinnedLikelihood::modelCountsSpectrum(const std::string & srcname) const {
    if (!m_modelIsCurrent) {
      computeModelMap_internal();
    }

    NpredValue(srcname);
    SourceMap * srcMap = getSourceMap(srcname, false);
    return srcMap->drm_cache()->meas_counts();
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
	
	std::map<std::string, SourceMap *>::const_iterator itrFind = m_srcMaps.find(srcIt->first);
	if ( itrFind == m_srcMaps.end() ) {
	  throw std::runtime_error("BinnedLikelihood::fixedModelUpdated: "
				   "inconsistent m_srcMaps.");
	}
	if ( itrFind->second->spectrum_changed() ) {
	  return true;
	}	
      }
    }
    return false;
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
      }
    }
  
    computeFixedCountsSpectrum();
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
      m_srcMaps[srcName] = srcMap;
    } else {
      srcMap = srcMapIt->second;
    }
    addSourceWts(m_fixedModelWts, srcName, srcMap);
    double xi(1);

    bool use_edisp_val = use_edisp(srcName);
    const Drm_Cache* drm_cache = srcMap->drm_cache();
  
    int kref(-1);
    for (size_t k(0); k < m_energies.size()-1; k++) {
      optimizers::dArg ee(m_energies[k]);
      if (use_edisp_val) {
	xi = drm_cache->get_correction(k,kref);
	if ( kref >= 0 ) {
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
    srcMap->clear_model();
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
    if (srcMap == 0) {
      throw std::runtime_error("SourceMap cannot be created for " + srcName);
    }
    m_srcMaps[srcName] = srcMap;

    bool use_edisp_val = use_edisp(srcName);
    const Drm_Cache* drm_cache = srcMap->drm_cache();
  
    // Subtract the contribution to the summed Npred spectrum.
    double xi(1);
    int kref(-1);
    for (size_t k(0); k < m_energies.size(); k++) {
      optimizers::dArg ee(m_energies[k]);
      if (use_edisp_val) {
	xi = drm_cache->get_correction(k,kref);
	if ( kref >= 0 ) {
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
      npred = computeModelMap_internal();
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
      double srcProb = pixelCounts(emin, emax, 
				   modelWts[j].first,
				   modelWts[j].second,
				   m_log_energy_ratios[k])/m_model[j];
      size_t indx = m_filledPixels[j];
      counts_spectrum[k - kmin] += srcProb*m_dataMap.data()[indx];
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
      m_drm = new Drm(m_dataMap.refDir().ra(), m_dataMap.refDir().dec(), 
		      observation(), m_dataMap.energies());
    }
    return *m_drm;
  }

  void BinnedLikelihood::saveWeightsMap(bool replace) const {
    if ( m_weightMap == 0 ) return;
    bool has_weights = fileHasSourceMap("__weights__", m_srcMapsFile);    
    if ( has_weights ) {
      if ( replace ) {
	FileUtils::replace_image_from_float_vector(m_srcMapsFile,"__weights__",
						   m_dataMap,m_weightMap->model(),false);
      } else {
	// just leave it be;
	;
      }
    } else {
      	FileUtils::append_image_from_float_vector(m_srcMapsFile,"__weights__",
						  m_dataMap,m_weightMap->model(),false);	
    }
  }


  void BinnedLikelihood::fillWeightedCounts() {
    delete m_weightedCounts;
    m_weightedCounts = m_dataMap.clone();
    // FIXME, this would be more efficient if CountsMapBase had a multiply by function
    size_t ne = m_energies.size();
    size_t npix = m_pixels.size();
    std::vector<float> wts(ne*npix);
    for ( size_t j(0); j < npix; j++ ) {
      for (size_t k(0); k < ne-1; k++) {
	size_t idx = k*npix +j;
	double w = m_weightMap->model()[idx];
	bool is_null = w <= 0 || m_dataMap.data()[idx] <= 0;
	w = is_null ? 0 : w;
	wts[idx] = m_dataMap.data()[idx] * w;
      }
    }
    m_weightedCounts->setImage(wts);
  }

  double BinnedLikelihood::spectrum(const Source * src, double energy)  {
    const optimizers::Function & spectrum(src->spectrum());
    optimizers::dArg eArg(energy);
    return spectrum(eArg);
  }


  void BinnedLikelihood::log_energy_ratios(const std::vector<double>& energies,
					   std::vector<double>& log_ratios) {
    log_ratios.resize(energies.size() -1);
    for ( size_t i(0); i < log_ratios.size(); i++ ) {
      log_ratios[i] = std::log(energies[i+1] / energies[i]);
    }
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


  void BinnedLikelihood::addSourceWts_static(std::vector<std::pair<double, double> > & modelWts,
					     SourceMap& srcMap,
					     size_t npix,
					     const std::vector<unsigned int>& filledPixels,
					     const Drm_Cache* drm_cache,
					     bool use_edisp_val,
					     bool subtract) {
    double my_sign(1.);
    if (subtract) {
      my_sign = -1.;
    }
    int kref(-1);
    const std::vector<double> & spec = srcMap.specVals();
    const std::vector<float> & model = srcMap.model();
    for (size_t j(0); j < filledPixels.size(); j++) {
      size_t jmin(filledPixels.at(j));
      size_t jmax(jmin + npix);
      size_t k(jmin/npix);
      if (use_edisp_val) {
	double xi = drm_cache->get_correction(k,kref);
	if ( kref < 0 ) {
	  modelWts[j].first += my_sign*model[jmin]*spec[k]*xi;
	  modelWts[j].second += my_sign*model[jmax]*spec[k+1]*xi;
	} else {
	  size_t ipix(jmin % npix);	
	  size_t jref = kref*npix + ipix;
	  modelWts[j].first += (my_sign*model[jref]*spec[kref]*xi);
	  modelWts[j].second += (my_sign*model[jref+npix]*spec[kref+1]*xi);
	}
      } else {
	modelWts[j].first += my_sign*model[jmin]*spec[k];
	modelWts[j].second += my_sign*model[jmax]*spec[k+1];
      }
    }  
  }



void BinnedLikelihood::replaceSourceMap(const std::string & srcName,
                                        const std::string & fitsFile) const {

  SourceMap & srcMap = sourceMap(srcName);
  srcMap.setFilename(fitsFile);
  FileUtils::replace_image_from_float_vector(fitsFile,srcName,m_dataMap,
					     srcMap.model(),true);  
}

void BinnedLikelihood::appendSourceMap(const std::string & srcName,
                                       const std::string & fitsFile) const {

  SourceMap & srcMap = sourceMap(srcName);
  srcMap.setFilename(fitsFile);
  FileUtils::append_image_from_float_vector(fitsFile,srcName,m_dataMap,
					    srcMap.model(),true);
}


double BinnedLikelihood::pixelCounts(double emin, double emax, double y1, double y2, double log_ratio) const {
  return m_config.use_linear_quadrature() ? 
    FitUtils::pixelCounts_linearQuad(emin,emax,y1,y2) :
    FitUtils::pixelCounts_loglogQuad(emin,emax,y1,y2,log_ratio);
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


void BinnedLikelihood::computeFixedCountsSpectrum() {
  m_fixed_counts_spec.clear();
  for (size_t k(0); k < m_energies.size()-1; k++) {
    m_fixed_counts_spec.push_back(pixelCounts(m_energies[k], 
					      m_energies[k+1],
					      m_fixedNpreds[k],
					      m_fixedNpreds[k+1],
					      m_log_energy_ratios[k]));
  }
}
  

void BinnedLikelihood::addSourceWts(std::vector<std::pair<double, double> > & modelWts,
				    const std::string & srcName,
				    SourceMap * srcMap,
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
  SourceMap * sourceMap(srcMap);
  if ( sourceMap == 0 ) {
    std::map<std::string, SourceMap *>::iterator srcIt 
      = m_srcMaps.find(srcName);
    if (srcIt == m_srcMaps.end()) {
      throw std::runtime_error("BinnedLikelihood::addSourceWts: "
			       "source not found in m_srcMaps.");
    }
    sourceMap = srcIt->second;
  } 
  
  sourceMap->setSpectralValues(m_energies);

  bool use_edisp_val = use_edisp(srcName);
  const Drm_Cache* drm_cache = sourceMap->drm_cache();
  
  addSourceWts_static(modelWts,*sourceMap,m_pixels.size(),
		      m_filledPixels,drm_cache,use_edisp_val,subtract);
}
  
  
void BinnedLikelihood::identifyFilledPixels() {
  if ( m_weightMap ) {
    fillWeightedCounts();
  }
  const std::vector<float> & data = m_weightMap == 0 ? m_dataMap.data() : m_weightedCounts->data();
  m_filledPixels.clear();
  for (unsigned int i = 0; i < data.size(); i++) {
    if (data.at(i) > 0) {
      m_filledPixels.push_back(i);
    }
  }
}
  

void BinnedLikelihood::updateCorrectionFactors(const std::string & srcName, 
					       SourceMap & sourceMap) const {

  const Source * src(const_cast<BinnedLikelihood *>(this)->getSource(srcName));

  Drm_Cache* drm_cache = const_cast<Drm_Cache*>(sourceMap.drm_cache());
  if ( drm_cache != 0 ) {
    Drm* the_drm(0);
    if ( use_edisp(srcName) ) {
      BinnedLikelihood* nct = const_cast<BinnedLikelihood*>(this);
      the_drm = &(nct->drm());
    }
    drm_cache->update(the_drm,sourceMap,energies());
  }
}

} // namespace Likelihood
