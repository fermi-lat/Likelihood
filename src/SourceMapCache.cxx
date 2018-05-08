/**
 * @file BinnedLikelihood.cxx
 * @brief Functionality to deal with source maps extracted from BinnedLikelihood
 * @author E. Charles, (from BinnedLikelihood by J. Chiang)
 *
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMapCache.cxx,v 1.10 2017/09/29 01:38:03 echarles Exp $
 */


#include "Likelihood/SourceMapCache.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "tip/Extension.h"
#include "st_stream/StreamFormatter.h"

#include "Likelihood/BinnedCountsCache.h"
#include "Likelihood/BinnedConfig.h"
#include "Likelihood/CountsMapBase.h"
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/Drm.h"
#include "Likelihood/FitUtils.h"
#include "Likelihood/FileUtils.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/WeightMap.h"
#include "Likelihood/PSFUtils.h"

#define ST_DLL_EXPORTS
#include "Likelihood/SourceMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/SourceModel.h"

namespace Likelihood {

  SourceMapCache::SourceMapCache(const BinnedCountsCache& dataCache,
				 const Observation & observation,
				 const std::string & srcMapsFile,
				 const BinnedLikeConfig& config,
				 const Drm* drm)
    : m_dataCache(dataCache),
      m_observation(observation),
      m_drm(drm),
      m_srcMapsFile(srcMapsFile),
      m_config(config) {
  }

  SourceMapCache::SourceMapCache(const SourceMapCache& other)
    : m_dataCache(other.m_dataCache),
      m_observation(other.m_observation),
      m_drm(other.m_drm),
      m_srcMapsFile(m_srcMapsFile),
      m_config(m_config) {
  }

  SourceMapCache::~SourceMapCache() {
    std::map<std::string, SourceMap *>::iterator srcMap(m_srcMaps.begin());
    for ( ; srcMap != m_srcMaps.end(); ++srcMap) {
      delete srcMap->second;
    }
  }

  CountsMapBase * SourceMapCache::createCountsMap(const std::vector<const Source*>& srcs, 
						  CountsMapBase & dataMap) const {
    std::vector<float> map;
    computeModelMap(srcs,map);
    dataMap.setImage(map);
    return &dataMap;
  } 

  CountsMapBase * SourceMapCache::createCountsMap(const std::vector<const Source*>& srcs) const {
    std::vector<float> map;
    computeModelMap(srcs,map);
    CountsMapBase * modelMap = m_dataCache.countsMap().clone();         
    modelMap->setImage(map);
    return modelMap;
  }

  double SourceMapCache::NpredValue(const Source& src, size_t kmin, size_t kmax, bool weighted) const {
    // call the other version of NpredValue
    SourceMap* srcMap = getSourceMap(src,false);
    return NpredValue(src, *srcMap, kmin, kmax, weighted);
  }

  // This version forces the recalculation of Npred, whether the source is
  // fixed or not. It is also called from buildFixedModelWts.
  double SourceMapCache::NpredValue(const Source& src,
				    SourceMap & sourceMap,
				    size_t kmin, size_t kmax,
				    bool weighted) const {

    sourceMap.setSpectralValues(m_dataCache.energies());
    updateCorrectionFactors(src,sourceMap);
    double value = sourceMap.summed_counts(kmin,kmax,
					   use_edisp(&src),weighted);
    return value;
  }

  bool SourceMapCache::hasSourceMap(const std::string & name) const {
    if (m_srcMaps.find(name) == m_srcMaps.end())
      return false;
    else
      return true;
  }

  SourceMap * SourceMapCache::getSourceMap(const Source& src,
					   bool verbose,
					   const BinnedLikeConfig* config) const {

    const std::string& srcName = src.getName();
    SourceMap* srcMap(0);
    const Drm* the_drm = use_edisp(&src) ? m_drm : 0;
 
    // Check to see if we already have the map
    std::map<std::string, SourceMap *>::iterator itrFind = m_srcMaps.find(srcName);
    if ( itrFind != m_srcMaps.end() ) {
      srcMap = itrFind->second;
      srcMap->setSource(src);
      srcMap->update_drm_cache(the_drm);
    } else {
      // Check to see if the map is in the file
      if ( m_config.load_existing_srcmaps() && FileUtils::fileHasExtension(m_srcMapsFile, srcName)) {
	srcMap = new SourceMap(m_srcMapsFile, src, &m_dataCache, 
			       m_observation, m_dataCache.weightMap(), the_drm, m_config.save_all_srcmaps());
      } else {
	switch ( src.srcType() ) {
	case Source::Point:
	  if  ( m_config.computePointSources() ) {
	    srcMap = createSourceMap(src, config);
	  }
	  break;
	case Source::Diffuse:
	case Source::Composite:
	  srcMap = createSourceMap(src, config);
	  break;
	default:
	  throw std::runtime_error("SourceMapCache::getSourceMap unknown source type for " + srcName);
	}
      }
      m_srcMaps[srcName] = srcMap;
    }
    return srcMap;
  }


  SourceMap * SourceMapCache::createSourceMap(const Source& src, const BinnedLikeConfig* config) const {
    const BinnedLikeConfig& the_config = config == 0 ? m_config : *config;
    const Drm* the_drm = use_edisp(&src) ? m_drm : 0;
    return new SourceMap(src, &m_dataCache, m_observation, the_config.psf_integ_config(), 
			 the_drm, m_dataCache.weightMap(), the_config.save_all_srcmaps() );
  }
  
  void SourceMapCache::eraseSourceMap(const std::string & srcName) {
    delete m_srcMaps[srcName];
    m_srcMaps.erase(srcName);
  }


  void SourceMapCache::insertSourceMap(const std::string & srcName,
				       SourceMap& srcMap) {
    std::map<std::string, SourceMap *>::iterator itr = m_srcMaps.find(srcName);
    if ( itr == m_srcMaps.end() ) {
      m_srcMaps[srcName] = &srcMap;
    } else {
      if ( itr->second != &srcMap ) {
	throw std::runtime_error("SourceMapCache already has a Source " + srcName + " in cache");
      }
    }
  }

  SourceMap* SourceMapCache::removeSourceMap(const std::string & srcName) {
    if ( ! hasSourceMap(srcName) ) {
      throw std::runtime_error("SourceMapCache does not have " + srcName + " in cache");
    }
    SourceMap* srcMap = m_srcMaps[srcName];
    m_srcMaps.erase(srcName);
    return srcMap;
  }


  void SourceMapCache::loadSourceMaps(const  std::vector<const Source*>& srcs,
				      bool recreate, bool saveMaps) {
    
 
    std::vector<tip::Extension*> hdus;

    for ( std::vector<const Source*>::const_iterator itr = srcs.begin();
	  itr != srcs.end(); itr++ ) {
    
      const Source* src = *itr;
      const std::string& name = src->getName();
      if (m_srcMaps.find(name) == m_srcMaps.end() || recreate) {
	loadSourceMap(*src,recreate);
      }
    
      if(saveMaps) {
	tip::Extension* ptr(0);
	if (FileUtils::fileHasExtension(m_srcMapsFile,name)) {
	  ptr = replaceSourceMap(*src, m_srcMapsFile);
	  delete ptr;
	  // hdus.push_back(ptr);
	} else {
	  ptr = appendSourceMap(*src, m_srcMapsFile);
	  delete ptr;
	  // hdus.push_back(ptr);
	}
      }
    }

    for ( std::vector<tip::Extension*>::iterator itrDel = hdus.begin();
	  itrDel != hdus.end(); itrDel++ ) {
      delete *itrDel;
    }
    
  }


  void SourceMapCache::loadSourceMap(const Source& src, bool recreate, const BinnedLikeConfig* config) {
    
    const std::string& srcName = src.getName();  
    if(!(src.getType() == "Diffuse" || m_config.computePointSources() ))
      return;
  
    std::map<std::string, SourceMap *>::iterator mapIt = m_srcMaps.find(srcName);
  
    //if( mapIt != m_srcMaps.end() ) {
    //  delete m_srcMaps[srcName];
    //  m_srcMaps.erase(srcName);
    //}
  
    SourceMap * srcMap = mapIt != m_srcMaps.end() ? mapIt->second : 0;
  
    if(recreate) {
      delete srcMap;
      srcMap = createSourceMap(src, config);
      m_srcMaps[srcName] = srcMap;    
    } else {
      srcMap = getSourceMap(src, true, config);      
    }
      
  }
  

  void SourceMapCache::setSourceMapImage(const Source& src,
					 const std::vector<float>& image) {
    // Create the source map if it doesn't exist
    const std::string& name = src.getName();
    if (m_srcMaps.find(name) == m_srcMaps.end()) {
      m_srcMaps[name] = getSourceMap(src);
    }
    m_srcMaps[name]->setImage(image);
  }



  void SourceMapCache::saveSourceMaps(const std::string & filename,
				      const std::vector<const Source*>& srcs,
				      bool replace) {
    if (filename != "") {
      m_srcMapsFile = filename;
    }

    std::vector<tip::Extension*> hdus;
    st_stream::StreamFormatter formatter("SourceMapCache",
					 "saveSourceMaps", 4);

    unsigned int i(0);
    unsigned int ifile(0);
    std::string outFileName = m_srcMapsFile;
    for ( std::vector<const Source*>::const_iterator itr = srcs.begin();
	  itr != srcs.end(); itr++, i++ ){
      const Source* src = *itr;
      tip::Extension* ptr(0);  

      if ( srcs.size() > 499 ) {
	if ( i % 500 == 0 ) {
	  if ( i > 0 ) {
	    ifile += 1;
	  }
	  std::ostringstream ofile;
	  ofile << m_srcMapsFile << '_' << ifile << ".fits";
	  outFileName = ofile.str();
	}
      }	      
      if (m_srcMaps.count(src->getName())) {
	if (FileUtils::fileHasExtension(outFileName, src->getName()) ) {
	    if ( replace ) {
	      ptr = replaceSourceMap(*src, outFileName);
	      delete ptr;
	      //hdus.push_back(ptr);
	    } else {
	    }
	} else {
	  formatter.info() << "appending map for " 
			   << src->getName() << std::endl;
	  ptr = appendSourceMap(*src, outFileName);
	  delete ptr;
	  // hdus.push_back(ptr);
	}
      }
      
    }
      
    for ( std::vector<tip::Extension*>::iterator itrDel = hdus.begin();
	  itrDel != hdus.end(); itrDel++ ) {
      delete *itrDel;
    }     
  }


  void SourceMapCache::saveSourceMap_partial(const std::string & filename,
					      const Source& source,
					      int kmin, int kmax,
					      bool replace) {
    st_stream::StreamFormatter formatter("SourceMapCache",
					 "saveSourceMaps_partial", 4);
    tip::Extension* ptr(0);            
    if (FileUtils::fileHasExtension(filename, source.getName()) ) {
      if ( replace ) {	  
	ptr = replaceSourceMap_partial(source, filename, kmin, kmax);
	delete ptr;
      } 
    } else {
      formatter.warn() << "appending map for " 
		       << source.getName() << std::endl;
      ptr = appendSourceMap_partial(source, filename, kmin, kmax);
      delete ptr;
    }   
  }

  
  void SourceMapCache::computeModelMap(const Source & src, 
				       std::vector<float> & modelMap,
				       bool use_mask) const {
    modelMap.resize(m_dataCache.data_map_size(), 0);      
    bool hasMap = hasSourceMap(src.getName());
    SourceMap* srcMap = getSourceMap(src);
    updateCorrectionFactors(src,*srcMap);
    updateModelMap(modelMap, src, srcMap, use_mask);
    if( !hasMap ) {
      delete srcMap;
    }
  }

  void SourceMapCache::computeModelMap(const std::vector<const Source*>& srcs, 
				       std::vector<float> & modelMap,
				       bool use_mask) const {
    modelMap.resize(m_dataCache.data_map_size(), 0);  
    for ( std::vector<const Source*>::const_iterator itr = srcs.begin();
	  itr != srcs.end(); itr++ ) {
      const Source* src = *itr;
      bool hasMap = hasSourceMap(src->getName());
      SourceMap* srcMap = getSourceMap(*src);
      updateCorrectionFactors(*src, *srcMap);
      updateModelMap(modelMap, *src, srcMap, use_mask);
      if( !hasMap && ! m_config.save_all_srcmaps() ) {
	delete srcMap;
      }
    }
  }


  void SourceMapCache::updateModelMap(std::vector<float> & modelMap,
				      const Source& src, 
				      SourceMap * srcMap,
				      bool use_mask) const {
    size_t npix = m_dataCache.num_pixels();
    size_t kmin = 0;
    size_t kmax = m_dataCache.num_ebins();
    const std::string& name = srcMap->name();

    double np = NpredValue(src,kmin,kmax); // This computes the convolved spectrum

    const std::vector<double> & specVals = srcMap->specVals();  
    const WeightMap* mask = use_mask ? srcMap->weights() : 0;
  
    bool use_edisp_val = use_edisp(&src);
    const Drm_Cache* drm_cache = srcMap->drm_cache();

    int kref(-1);
    for (size_t j(0); j < npix; j++) {
      for (size_t k(0); k < m_dataCache.num_ebins(); k++) {
	double emin(m_dataCache.energies().at(k));
	double emax(m_dataCache.energies().at(k+1));
	size_t jmin(k*npix + j);
	size_t jmax(jmin + npix);
	// EAC, skip masked pixels
	if ( mask && 
	     ( mask->model().at(jmin) <= 0. ) ) continue;
	double wt1(0);
	double wt2(0);	
	double xi(1.);
	if ( use_edisp_val ) {
	  xi = drm_cache->get_correction(k,kref);
	  if ( kref < 0 ) {
	    // Correction factor is for this energy, use it
	    wt1 = specVals[k]*(*srcMap)[jmin];
	    wt2 = specVals[k+1]*(*srcMap)[jmax];
	  } else {
	    // Correction factor is for different energy bin, use kref
	    size_t ipix(jmin % npix);
	    size_t jref = kref*npix + ipix;
	    wt1 = specVals[k]*(*srcMap)[jref];
	    wt2 = specVals[k+1]*(*srcMap)[jref+npix];
	  }
	} else {
	  wt1 = specVals[k]*(*srcMap)[jmin];
	  wt2 = specVals[k+1]*(*srcMap)[jmax];
	}
	modelMap[jmin] += xi*FitUtils::pixelCounts_loglogQuad(emin, emax, wt1, wt2, m_dataCache.log_energy_ratios()[k]);
      }
    }
  }


  void SourceMapCache::getNpreds(const Source& src,
				 std::vector<double> & npreds) const {
    
    SourceMap * srcMap = getSourceMap(src, false);
    npreds = srcMap->npreds();    
  }


  const std::vector<double>& SourceMapCache::modelCountsSpectrum(const Source& src,
								 bool weighted) const {
    SourceMap * srcMap = getSourceMap(src, false);
    // This forces updating the source map if needed
    const Drm_Cache* drm_cache = srcMap->drm_cache(true);
    return weighted ? drm_cache->meas_counts_wt() : drm_cache->meas_counts();
  }


  void SourceMapCache::setWeightsMap(const ProjMap* wmap) {
    BinnedCountsCache& nc_dataCache = const_cast<BinnedCountsCache&>(m_dataCache);
    nc_dataCache.setWeightsMap(wmap, m_observation);
    const WeightMap* wwmap = nc_dataCache.weightMap();
    for ( std::map<std::string, SourceMap *>::iterator itr = m_srcMaps.begin();
	  itr != m_srcMaps.end(); itr++ ) {
      SourceMap* smap = itr->second;
      smap->setWeights(wwmap);
    }
  }

  void SourceMapCache::fillSummedSourceMap(const std::vector<const Source*>& sources, 
					   std::vector<float>& model, 
					   int kmin, int kmax) {
    for ( std::vector<const Source*>::const_iterator srcIt = sources.begin();  
	  srcIt != sources.end(); ++srcIt) {
      SourceMap* srcMap = getSourceMap(*(*srcIt), false);
      srcMap->addToVector(model,true,kmin,kmax);
    }
  }
  
  void SourceMapCache::fillSingleSourceMap(const Source& src,
					   std::vector<float>& model, 
					   FileUtils::SrcMapType& mapType,
					   int kmin, int kmax) const {
    st_stream::StreamFormatter formatter("fillSingleSourceMap", "", 2);

    const MeanPsf* meanPsf(0);
    if ( src.srcType() == Source::Point && ( ! m_config.psf_integ_config().use_single_psf() ) ) {
      meanPsf = PSFUtils::build_psf(src, m_dataCache.countsMap(),m_observation);
    }
      
    PSFUtils::makeModelMap(src, m_dataCache,
			   meanPsf==0 ? m_observation.meanpsf() : *meanPsf,
			   m_observation.bexpmap_ptr(),
			   m_config.psf_integ_config(),
			   m_srcMapsFile,m_drm,
			   formatter,
			   model,mapType,kmin,kmax);
    delete meanPsf;
  }

  bool SourceMapCache::use_edisp(const Source* src) const {
    bool retVal = m_config.use_edisp();
    retVal &=  src != 0 ? src->use_edisp() : true;
    return retVal;
  }

  size_t SourceMapCache::memory_size() const {
    size_t sum(0);
    for ( std::map<std::string, SourceMap *>::const_iterator itr = m_srcMaps.begin();
	  itr != m_srcMaps.end(); itr++ ) {
      sum += itr->second->memory_size();
    }
    return sum;
  }
  

  void SourceMapCache::addSourceWts_static(std::vector<std::pair<double, double> > & modelWts,
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
    for (size_t j(0); j < filledPixels.size(); j++) {
      size_t jmin(filledPixels.at(j));
      size_t jmax(jmin + npix);
      size_t k(jmin/npix);
      if (use_edisp_val) {
	double xi = drm_cache->get_correction(k,kref);
	if ( kref < 0 ) {
	  modelWts[j].first += my_sign*srcMap[jmin]*spec[k]*xi;
	  modelWts[j].second += my_sign*srcMap[jmax]*spec[k+1]*xi;
	} else {
	  size_t ipix(jmin % npix);	
	  size_t jref = kref*npix + ipix;
	  modelWts[j].first += (my_sign*srcMap[jref]*spec[kref]*xi);
	  modelWts[j].second += (my_sign*srcMap[jref+npix]*spec[kref+1]*xi);
	}
      } else {
	modelWts[j].first += my_sign*srcMap[jmin]*spec[k];
	modelWts[j].second += my_sign*srcMap[jmax]*spec[k+1];
      }
    }  
  }

  tip::Extension* SourceMapCache::replaceSourceMap(const Source & src,
						   const std::string & fitsFile) const {
    
    SourceMap* srcMap = getSourceMap(src,false);
    if ( srcMap->cached_model().size() == 0 &&
	 srcMap->cached_sparse_model().size() == 0 ) {
      srcMap->model();
    }
    srcMap->setFilename(fitsFile);    
    srcMap->setModelIsLocal(false);
    switch ( srcMap->mapType() ) {
    case FileUtils::HPX_Sparse:
      return FileUtils::replace_image_from_sparse_vector_healpix(fitsFile,src.getName(),
								 static_cast<const CountsMapHealpix&>(m_dataCache.countsMap()),
								 srcMap->cached_sparse_model(),true);
      break;
    case FileUtils::WCS:
    case FileUtils::HPX_AllSky:
    case FileUtils::HPX_Partial:
    default:
      return FileUtils::replace_image_from_float_vector(fitsFile,src.getName(),m_dataCache.countsMap(),
							srcMap->cached_model(),true);
    }
  }

  tip::Extension* SourceMapCache::appendSourceMap(const Source & src,
						  const std::string & fitsFile) const {
    
    SourceMap* srcMap = getSourceMap(src,false);
    if ( srcMap->cached_model().size() == 0 &&
	 srcMap->cached_sparse_model().size() == 0 ) {
      srcMap->model();
    }
    srcMap->setFilename(fitsFile);
    srcMap->setModelIsLocal(false);
    switch ( srcMap->mapType() ) {
    case FileUtils::HPX_Sparse:
      return FileUtils::append_image_from_sparse_vector_healpix(fitsFile,src.getName(),
								static_cast<const CountsMapHealpix&>(m_dataCache.countsMap()),
								srcMap->cached_sparse_model(),true);
    case FileUtils::WCS:
    case FileUtils::HPX_AllSky:
    case FileUtils::HPX_Partial:
    default:
      return FileUtils::append_image_from_float_vector(fitsFile,src.getName(),m_dataCache.countsMap(),
						       srcMap->cached_model(),true);
    }
  }
  

  

  tip::Extension* SourceMapCache::replaceSourceMap_partial(const Source & src,
							   const std::string & fitsFile,
							   int kmin, int kmax) const {
    
    std::vector<float> model;
    FileUtils::SrcMapType mapType = FileUtils::Unknown;
    fillSingleSourceMap(src, model, mapType, kmin, kmax);
  
    float sum(0);
    FitUtils::sumVector(model.begin(),model.end(),sum);

    SparseVector<float> sparse;
    switch ( mapType ) {
    case FileUtils::HPX_Sparse:
      sparse.resize(model.size());
      sparse.fill_from_vect(model);
      return FileUtils::replace_image_from_sparse_vector_healpix(fitsFile,src.getName(),
								 static_cast<const CountsMapHealpix&>(m_dataCache.countsMap()),
								 sparse,true,
								 kmin,kmax);
      break;
    case FileUtils::WCS:
    case FileUtils::HPX_AllSky:
    case FileUtils::HPX_Partial:
    default:
      return FileUtils::replace_image_from_float_vector(fitsFile,src.getName(),m_dataCache.countsMap(),
							model,true,kmin,kmax);
    }
  }

  tip::Extension* SourceMapCache::appendSourceMap_partial(const Source & src,
							  const std::string & fitsFile,
							  int kmin, int kmax) const {
							    
    std::vector<float> model;
    FileUtils::SrcMapType mapType = FileUtils::Unknown;
    fillSingleSourceMap(src, model, mapType, kmin, kmax);

    SparseVector<float> sparse;

    switch ( mapType ) {
    case FileUtils::HPX_Sparse:
      sparse.resize(model.size());
      sparse.fill_from_vect(model);
      return FileUtils::append_image_from_sparse_vector_healpix(fitsFile,src.getName(),
								static_cast<const CountsMapHealpix&>(m_dataCache.countsMap()),
								sparse,true,kmin,kmax);
    case FileUtils::WCS:
    case FileUtils::HPX_AllSky:
    case FileUtils::HPX_Partial:
    default:
      return FileUtils::append_image_from_float_vector(fitsFile,src.getName(),m_dataCache.countsMap(),
						       model,true,kmin,kmax);
    }
  }

  void SourceMapCache::addSourceWts(std::vector<std::pair<double, double> > & modelWts,
				    const Source & src,
				    SourceMap * srcMap,
				    bool subtract,
				    bool latchParams) const {
    if (modelWts.size() != m_dataCache.nFilled() ) {
      throw std::runtime_error("SourceMapCache::addSourceWts: "
			       "modelWts size does not match "
			       "number of filled pixels.");
    }
    SourceMap * sourceMap = getSourceMap(src);
    sourceMap->setSpectralValues(m_dataCache.energies(),latchParams);
    
    bool use_edisp_val = use_edisp(&src);
    const Drm_Cache* drm_cache = sourceMap->drm_cache();
    
    addSourceWts_static(modelWts,*sourceMap,m_dataCache.num_pixels(),
			m_dataCache.filledPixels(),drm_cache,use_edisp_val,subtract);
  }
  
  
  
  void SourceMapCache::updateCorrectionFactors(const Source & src,
					       SourceMap & sourceMap) const {
    
    Drm_Cache* drm_cache = const_cast<Drm_Cache*>(sourceMap.drm_cache());
    if ( m_drm == 0 ) {
      throw std::runtime_error("No DRM object");
    }
    const Drm* the_drm = use_edisp(&src) ? m_drm : 0;
    sourceMap.update_drm_cache(the_drm, true);
  }
  
} // namespace Likelihood
