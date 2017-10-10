/**
 * @file BinnedCountsCache.cxx
 * @brief  
 * @author 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedCountsCache.cxx,v 1.3 2017/10/06 01:31:00 echarles Exp $
 */

#include "Likelihood/BinnedCountsCache.h"

#include <stdexcept>

#include "Likelihood/WeightMap.h"
#include "Likelihood/FileUtils.h"
#include "Likelihood/ProjMap.h"

namespace Likelihood {

  BinnedCountsCache::BinnedCountsCache(CountsMapBase & dataMap,
				       const Observation & observation,
				       const ProjMap* weightMap,
				       const std::string& srcMapsFile,
				       bool overwriteWeights)
    :m_dataMap(dataMap), 
     m_numPixels(dataMap.pixels().size()),
     m_weightMap_orig(weightMap),
     m_weightMap(0),
     m_weightedCounts(0) {
    
    log_energy_ratios(m_dataMap.energies(),m_log_energy_ratios);
  
    if ( weightMap != 0 ) {
      if ( FileUtils::fileHasExtension(srcMapsFile,"__weights__") && !overwriteWeights ) {
	st_stream::StreamFormatter formatter("BinnedLikelihood","", 2);
	formatter.warn() << "Reading existing weights map from file " << srcMapsFile << std::endl;
	m_weightMap = new WeightMap(srcMapsFile, &dataMap, observation);
      } else {
	m_weightMap = new WeightMap(*weightMap,&dataMap,observation,true);
      }
    }
    identifyFilledPixels();  
    computeCountsSpectrum();
  }

  BinnedCountsCache::~BinnedCountsCache() {
    delete m_weightMap;
    delete m_weightedCounts;
  }


  void BinnedCountsCache::setCountsMap(const std::vector<float> & counts) {
    if(counts.size() != m_dataMap.data().size())
      throw std::runtime_error("Wrong size for input counts map.");
    m_dataMap.setImage(counts);
    identifyFilledPixels();
    computeCountsSpectrum();
  }

  void BinnedCountsCache::setWeightsMap(const ProjMap* wmap, 
					const Observation & observation) {
    if ( wmap == m_weightMap_orig) return;
    delete m_weightMap_orig;
    delete m_weightMap;
    m_weightMap_orig = wmap;
    if ( wmap == 0 ) {
      m_weightMap = 0;
    } else {
      m_weightMap = new WeightMap(*wmap,&m_dataMap,observation,true);
    }
    identifyFilledPixels();
    computeCountsSpectrum();
  }
  

  tip::Extension* BinnedCountsCache::saveWeightsMap(const std::string& srcMapsFile, bool replace) const {
    if ( m_weightMap == 0 ) return 0;
    bool has_weights = FileUtils::fileHasExtension(srcMapsFile,"__weights__");    
    if ( has_weights ) {
      if ( replace ) {
	return FileUtils::replace_image_from_float_vector(srcMapsFile,"__weights__",
							  m_dataMap,m_weightMap->model(),false);
      } else {
	// just leave it be;
	;
      }
    } else {
      return FileUtils::append_image_from_float_vector(srcMapsFile,"__weights__",
						       m_dataMap,m_weightMap->model(),false);	
    }
    return 0;
  }

  void BinnedCountsCache::fillWeightedCounts() {
    if ( m_weightMap == 0 ) {
      delete m_weightedCounts;
      m_weightedCounts = 0;
      return;
    }
    delete m_weightedCounts;
    m_weightedCounts = m_dataMap.clone();
    // FIXME, this would be more efficient if CountsMapBase had a multiply by function
    size_t ne = num_ebins();
    size_t npix = num_pixels();
    std::vector<float> wts(ne*npix);
    for ( size_t j(0); j < npix; j++ ) {
      for (size_t k(0); k < ne; k++) {
	size_t idx = k*npix +j;
	double w = m_weightMap->model()[idx];
	bool is_null = w <= 0 || m_dataMap.data()[idx] <= 0;
	w = is_null ? 0 : w;
	wts[idx] = m_dataMap.data()[idx] * w;
      }
    }
    m_weightedCounts->setImage(wts);
  }

  void BinnedCountsCache::identifyFilledPixels() {
    fillWeightedCounts();
    const std::vector<float> & the_data = data(has_weights());
    m_filledPixels.clear();
    for (unsigned int i = 0; i < the_data.size(); i++) {
      if (the_data.at(i) > 0) {
	m_filledPixels.push_back(i);
      }
    }
  }

  void BinnedCountsCache::log_energy_ratios(const std::vector<double>& energies,
					   std::vector<double>& log_ratios) {
    log_ratios.resize(energies.size() -1);
    for ( size_t i(0); i < log_ratios.size(); i++ ) {
      log_ratios[i] = std::log(energies[i+1] / energies[i]);
    }
  }
  
  void BinnedCountsCache::computeCountsSpectrum() {
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

  void BinnedCountsCache::computeCountsSpectrum_wcs() {
    m_countsSpectrum.clear();
    m_countsSpectrum_wt.clear();
    size_t nx(m_dataMap.imageDimension(0));
    size_t ny(m_dataMap.imageDimension(1));
    size_t nz(m_dataMap.imageDimension(2));
    size_t indx(0);
    for (size_t k = 0; k < nz; k++) {
      double ntot(0);
      double ntot_wt(0);
      for (size_t j = 0; j < ny; j++) {
	for (size_t i = 0; i < nx; i++) {
	  ntot += data(false)[indx];
	  ntot_wt += data(has_weights())[indx];
	  indx++;
	}
      }
      m_countsSpectrum.push_back(ntot);
      m_countsSpectrum_wt.push_back(ntot_wt);
    }
  }

  void BinnedCountsCache::computeCountsSpectrum_healpix() {
    m_countsSpectrum.clear();
    m_countsSpectrum_wt.clear();
    size_t nx(m_dataMap.imageDimension(0));
    size_t ny(m_dataMap.imageDimension(1));
    size_t indx(0);
    for (size_t k = 0; k < ny; k++) {
      double ntot(0);
      double ntot_wt(0);
      for (size_t i = 0; i < nx; i++) {
	ntot += data(false)[indx];
	ntot_wt += data(has_weights())[indx];
	indx++;
      }
      m_countsSpectrum.push_back(ntot);
      m_countsSpectrum_wt.push_back(ntot_wt);
    }
  }

} // namespace Likelihood
