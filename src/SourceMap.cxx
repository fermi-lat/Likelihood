/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.145 2017/10/07 01:29:20 echarles Exp $
 */

#include <cmath>

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "astro/SkyProj.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "st_facilities/Util.h"

#include "Likelihood/Accumulator.h"
#include "Likelihood/BinnedExposure.h"
#include "Likelihood/BinnedCountsCache.h"
#include "Likelihood/CompositeSource.h"
#include "Likelihood/CountsMapBase.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/ConvolveHealpix.h"
#include "Likelihood/AppHelpers.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Drm.h"
#include "Likelihood/FitUtils.h"
#include "Likelihood/FileUtils.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/PSFUtils.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/Source.h"
#define ST_DLL_EXPORTS
#include "Likelihood/SourceMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/SpatialFunction.h"

#include "Likelihood/WcsMap2.h"
#include "Likelihood/WeightMap.h"
#include "Likelihood/HealpixProjMap.h"


namespace Likelihood { 


void SourceMap::fill_sparse_model(const std::vector<float>& vect,
				  SparseVector<float>& sparse) {
  sparse.resize(vect.size());
  sparse.fill_from_vect(vect);
}

void SourceMap::fill_full_model(const SparseVector<float>& sparse,
				std::vector<float>& vect) {
  sparse.fill_vect(vect);
 }


SourceMap::SourceMap(const Source& src, 
		     const BinnedCountsCache * dataCache,
                     const Observation & observation, 
                     const BinnedLikeConfig config,
		     const Drm& drm,
		     const WeightMap* weights,
		     bool save_model)
   : m_src(&src), 
     m_name(src.getName()),
     m_srcType(src.getType()),
     m_dataCache(dataCache),
     m_observation(observation),
     m_meanPsf(0),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_config(config),
     m_drm(&drm),
     m_edisp_val(src.use_edisp() ? config.edisp_val() : 0),
     m_edisp_bins(src.use_edisp() ? config.edisp_bins() : 0),
     m_edisp_offset(m_edisp_bins - drm.edisp_bins()),
     m_weights(weights),
     m_mapType(FileUtils::Unknown),
     m_save_model(save_model),     
     m_model_is_local(true),
     m_drm_cache(0) {


   set_energies();

   int status = make_model();
   if ( status != 0 ) {
     throw std::runtime_error("SourceMap construction failed");
   }
   m_drm_cache = new Drm_Cache(*m_drm, *this);
   m_loaded = true;
}

SourceMap::SourceMap(const std::string & sourceMapsFile,
                     const Source& src, 
		     const BinnedCountsCache * dataCache,
		     const Observation & observation,
                     const BinnedLikeConfig config,
		     const Drm& drm,
		     const WeightMap* weights,
		     bool save_model)
  : m_src(&src),
    m_name(src.getName()),
    m_dataCache(dataCache),
    m_observation(observation),
    m_config(config),
    m_meanPsf(0),
    m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
    m_weights(weights),
    m_mapType(FileUtils::Unknown),
    m_save_model(save_model),
    m_model_is_local(true),
    m_drm(&drm),
    m_edisp_val(src.use_edisp() ? config.edisp_val() : 0),
    m_edisp_bins(0),    
    m_edisp_offset(0),
    m_drm_cache(0) {

    int status = readModel(sourceMapsFile);

  if ( status != 0 ) {
    // throw std::runtime_error("SourceMap construction failed to read model");
  }
  m_drm_cache = new Drm_Cache(*m_drm, *this);
}


SourceMap::SourceMap(const SourceMap& other)
  :m_src(other.m_src),
   m_name(other.m_name),
   m_filename(other.m_filename),
   m_srcType(other.m_srcType),
   m_dataCache(other.m_dataCache),
   m_observation(other.m_observation),
   m_meanPsf( other.m_meanPsf != 0 ? new MeanPsf(*other.m_meanPsf) : 0),
   m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
   m_config(other.m_config),
   m_drm(other.m_drm),
   m_edisp_val(other.m_edisp_val),
   m_edisp_bins(other.m_edisp_bins),
   m_edisp_offset(other.m_edisp_offset),
   m_weights(other.m_weights),
   m_save_model(other.m_save_model),
   m_model_is_local(other.m_model_is_local),
   m_energies(other.m_energies),
   m_logEnergyRatios(other.m_logEnergyRatios),
   m_model(other.m_model),
   m_sparseModel(other.m_sparseModel),
   m_mapType(other.m_mapType),
   m_specVals(other.m_specVals),
   m_specWts(other.m_specWts),
   m_modelPars(other.m_modelPars),
   m_latchedModelPars(other.m_latchedModelPars),
   m_derivs(other.m_derivs),
   m_npreds(other.m_npreds),
   m_weighted_npreds(other.m_weighted_npreds),
   m_drm_cache(other.m_drm_cache->clone()),
   m_dataCleared(other.m_dataCleared),
   m_loaded(other.m_loaded) {

}


SourceMap::~SourceMap() {
   delete m_formatter;
   delete m_meanPsf; 
   delete m_drm_cache;
}

float SourceMap::operator[](size_t idx) const {
  return m_mapType == FileUtils::HPX_Sparse ? find_value(idx) : m_model[idx];
}


void SourceMap::sparsify_model(bool clearFull) {
  fill_sparse_model(m_model,m_sparseModel);
  if ( clearFull ) {
    // This deallocates the memory used by the model
    // in C++-11 there is a function shrink_to_fit that we could use.
    std::vector<float> nullVect;
    m_model.swap(nullVect);
  }
}
  

void SourceMap::expand_model(bool clearSparse) {
  fill_full_model(m_sparseModel,m_model);
  if ( clearSparse ) {
    m_sparseModel.clear();
  }
}


void SourceMap::computeNpredArray() {

   if ( m_mapType == FileUtils::HPX_Sparse && m_model.size() == 0 ) {
     return computeNpredArray_sparse();     
   }

   if ( m_model.size() == 0 ) {
     // The model was clear, re-make it
     // Note that this call will also call computeNpredArray,
     // so we can return now
     const std::vector<float>& dummy = model(true);
     return;
   }

   // For the npreds we are going to loop over the number of energies in the source map
   size_t ne = n_energies();
   // For the weighted npreds we are going to loop over the number of measure energy bins
   size_t nw = m_dataCache->num_ebins();

   // For the weighted npreds, we loop from edisp_bins() to edisp_bins() + nw
   size_t min_ebin = edisp_bins();
   size_t max_ebin = min_ebin + nw;

   m_npreds.clear();
   m_npreds.resize(ne);
   m_weighted_npreds.clear();
   m_weighted_npreds.resize(nw);

   size_t npix = m_dataCache->num_pixels();
   size_t kmin_edisp(0.);
   size_t kmax_edisp(0.);            
   size_t kmeas(0);

   // k is the index in the output vector
   for (size_t k(0); k < ne; k++) {      

     bool is_meas_bin(false);

     if ( k >= min_ebin && k < max_ebin ) {
       is_meas_bin = true;
       kmeas = k - edisp_bins();
       FitUtils::get_edisp_range(*this, kmeas, kmin_edisp, kmax_edisp);
       if ( edisp_val() <= 0 ) {
	 kmin_edisp += edisp_bins();
	 kmax_edisp += edisp_bins();
       }
       m_weighted_npreds[kmeas].resize(kmax_edisp-kmin_edisp, std::make_pair(0., 0.));
     }     

     for (size_t j(0); j < npix; j++) {
       // This in the index in the model
       size_t indx_0(k*npix + j);
       double model_0 = m_model[indx_0];
       m_npreds[k] += model_0;
       // If this energy bin is not in the central range, we are done with this pixel
       if ( ! is_meas_bin ) continue;
       // This is the index in the weights
       size_t indx_w(kmeas*npix + j);
       double weight_val = m_weights != 0 ? m_weights->model()[indx_w] : 1.;       
       size_t idx_k(0);
       size_t idx_k0(kmin_edisp*npix + j);
       for ( size_t kk(kmin_edisp); kk < kmax_edisp; kk++, idx_k++, idx_k0 += npix ) {
	 size_t idx_k1 = idx_k0 + npix;
	 m_weighted_npreds.at(kmeas).at(idx_k).first += weight_val *  m_model[idx_k0];
	 m_weighted_npreds.at(kmeas).at(idx_k).second += weight_val *  m_model[idx_k1];
       }
     }   
   }
}

void SourceMap::computeNpredArray_sparse() {
  
   // For the npreds we are going to loop over the number of energies in the source map
   size_t ne = n_energies();
   // For the weighted npreds we are going to loop over the number of measure energy bins
   size_t nw = m_dataCache->num_ebins();

   // For the weighted npreds, we loop from edisp_bins() to edisp_bins() + 
   size_t min_ebin = edisp_bins();
   size_t max_ebin = min_ebin + nw;

   m_npreds.clear();
   m_npreds.resize(ne);
   m_weighted_npreds.clear();
   m_weighted_npreds.resize(nw);

   size_t npix = m_dataCache->num_pixels();
   size_t k(0);
   size_t kmin_edisp(0.);
   size_t kmax_edisp(0.);     
   size_t kmeas(0);
   bool is_meas_bin(false);

   for (size_t k=0; k < ne; k++) {
     if ( k >= min_ebin && k < max_ebin ) {
       kmeas = k - edisp_bins();
       FitUtils::get_edisp_range(*this, kmeas, kmin_edisp, kmax_edisp);
       if ( edisp_val() < 0 ) {
	 kmin_edisp += edisp_bins();
	 kmax_edisp += edisp_bins();
       } else if ( edisp_val() > 0 ) {
	 kmin_edisp += edisp_offset();
	 kmax_edisp += edisp_offset();
       }
       m_weighted_npreds[kmeas].resize(kmax_edisp-kmin_edisp, std::make_pair(0., 0.));
     }     
   }

   for ( SparseVector<float>::iterator itr = m_sparseModel.begin(); itr != m_sparseModel.end(); itr++ ) {
     size_t indx_0 = itr->first;
     size_t k_true = indx_0 / npix;
     size_t j = indx_0 % npix;
     double model_0 = itr->second;
     m_npreds[k_true] += model_0;
     if ( k_true < min_ebin || k_true >= max_ebin ) {
       continue;
     }
     size_t kmeas = k_true - edisp_bins();
     size_t indx_w(kmeas*npix + j);

     double weight_val = m_weights != 0 ? m_weights->model()[indx_w] : 1.;
     FitUtils::get_edisp_range(*this, kmeas, kmin_edisp, kmax_edisp);
     size_t idx_k(0);
     size_t idx_k0(kmin_edisp*npix + j);
     for ( size_t kk(kmin_edisp); kk < kmax_edisp; kk++, idx_k++, idx_k0 += npix ) {
       size_t idx_k1 = idx_k0 + npix;
       m_weighted_npreds.at(kmeas).at(idx_k).first += weight_val * m_sparseModel[idx_k0];
       m_weighted_npreds.at(kmeas).at(idx_k).second += weight_val * m_sparseModel[idx_k1];
     }
   }
}

  void SourceMap::applyPhasedExposureMap() {
    if (!m_observation.have_phased_expmap()) {
      return;
    }
    if ( m_mapType == FileUtils::HPX_Sparse ) {
      return applyPhasedExposureMap_sparse();
    }

    const ProjMap * phased_expmap = &(m_observation.phased_expmap());
    const std::vector<Pixel> & pixels = m_dataCache->countsMap().pixels();
    for (size_t k(0); k < n_energies(); k++) {
      std::vector<Pixel>::const_iterator pixel(pixels.begin());
      for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
	size_t indx(k*pixels.size() + j);
	m_model[indx] *= phased_expmap->operator()(pixel->dir(),
						   m_energies[k]);
      }
    }
  }

  void SourceMap::applyPhasedExposureMap_sparse() {
    if (!m_observation.have_phased_expmap()) {
      return;
    }
    const ProjMap * phased_expmap = &(m_observation.phased_expmap());
    const std::vector<Pixel> & pixels = m_dataCache->countsMap().pixels();
   
    for ( SparseVector<float>::iterator itr = m_sparseModel.begin(); itr != m_sparseModel.end(); itr++ ) {
      size_t indx = itr->first;
      size_t j = indx % pixels.size();
      size_t k = indx / pixels.size();
      const Pixel& pixel = pixels[j];
      itr->second *= phased_expmap->operator()(pixel.dir(), m_energies[k]);
    }
  }

  void SourceMap::setSource(const Source& src) {
    if ( m_src == &src ) {
      if(printDebug) std::cout << "SourceMap::setSource for " << m_name << " - m_loaded=" << std::boolalpha << m_loaded << ", m_save_model = " << m_save_model
                << ", m_model.size()=" << m_model.size() <<  ", m_sparseModel.size()=" << m_sparseModel.size() << std::endl;
      if (!m_loaded || (m_save_model && m_model.size() == 0 && m_sparseModel.size() == 0) ) {
        if(printDebug) std::cout << "SourceMap::setSource for " << m_name << " - passed - m_loaded=" << std::boolalpha << m_loaded << ", m_save_model = " << m_save_model
                  << ", m_model.size()=" << m_model.size() <<  ", m_sparseModel.size()=" << m_sparseModel.size() << std::endl;
        getSourceData();
      } else {
	      return;
      }
    }
    m_src = &src;
    resetSourceData();
  }


  const Drm_Cache* SourceMap::update_drm_cache(const Drm* drm, bool force) {
    bool changed = m_drm != drm;
    m_drm = drm;
    return drm_cache(force || changed);
  }


  void SourceMap::setSpectralValues(bool latch_params ) {
    if ( m_src == 0 ) return;
    FitUtils::extractSpectralVals(*m_src,m_energies,m_specVals);
    FitUtils::get_spectral_weights(m_specVals,m_energies,m_logEnergyRatios,m_specWts);
    m_modelPars.clear();
    m_src->spectrum().getParamValues(m_modelPars);
    if ( latch_params ) {
      m_latchedModelPars.resize(m_modelPars.size());
      std::copy(m_modelPars.begin(),m_modelPars.end(),m_latchedModelPars.begin());
    }
  }

  void SourceMap::setSpectralDerivs(const std::vector<std::string>& paramNames) {
    if ( m_src == 0 ) return;
    FitUtils::extractSpectralDerivs(*m_src,m_energies,paramNames,m_derivs);
  }


  const std::vector<double>& SourceMap::counts_spectra(int edisp_val,
						       bool use_weighted) const {
     const std::vector<double>& the_vect = use_weighted ? 
       ( edisp_val == 0 ? m_drm_cache->true_counts_wt() : m_drm_cache->meas_counts_wt() ) :
       ( edisp_val == 0 ? m_drm_cache->true_counts() : m_drm_cache->meas_counts() );
     return the_vect;
   }


bool SourceMap::spectrum_changed() const {
  if ( m_src == 0 ) return true;
  std::vector<double> parValues;
  m_src->spectrum().getParamValues(parValues);
  if ( parValues.size() != m_latchedModelPars.size() ) {
    return true;
  }
  for (size_t j(0); j < parValues.size(); j++) {
    if (parValues[j] != m_latchedModelPars[j] ) {
      return true;
    }
  }
  return false;
}


std::vector<float> & SourceMap::model(bool force) {
  if ( m_model.size() == 0 || force ) {
    if ( m_filename.size() > 0 && FileUtils::fileHasExtension(m_filename, m_name) ) {
      readModel(m_filename);
    } else {
      int status = make_model();
      if ( status != 0 ) {
	throw std::runtime_error("SourceMap model");
      }
    }
  }
  return m_model;
}
  

const std::vector<double> & SourceMap::specVals(bool force, bool latch_params) {
  if ( m_specVals.size() == 0 || force ) {
    // Decide if we need the full or truncated version of the energy vector
    setSpectralValues(latch_params);
  }
  return m_specVals;
}

const std::vector<std::pair<double,double> >&  SourceMap::specWts(bool force) {
  if ( m_specWts.size() == 0 || force ) {
    setSpectralValues();
  }
  return m_specWts;
}

const std::vector<std::vector<double> >& SourceMap::specDerivs(const std::vector<std::string>& paramNames, bool force) {
  if ( m_specVals.size() == 0 || force ) {
    setSpectralDerivs(paramNames);
  }
  return m_derivs;
}


const std::vector<double> & SourceMap::npreds(bool force) {
  if ( m_npreds.size() == 0 || force ) {
    computeNpredArray();
  }
  return m_npreds;
}
  
const std::vector<std::vector<std::pair<double,double> > >& SourceMap::weighted_npreds(bool force) {
  if ( m_weighted_npreds.size() == 0 || force ) {
    computeNpredArray();
  }
  return m_weighted_npreds;
}  

const Drm_Cache* SourceMap::drm_cache(bool force) {
  if ( m_drm_cache == 0) {
    m_drm_cache = new Drm_Cache(*m_drm, *this);
  } else if ( force ) {
    m_drm_cache->update(*m_drm, *this);  
  }
  return m_drm_cache;
}
 
void SourceMap::addToVector(std::vector<float>& vect, bool includeSpec, int kmin, int kmax) {
  switch ( m_mapType ) {
  case FileUtils::HPX_Sparse:
    return addToVector_sparse(vect,includeSpec,kmin,kmax);
    break;
  case FileUtils::WCS:
  case FileUtils::HPX_AllSky:
  case FileUtils::HPX_Partial:
  default:
    return addToVector_full(vect,includeSpec,kmin,kmax);
  }
}
   
void SourceMap::subtractFromVector(std::vector<float>& vect, bool includeSpec, int kmin, int kmax){
  switch ( m_mapType ) {
  case FileUtils::HPX_Sparse:
    return subtractFromVector_sparse(vect,includeSpec,kmin,kmax);
    break;
  case FileUtils::WCS:
  case FileUtils::HPX_AllSky:
  case FileUtils::HPX_Partial:
  default:
    return subtractFromVector_full(vect,includeSpec,kmin,kmax);
  }
}

double SourceMap::summed_counts(size_t kmin, size_t kmax,
				int edisp_val,
				bool use_weighted) {
    
  if ( m_drm_cache == 0 ) {
    throw std::runtime_error("Called summed counts on a SourceMap that does not have a Drm_Cache ");
  }
 
  const std::vector<double>& the_vect = counts_spectra(edisp_val, use_weighted);
 
  double retVal(0.);
  for ( size_t k(kmin); k < kmax; k++ ) {
    retVal += the_vect[k];
  }
  return retVal;
}

void SourceMap::setImage(const std::vector<float>& model) {
  if(model.size() != m_model.size() && m_model.size() != 0) {
    throw std::runtime_error("Wrong size for input model map.");
  }
  m_model = model;
  m_model_is_local = true;
  applyPhasedExposureMap();
  computeNpredArray();
}

void SourceMap::setWeights(const WeightMap* weights) {
  m_weights = weights;
  computeNpredArray();
}


size_t SourceMap::memory_size() const {
  size_t retVal(0);
  retVal += sizeof(*this);
  retVal += m_name.capacity();
  retVal += m_filename.capacity();
  retVal += m_srcType.capacity();
  retVal += sizeof(*m_formatter);
  retVal += sizeof(float)*m_model.capacity();
  retVal += sizeof(std::pair<size_t,float>)*m_sparseModel.capacity();
  retVal += sizeof(double)*m_modelPars.capacity();
  retVal += sizeof(double)*m_npreds.capacity();
  retVal += sizeof(std::pair<double,double>)*m_weighted_npreds.capacity();
  if ( m_drm_cache != 0 ) {
    retVal += m_drm_cache->memory_size();
  }
  return retVal;
}


void SourceMap::test_sparse(const std::string& prefix) const {
  for ( SparseVector<float>::const_iterator itr = m_sparseModel.begin();
	itr != m_sparseModel.end(); itr++ ) {
    if ( itr->first >= m_sparseModel.size() ) {
    }
  }
}


int SourceMap::readModel(const std::string& filename) {
  m_model.clear();
  m_filename = filename;
  m_model_is_local = false;
  //m_dataCleared = false;

  m_specVals.clear();
  m_specWts.clear();
  m_modelPars.clear();
  m_derivs.clear();
  m_npreds.clear();
  m_weighted_npreds.clear();
  int status(0);
  switch ( m_dataCache->countsMap().projection().method()  ) {
  case astro::ProjBase::WCS:
    status = readImage(m_filename);
    break;
  case astro::ProjBase::HEALPIX:
    status = readTable_healpix(m_filename);
    break;
  default:
    status = -1;
    break;
  }

  if ( status != 0 ) {
    std::string errMsg("SourceMap failed to read source map: ");
    errMsg += m_filename;
    errMsg += "::";
    errMsg += m_name;
    errMsg += ".  To match data file: ";
    errMsg += m_dataCache->countsMap().filename();
    throw std::runtime_error(errMsg);
  }

  applyPhasedExposureMap();
  setSpectralValues();
  computeNpredArray();
  m_loaded = true;
    
  return status;
}

int SourceMap::readImage(const std::string& sourceMapsFile) {
  m_mapType = FileUtils::get_src_map_type(sourceMapsFile,m_name);
  int nEvals(0.);
  FileUtils::read_fits_image_to_float_vector(sourceMapsFile,m_name,m_model,nEvals);
  m_edisp_bins = size_t((nEvals - m_dataCache->num_energies())/2);
  m_edisp_offset = m_edisp_bins - m_drm->edisp_bins();

  size_t n_edisp_needed = m_edisp_val > 0 ? m_edisp_val : 0;
  if ( m_edisp_bins < n_edisp_needed ) {
    std::ostringstream message;
    message << "SourceMap::readImage "
	    << "not enough energy bins to apply energy dispersion for source " 
	    << m_src->getName() << ' ' << m_edisp_bins << " < " << n_edisp_needed;    
    throw std::runtime_error(message.str());
  }

  m_energies.resize(m_dataCache->num_energies());
  std::copy(m_dataCache->energies().begin(),m_dataCache->energies().end(),m_energies.begin());
  FitUtils::expand_energies(m_energies, m_edisp_bins);
  FitUtils::log_energy_ratios(m_energies, m_logEnergyRatios);
  return 0;
}

int SourceMap::readTable_healpix(const std::string& sourceMapsFile) {
  m_mapType = FileUtils::get_src_map_type(sourceMapsFile,m_name);
  int status(0);
  int nEnergy(0);
  switch (m_mapType) {
  case FileUtils::HPX_AllSky:
  case FileUtils::HPX_Partial:
    // In either of these two cases we simple read the vector.  
    // If this is a partial-sky mapping, the projection will 
    // take care of doing the remapping
    status = FileUtils::read_healpix_table_to_float_vector(sourceMapsFile,m_name,m_model,nEnergy);
    status = m_model.size() > 0 ? status : -1;
    break;
  case FileUtils::HPX_Sparse:
    // In this case we read the map.
    /// FIXME, we should have a better way of getting this...
    m_sparseModel.resize( m_dataCache->source_map_size() );
    status = FileUtils::read_healpix_table_to_sparse_vector(sourceMapsFile,m_name,
							    m_dataCache->num_pixels(),
							    m_sparseModel,nEnergy);
    status = m_sparseModel.size() > 0 ? status : -1;  
    break;
  default:
    // Either unknown or WCS based.  This is an error in either case.
    return -2;
  }
  m_edisp_bins = int((nEnergy -  m_dataCache->num_energies())/2);
  m_energies.resize(m_dataCache->num_energies());
  std::copy(m_dataCache->energies().begin(),m_dataCache->energies().end(),m_energies.begin());
  FitUtils::expand_energies(m_energies, m_edisp_bins);
  FitUtils::log_energy_ratios(m_energies, m_logEnergyRatios);
  return status;
}


int SourceMap::make_model() {
  if ( m_src == 0 ) return -1;
  
  m_model_is_local = true;
  m_model.clear();
  m_specVals.clear();
  m_specWts.clear();
  m_modelPars.clear();
  m_derivs.clear();
  m_npreds.clear();
  m_weighted_npreds.clear();

  int status(0);

  if ( m_src->srcType() == Source::Point && 
       !m_config.psf_integ_config().use_single_psf() &&
       m_meanPsf == 0 ) {
    m_meanPsf = PSFUtils::build_psf(*m_src, m_dataCache->countsMap(), m_energies, m_observation);
  }

  status = PSFUtils::makeModelMap(*m_src, *m_dataCache, 
				  m_energies,
				  m_meanPsf==0 ? m_observation.meanpsf() : *m_meanPsf,
				  m_observation.bexpmap_ptr(),
				  m_config.psf_integ_config(),
				  m_filename, m_drm,
				  *m_formatter, m_model, m_mapType);

  if ( status != 0 ) { 
    return status;
  }

  // Sparsify the model, if needed
  if ( m_mapType == FileUtils::HPX_Sparse ) {
    sparsify_model();
  }

  applyPhasedExposureMap();
  computeNpredArray();
  setSpectralValues();

  return status;
}

void SourceMap::addToVector_full(std::vector<float>& vect, bool includeSpec, int kmin, int kmax) const {
  SourceMap* nc_this = const_cast<SourceMap*>(this);
  const std::vector<float>& m = const_cast<SourceMap*>(this)->model();
  size_t npix = m_dataCache->num_pixels();
  kmax = kmax < 0 ? n_energies() : kmax;
  size_t ne = kmax - kmin;
  size_t nval = ne*npix;
  if ( vect.size() != nval ) {
    throw std::runtime_error("SourceMap::addToVector_full model size != vector size");
  }
  if ( includeSpec && m_specVals.size() != n_energies() ) {
    throw std::runtime_error("SourceMap::addToVector_full spectrum size != number of energy layers");
  }
  std::vector<float>::const_iterator itr_in = m.begin() + (npix*kmin);
  std::vector<float>::iterator itr_out = vect.begin();
  for ( size_t ie(kmin); ie != kmax; ie++ ) {
    double factor = includeSpec ? m_specVals[ie] : 1.;
    for ( size_t ipix(0); ipix < npix; ipix++, itr_in++, itr_out++ ) {
      *itr_out += *itr_in * factor;
    }
  }
}
  
void SourceMap::addToVector_sparse(std::vector<float>& vect, bool includeSpec, int kmin, int kmax) const {
  size_t npix = m_dataCache->num_pixels();
  kmax = kmax < 0 ? n_energies() : kmax;
  size_t ne = kmax - kmin;
  size_t nval = ne*npix;

  if ( vect.size() != nval ) {
    throw std::runtime_error("SourceMap::addToVector_sparse model size != vector size");
  }
  if ( includeSpec ) {
    SourceMap* nct = const_cast<SourceMap*>(this);
    nct->setSpectralValues();
  }
  SparseVector<float>::const_iterator itr = m_sparseModel.lower_bound(kmin*npix);
  SparseVector<float>::const_iterator itr_end = m_sparseModel.upper_bound(kmax*npix);
  size_t offset = kmin*npix;
  for ( ; itr != itr_end; itr++) {
    int ie = includeSpec ? ( itr->first / npix ) : -1;    
    if ( ie >= m_specVals.size() ) {
      std::string errMsg("Sparse map index is outside of bounds.");
      errMsg += "This usually indicates that one of the pixels in the map was corrupted with writing the map";    
      throw std::runtime_error(errMsg);
    }
    double factor = ie >= 0 ? m_specVals[ie] : 1.;
    vect[itr->first-offset] += itr->second * factor;
  }
}

void SourceMap::subtractFromVector_full(std::vector<float>& vect, bool includeSpec, int kmin, int kmax) const {
  const std::vector<float>& m = const_cast<SourceMap*>(this)->model();
  size_t npix = m_dataCache->num_pixels();
  kmax = kmax < 0 ? n_energies() : kmax;
  size_t ne = kmax - kmin;
  size_t nval = ne*npix;
  if ( vect.size() != nval ) {
    throw std::runtime_error("SourceMap::subtractFromVector_full model size != vector size");
  }
  if ( includeSpec ) {
    SourceMap* nct = const_cast<SourceMap*>(this);
    nct->setSpectralValues();
  }
  std::vector<float>::const_iterator itr_in = m_model.begin() + (npix*kmin);
  std::vector<float>::iterator itr_out = vect.begin();
  for ( size_t ie(kmin); ie != kmax; ie++ ) {
    double factor = includeSpec ? m_specVals[ie] : 1.;
    for ( size_t ipix(0); ipix < npix; ipix++, itr_in++, itr_out++ ) {
      *itr_out -= *itr_in * factor;
    }
  }
}
 

void SourceMap::subtractFromVector_sparse(std::vector<float>& vect, bool includeSpec, int kmin, int kmax) const {
  size_t npix = m_dataCache->num_pixels();
  kmax = kmax < 0 ? n_energies() : kmax;
  size_t ne = kmax - kmin;
  size_t nval = ne*npix;

  if ( vect.size() != nval ) {
    throw std::runtime_error("SourceMap::subtractFromVector_sparse model size != vector size");
  }
  if ( includeSpec && m_specVals.size() != n_energies() ) {
    throw std::runtime_error("SourceMap::subtractFromVector_sparse spectrum size != number of energy layers");
  }
  SparseVector<float>::const_iterator itr  = m_sparseModel.lower_bound(kmin*npix);
  SparseVector<float>::const_iterator itr_end = m_sparseModel.upper_bound(kmax*npix);
  size_t offset = kmin*npix;
  for ( ; itr != itr_end; itr++) {
    int ie = includeSpec ? ( itr->first / npix ) : -1;
    double factor = ie >= 0 ? m_specVals[ie] : 1.;
    vect[itr->first-offset] -= itr->second * factor;
  }
} 

  void SourceMap::set_energies(bool reload) {
    if (!reload){
      m_edisp_val = m_src->use_edisp() ? m_config.edisp_val() : 0;
    }
    m_edisp_offset = m_edisp_bins - m_drm->edisp_bins(),
    m_energies.resize(m_dataCache->num_energies());
    std::copy(m_dataCache->energies().begin(),m_dataCache->energies().end(),m_energies.begin());
    FitUtils::expand_energies(m_energies, m_edisp_bins);
    FitUtils::log_energy_ratios(m_energies, m_logEnergyRatios);
  }

  void SourceMap::reloadIfCleared(){
    if (m_dataCleared){
      if(printDebug) std::cout << "SourceMap::reloadIfCleared() - reloading for " << m_name << std::endl;
      getSourceData(false);
      m_dataCleared = false;
    }
    resetSourceData(false);
  }

  void SourceMap::resetSourceData(bool reload){
    set_energies(reload);
    m_specVals.clear();
    m_specWts.clear();
    m_modelPars.clear();
    m_derivs.clear();
    m_npreds.clear();
    m_weighted_npreds.clear();  
    m_model_is_local = false;
  }

  void SourceMap::getSourceData(bool reload){
    m_loaded = true;
    if (m_filename.size() > 0) {
      if(printDebug) std::cout << "sourceMap::getSourceData() - reading data for " << m_name << std::endl;
      readModel(m_filename);
    } else {
      if(printDebug) std::cout << "sourceMap::getSourceData() - building data for " << m_name << std::endl;
      set_energies(reload);
      make_model();
      return;
    }
  }

} // Likelihood
