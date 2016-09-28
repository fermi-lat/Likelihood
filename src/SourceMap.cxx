/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.128 2016/09/22 01:38:09 echarles Exp $
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

#include "Likelihood/BinnedExposure.h"
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
		     const CountsMapBase * dataMap,
                     const Observation & observation, 
                     const PsfIntegConfig & psf_config,
		     const Drm* drm,
		     const WeightMap* weights,
		     bool save_model)
   : m_src(&src), 
     m_name(src.getName()),
     m_srcType(src.getType()),
     m_dataMap(dataMap),
     m_observation(observation),
     m_meanPsf(0),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_psf_config(psf_config),
     m_drm(drm),
     m_weights(weights),
     m_mapType(FileUtils::Unknown),
     m_save_model(save_model),
     m_drm_cache(0) {
   
   int status = make_model();
   if ( status != 0 ) {
     throw std::runtime_error("SourceMap construction failed");
   }
   m_drm_cache = new Drm_Cache(m_drm,*this,dataMap->energies());  
}

SourceMap::SourceMap(const std::string & sourceMapsFile,
                     const Source& src, 
		     const CountsMapBase * dataMap,
		     const Observation & observation,
		     const WeightMap* weights,
		     const Drm* drm,
		     bool save_model) 
  : m_src(&src),
    m_name(src.getName()),
    m_dataMap(dataMap),
    m_observation(observation),
    m_meanPsf(0),
    m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
    m_weights(weights),
    m_mapType(FileUtils::Unknown),
    m_save_model(save_model),
    m_drm(drm),
    m_psf_config(),
    m_drm_cache(0) {


  int status = readModel(sourceMapsFile);
  if ( status != 0 ) {
    throw std::runtime_error("SourceMap construction failed to read model");
  }
  m_drm_cache = new Drm_Cache(m_drm,*this,m_dataMap->energies());
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
   const std::vector<Pixel> & pixels(m_dataMap->pixels());
   
   std::vector<double> energies;
   switch ( m_dataMap->projection().method() ) {
   case astro::ProjBase::WCS:
     m_dataMap->getAxisVector(2, energies);
     break;
   case astro::ProjBase::HEALPIX:
     m_dataMap->getAxisVector(1, energies);
     break;
   }
   

   // If we are making the __weights__ source map,
   // the output npreds is one less than the input map size
   // since we have the energies at the bin edges, but we need the
   // weights integrated over the bins
   size_t ne = energies.size();
   // The number of weights is always number of energy bins
   size_t nw = energies.size() -1 ;
   m_npreds.clear();
   m_npreds.resize(ne);
   m_npred_weights.clear();
   m_npred_weights.resize(nw, std::make_pair<double,double>(0.,0.));

   bool expanded = false;
   if ( m_mapType == FileUtils::HPX_Sparse && m_model.size() == 0 ) {
     expanded = true;
     expand_model(false);
   }

   for (size_t k(0); k < ne; k++) {

      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      double w_0_sum(0.);
      double w_1_sum(0.);

      for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
	 size_t indx(k*pixels.size() + j);
         size_t indx_0 = k > 0 ? indx  - pixels.size() : indx;
	 size_t indx_1 = k < energies.size()-1 ? indx : indx - pixels.size();
	 double addend = m_model.at(indx);
         m_npreds[k] += addend;
	 double w_0_addend = m_weights != 0 ? ( m_weights->model()[indx_0]*addend ) : addend;
	 double w_1_addend = m_weights != 0 ? ( m_weights->model()[indx_1]*addend ) : addend;
	 w_0_sum += w_0_addend;
	 w_1_sum += w_1_addend;
      }
      double w_0 = m_weights != 0 ? (m_npreds[k] > 0 ? w_0_sum / m_npreds[k] : 0.) : 1.0;
      double w_1 = m_weights != 0 ? (m_npreds[k] > 0 ? w_1_sum / m_npreds[k] : 0.) : 1.0;

      if ( k < nw ) {
	m_npred_weights[k].first = w_0;
      }
      if ( k > 0 ) {
	m_npred_weights[k-1].second = w_1;
      }
   }

   if ( expanded ) {
     m_model.clear();
   }

}


void SourceMap::applyPhasedExposureMap() {
   if (!m_observation.have_phased_expmap()) {
      return;
   }
   const ProjMap * phased_expmap = &(m_observation.phased_expmap());
   const std::vector<Pixel> & pixels(m_dataMap->pixels());
   std::vector<double> energies;
   m_dataMap->getEnergies(energies);
   for (size_t k(0); k < energies.size(); k++) {
      std::vector<Pixel>::const_iterator pixel(pixels.begin());
      for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
         size_t indx(k*pixels.size() + j);
         m_model.at(indx) *= phased_expmap->operator()(pixel->dir(),
                                                       energies[k]);
      }
   }
}


void SourceMap::setSource(const Source& src) {
  if ( m_src == &src ) {
    return;
  }
  m_src = &src;
  m_specVals.clear();
  m_modelPars.clear();
  m_derivs.clear();
  m_npreds.clear();
  m_npred_weights.clear();  
}

void SourceMap::setSpectralValues(const std::vector<double>& energies,
				    bool latch_params ) {
  if ( m_src == 0 ) return;
  FitUtils::extractSpectralVals(*m_src,energies,m_specVals);
  m_modelPars.clear();
  m_src->spectrum().getParamValues(m_modelPars);
  if ( latch_params ) {
    m_latchedModelPars.resize(m_modelPars.size());
    std::copy(m_modelPars.begin(),m_modelPars.end(),m_latchedModelPars.begin());
  }
}


void SourceMap::setSpectralDerivs(const std::vector<double>& energies,
				  const std::vector<std::string>& paramNames) {
  if ( m_src == 0 ) return;
  FitUtils::extractSpectralDerivs(*m_src,energies,paramNames,m_derivs);
}


bool SourceMap::spectrum_changed() const {
  if ( m_src == 0 ) return true;
  std::vector<double> parValues;
  m_src->spectrum().getParamValues(parValues);
  if ( parValues.size() != m_latchedModelPars.size() ) {
    return true;
  }
  for (size_t j(0); j < parValues.size(); j++) {
    if (parValues.at(j) != m_latchedModelPars.at(j) ) {
      return true;
    }
  }
  return false;
}


const std::vector<float> & SourceMap::model(bool force) {
  if ( m_model.size() == 0 || force ) {
    if ( m_filename.size() > 0 ) {        
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
  

const std::vector<double> & SourceMap::specVals(bool force) {
  if ( m_specVals.size() == 0 || force ) {
    setSpectralValues(m_dataMap->energies());
  }
  return m_specVals;
}


const std::vector<std::vector<double> >& SourceMap::specDerivs(const std::vector<std::string>& paramNames, bool force) {
  if ( m_specVals.size() == 0 || force ) {
    setSpectralDerivs(m_dataMap->energies(),paramNames);
  }
  return m_derivs;
}


const std::vector<double> & SourceMap::npreds(bool force) {
  if ( m_npreds.size() == 0 || force ) {
    computeNpredArray();
  }
  return m_npreds;
}
  
const std::vector<std::pair<double,double> > & SourceMap::npred_weights(bool force) {
  if ( m_npred_weights.size() == 0 || force ) {
    computeNpredArray();
  }
  return m_npred_weights;
}  

const Drm_Cache* SourceMap::drm_cache(bool force) {
  if ( m_drm_cache == 0) {
    m_drm_cache = new Drm_Cache(m_drm,*this,m_dataMap->energies());  
  } else if ( force ) {
    m_drm_cache->update(m_drm,*this,m_dataMap->energies());  
  }
  return m_drm_cache;
}

double SourceMap::summed_counts(size_t kmin, size_t kmax,
				bool use_edisp,
				bool use_weighted) {

  if ( m_drm_cache == 0 ) {
    throw std::runtime_error("SourceMap::summed_counts: no Drm_Cache");
  }  

  const std::vector<double>& counts_spec = use_edisp ? 
    ( use_weighted ? m_drm_cache->meas_counts_wt() : m_drm_cache->meas_counts() ) :
    ( use_weighted ? m_drm_cache->true_counts_wt() : m_drm_cache->true_counts() );
     
  double value(0.);
  for (size_t k(kmin); k < kmax; k++) {    
    value += counts_spec[k];
  }
  return value;    
}

void SourceMap::setImage(const std::vector<float>& model) {
  if(model.size() != m_model.size())
    throw std::runtime_error("Wrong size for input model map.");

  m_model = model;
  m_filename.clear();
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
  retVal += sizeof(std::pair<double,double>)*m_npred_weights.capacity();
  if ( m_drm_cache != 0 ) {
    retVal += m_drm_cache->memory_size();
  }
  return retVal;
}


int SourceMap::readModel(const std::string& filename) {
  m_model.clear();
  m_filename = filename;
  
  m_specVals.clear();
  m_modelPars.clear();
  m_derivs.clear();
  m_npreds.clear();
  m_npred_weights.clear();

  int status(0);
  switch ( m_dataMap->projection().method()  ) {
  case astro::ProjBase::WCS:
    status = readImage(m_filename);
    break;
  case astro::ProjBase::HEALPIX:
    status = readTable_healpix(m_filename);
    break;
  default:
    break;
  }
  if ( status != 0 ) {
    std::string errMsg("SourceMap failed to read source map: ");
    errMsg += m_filename;
    errMsg += ".  To match data file: ";
    errMsg += m_dataMap->filename();
    throw std::runtime_error(errMsg);
  }

  // FIXME, we could be more efficient about this
  if ( m_mapType == FileUtils::HPX_Sparse ) {
    expand_model(false);
  }

  applyPhasedExposureMap();
  computeNpredArray();
  setSpectralValues(m_dataMap->energies());

  // FIXME, we could be more efficient about this
  if ( m_mapType == FileUtils::HPX_Sparse ) {
    m_model.clear();
  }

  return status;
}

int SourceMap::readImage(const std::string& sourceMapsFile) {
  m_mapType = FileUtils::get_src_map_type(sourceMapsFile,m_name);
  FileUtils::read_fits_image_to_float_vector(sourceMapsFile,m_name,m_model);
  return 0;
}

int SourceMap::readTable_healpix(const std::string& sourceMapsFile) {
  m_mapType = FileUtils::get_src_map_type(sourceMapsFile,m_name);
  int status(0);
  switch (m_mapType) {
  case FileUtils::HPX_AllSky:
  case FileUtils::HPX_Partial:
    // In either of these two cases we simple read the vector.  
    // If this is a partial-sky mapping, the projection will 
    // take care of doing the remapping
    status = FileUtils::read_healpix_table_to_float_vector(sourceMapsFile,m_name,m_model);
    break;
  case FileUtils::HPX_Sparse:
    // In this case we read the map.
    /// FIXME, we should have a better way of getting this...
    m_sparseModel.resize( ( m_dataMap->data().size() / (m_dataMap->energies().size() -1) ) * m_dataMap->energies().size() );
    status = FileUtils::read_healpix_table_to_sparse_vector(sourceMapsFile,m_name,m_sparseModel);
    break;
  default:
    // Either unknown or WCS based.  This is an error in either case.
    return -1;
  }
  
  if ( status != 0 ) return status;


  return 0;

}


int SourceMap::make_model() {
  if ( m_src == 0 ) return -1;
  
  m_filename.clear();
  m_model.clear();
  m_specVals.clear();
  m_modelPars.clear();
  m_derivs.clear();
  m_npreds.clear();
  m_npred_weights.clear();

  int status(0);
  if (m_psf_config.verbose()) {
    m_formatter->warn() << "Generating SourceMap for " << m_name;
  }
  switch ( m_src->srcType() ) {
  case Source::Diffuse:
    status = PSFUtils::makeDiffuseMap(static_cast<const DiffuseSource&>(*m_src), *m_dataMap, 
				      m_observation.meanpsf(), m_observation.bexpmap(),
				      m_psf_config, *m_formatter, m_model, m_mapType);
    break;
  case Source::Point:
    m_meanPsf = m_psf_config.use_single_psf() ? 0 : PSFUtils::build_psf(*m_src,*m_dataMap,m_observation);
    status =  PSFUtils::makePointSourceMap(static_cast<const PointSource&>(*m_src), *m_dataMap, 
					   m_psf_config, m_meanPsf==0 ? m_observation.meanpsf() : *m_meanPsf, 
					   *m_formatter, m_model, m_mapType);
    break;
  default:
    throw std::runtime_error("Unrecognized source type");
  }
  
  if ( status != 0 ) { 
    return status;
  }

  if (m_psf_config.verbose()) {
    m_formatter->warn() << "!" << std::endl;
  }
  
  applyPhasedExposureMap();
  computeNpredArray();
  setSpectralValues(m_dataMap->energies());

  // FIXME, we could be more efficient about this
  if ( m_mapType == FileUtils::HPX_Sparse ) {
    sparsify_model();
  }

  return status;
}

} // Likelihood
