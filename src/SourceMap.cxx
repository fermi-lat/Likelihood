/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.123 2016/09/13 19:26:23 echarles Exp $
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

SourceMap::SourceMap(const Source& src, 
		     const CountsMapBase * dataMap,
                     const Observation & observation, 
                     const PsfIntegConfig & psf_config,
		     const Drm* drm,
		     const WeightMap* weights)
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
		     const Drm* drm) 
  : m_src(&src),
    m_name(src.getName()),
    m_dataMap(dataMap),
    m_observation(observation),
    m_meanPsf(0),
    m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
    m_weights(weights),
    m_drm(drm),
    m_psf_config(),
    m_drm_cache(0) {
  readModel(sourceMapsFile);
  applyPhasedExposureMap();
  computeNpredArray();
  setSpectralValues(m_dataMap->energies());
  m_drm_cache = new Drm_Cache(m_drm,*this,m_dataMap->energies());
}


SourceMap::~SourceMap() {
   delete m_formatter;
   delete m_meanPsf; 
   delete m_drm_cache;
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

   const std::vector<float> & mm = model();

   for (size_t k(0); k < ne; k++) {

      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      double w_0_sum(0.);
      double w_1_sum(0.);

      for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
	 size_t indx(k*pixels.size() + j);
         size_t indx_0 = k > 0 ? indx  - pixels.size() : indx;
	 size_t indx_1 = k < energies.size()-1 ? indx : indx - pixels.size();
	 double addend = mm.at(indx);
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
}

void SourceMap::makeProjectedMap(const ProjMap& weight_map, bool& extrapolated) {
   const std::vector<Pixel> & pixels(m_dataMap->pixels());
   std::vector<double> energy_edges;
   m_dataMap->getEnergies(energy_edges);

   std::vector<double> energies(energy_edges.size()-1,0.);
   for ( size_t ie(0); ie < energies.size(); ie++ ) {
     energies[ie] = sqrt(energy_edges[ie]*energy_edges[ie+1]);
   }

   m_model.resize(energies.size()*pixels.size());
   std::vector<Pixel>::const_iterator pixel(pixels.begin());
   extrapolated = false;
   for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
     bool in_map = weight_map.insideMap(pixel->dir());
     for (size_t k(0); k < energies.size(); k++) {      
       size_t indx(k*pixels.size() + j);
       if ( in_map ){
	 try {	   
	   m_model.at(indx) = weight_map.operator()(pixel->dir(),
						    energies[k]);
	 } catch (...) {
	   // Outside of energy bounds, set weight to 1.0 (FIXME, agree on convention)
	   m_model.at(indx) = 1.0;
	   extrapolated = true;
	 }	 
       } else {
	 // Outside of map, set weight to 1.0 (FIXME, agree on convention)
	 m_model.at(indx) = 1.0;
	 extrapolated = true;
       }
     }
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


void SourceMap::setSpectralValues(const std::vector<double>& energies) {
  if ( m_src == 0 ) return;
  FitUtils::extractSpectralVals(*m_src,energies,m_specVals);
  m_modelPars.clear();
  m_src->spectrum().getParamValues(m_modelPars);
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
  for (size_t j(0); j < parValues.size(); j++) {
    if (parValues.at(j) != m_modelPars.at(j) ) {
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
    m_drm_cache == new Drm_Cache(m_drm,*this,m_dataMap->energies());  
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

void SourceMap::readModel(const std::string& filename) {
  m_model.clear();
  m_filename = filename;
  bool ok(false);
  switch ( m_dataMap->projection().method()  ) {
  case astro::ProjBase::WCS:
    readImage(m_filename);
    ok = true;
    break;
  case astro::ProjBase::HEALPIX:
    readTable_healpix(m_filename);
    ok = true;
    break;
  default:
    break;
  }
  if ( !ok ) {
    std::string errMsg("SourceMap did not recognize CountsMapBase type at: ");
    errMsg += m_dataMap->filename();
    throw std::runtime_error(errMsg);
  }
}

void SourceMap::readImage(const std::string& sourceMapsFile) {
  FileUtils::read_fits_image_to_float_vector(sourceMapsFile,m_name,m_model);
}

void SourceMap::readTable_healpix(const std::string& sourceMapsFile) {
  FileUtils::read_healpix_table_to_float_vector(sourceMapsFile,m_name,m_model);
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
    
    status = PSFUtils::makeDiffuseMap(*m_src, *m_dataMap, m_observation.meanpsf(), m_observation.bexpmap(),
				      m_psf_config, *m_formatter, m_model);
    break;
  case Source::Point:
    m_meanPsf = PSFUtils::build_psf(*m_src,*m_dataMap,m_observation);
    status =  PSFUtils::makePointSourceMap(*m_src, *m_dataMap, m_psf_config, *m_meanPsf, *m_formatter, m_model);
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

  return status;
}

} // Likelihood
