/**
 * @file WeightMap.cxx
 * @brief Pixel-by-pixel weights for binned likelihood.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/WeightMap.cxx,v 1.123 2016/09/13 19:26:23 echarles Exp $
 */

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
#include "Likelihood/FitUtils.h"
#include "Likelihood/FileUtils.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/PSFUtils.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/Source.h"
#define ST_DLL_EXPORTS
#include "Likelihood/WeightMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/SpatialFunction.h"

#include "Likelihood/WcsMap2.h"
#include "Likelihood/HealpixProjMap.h"



namespace Likelihood {

  WeightMap::WeightMap(const ProjMap& weight_map,
		       const CountsMapBase * dataMap,
		       const Observation & observation,
		       bool verbose)
    : m_srcMapFile(""),
      m_dataMap(dataMap),
      m_observation(observation),
      m_formatter(new st_stream::StreamFormatter("WeightMap", "", 2)) {
    if (verbose) {
      m_formatter->warn() << "Generating WeightMap ";
    }
    
    bool extrapolated_weights(false);
    makeProjectedMap(weight_map,extrapolated_weights);
    // This is just to make sure that they are there.  
    // For this type of map they shouldn't be used for anything
    if (verbose) {
      m_formatter->warn() << "!" << std::endl;
    }
    if ( extrapolated_weights ) {
      m_formatter->warn() << "WeightMap::WeightMap(): "
			  << "CountsMap boundries are larger than input weights map.  "
			  << "Using weight 1 for all pixel and energies outside the input map. "
			  << std::endl;
    }
  }


  WeightMap::WeightMap(const std::string & sourceMapsFile,
		       const CountsMapBase * dataMap,
		       const Observation & observation) 
    : m_srcMapFile(sourceMapsFile),
      m_dataMap(dataMap),
      m_observation(observation),
      m_formatter(new st_stream::StreamFormatter("WeightMap", "", 2)){

    m_model.clear();
    bool ok(false);
    switch ( m_dataMap->projection().method()  ) {
    case astro::ProjBase::WCS:
      readImage(sourceMapsFile);
      ok = true;
      break;
    case astro::ProjBase::HEALPIX:
      readTable_healpix(sourceMapsFile);
      ok = true;
      break;
    default:
      break;
    }
  }


  WeightMap::~WeightMap() {
    delete m_formatter;
  }


  void WeightMap::makeProjectedMap(const ProjMap& weight_map, bool& extrapolated) {
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


  void WeightMap::setImage(const std::vector<float>& model) {
    if(model.size() != m_model.size())
      throw std::runtime_error("Wrong size for input model map.");
    
    m_model = model;
  }

  void WeightMap::readImage(const std::string& sourceMapsFile) {
    FileUtils::read_fits_image_to_float_vector(sourceMapsFile,"__weights__",m_model);
  }

  void WeightMap::readTable_healpix(const std::string& sourceMapsFile) {
    FileUtils::read_healpix_table_to_float_vector(sourceMapsFile,"__weights__",m_model);
  }

} // Likelihood
