/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.121 2016/08/05 21:04:44 echarles Exp $
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
#include "Likelihood/HealpixProjMap.h"

namespace {
   double my_acos(double mu) {
      if (mu > 1) {
         return 0;
      } else if (mu < -1) {
         return M_PI;
      } else {
         return acos(mu);
      }
   }
   
}

namespace Likelihood {

 

SourceMap::SourceMap(Source * src, const CountsMapBase * dataMap,
                     const Observation & observation, 
                     const PsfIntegConfig & psf_config,
		     const SourceMap* weights)
   : m_name(src->getName()),
     m_srcType(src->getType()),
     m_dataMap(dataMap),
     m_observation(observation),
     m_meanPsf(0),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_weights(weights),
     m_deleteDataMap(false),
     m_psf_config(psf_config),
     m_pixelOffset() {
   if (m_psf_config.verbose()) {
      m_formatter->warn() << "Generating SourceMap for " << m_name;
   }

   int status(0);
   switch ( src->srcType() ) {
   case Source::Diffuse:
      status = PSFUtils::makeDiffuseMap(*src, *dataMap, m_observation.meanpsf(), m_observation.bexpmap(),
					m_psf_config, *m_formatter, m_model);
      break;
   case Source::Point:
      m_meanPsf = PSFUtils::build_psf(*src,*dataMap,m_observation);
      status =  PSFUtils::makePointSourceMap(*src, *dataMap, m_psf_config, *m_meanPsf, *m_formatter, m_model);
      break;
   default:
      throw std::runtime_error("Unrecognized source type");
   }
   if ( status != 0 ) {
     throw std::runtime_error("SourceMap construction failed");
   }
   if (m_psf_config.verbose()) {
      m_formatter->warn() << "!" << std::endl;
   }
   applyPhasedExposureMap();
   computeNpredArray();
}

SourceMap::SourceMap(const ProjMap& weight_map,
		     const CountsMapBase * dataMap,
		     const Observation & observation,
		     bool verbose)
   : m_name("__weights__"),
     m_dataMap(dataMap),
     m_observation(observation),
     m_meanPsf(0),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_weights(0),
     m_deleteDataMap(false),
     m_psf_config(),
     m_pixelOffset() {
   if (verbose) {
      m_formatter->warn() << "Generating SourceMap for " << m_name;
   }

   bool extrapolated_weights(false);
   makeProjectedMap(weight_map,extrapolated_weights);
   // This is just to make sure that they are there.  
   // For this type of map they shouldn't be used for anything
   computeNpredArray(true);
   if (verbose) {
      m_formatter->warn() << "!" << std::endl;
   }
   if ( extrapolated_weights ) {
     m_formatter->warn() << "SourceMap::SourceMap(): "
			 << "CountsMap boundries are larger than input weights map.  "
			 << "Using weight 1 for all pixel and energies outside the input map. "
			 << std::endl;
   }
 }


SourceMap::SourceMap(const std::string & sourceMapsFile,
                     const std::string & srcName,
                     const Observation & observation,
		     const SourceMap* weights,
		     bool isWeights) 
   : m_name(srcName),
     m_dataMap(AppHelpers::readCountsMap(sourceMapsFile)),
     m_observation(observation),
     m_meanPsf(0),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_weights(weights),
     m_deleteDataMap(true),
     m_psf_config(),
     m_pixelOffset() {

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

   if ( !ok ) {
     std::string errMsg("SourceMap did not recognize CountsMapBase type at: ");
     errMsg += sourceMapsFile;
     throw std::runtime_error(errMsg);

   }
   applyPhasedExposureMap();
   computeNpredArray(isWeights);
}


SourceMap::~SourceMap() {
   if (m_deleteDataMap) {
      delete m_dataMap;
   }
   delete m_formatter;
   delete m_meanPsf; 
}


void SourceMap::computeNpredArray(bool isWeight) {
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
   size_t ne = isWeight ? energies.size() -1 : energies.size();
   // The number of weights is always number of energy bins
   size_t nw = energies.size() -1 ;
   m_npreds.clear();
   m_npreds.resize(ne);
   m_npred_weights.clear();
   m_npred_weights.resize(nw, std::make_pair<double,double>(0.,0.));

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


const std::vector<float> & SourceMap::model() const {
  return m_model;
}
  
const std::vector<double> & SourceMap::npreds() const {
  return m_npreds;
}
  
const std::vector<std::pair<double,double> > & SourceMap::npred_weights() const {
  return m_npred_weights;
}  

void SourceMap::setImage(const std::vector<float>& model) {
  if(model.size() != m_model.size())
    throw std::runtime_error("Wrong size for input model map.");

  m_model = model;
  applyPhasedExposureMap();
  computeNpredArray();
}

void SourceMap::setWeights(const SourceMap* weights) {
  m_weights = weights;
  computeNpredArray();
}

void SourceMap::readImage(const std::string& sourceMapsFile) {
  std::auto_ptr<const tip::Image>
  image(tip::IFileSvc::instance().readImage(sourceMapsFile, m_name));
  m_model.clear();
  image->get(m_model);
}

void SourceMap::readTable_healpix(const std::string& sourceMapsFile) {
  std::auto_ptr<const tip::Table> 
    table(tip::IFileSvc::instance().readTable(sourceMapsFile,m_name));

  // This is a bit tricky, basically all the data we care about
  // are in the columns called "CHANNELx"
  // Note also that tip work in lowercase
  std::vector<tip::FieldIndex_t> dataColumns;
  const tip::Table::FieldCont& colNames = table->getValidFields();
  for ( tip::Table::FieldCont::const_iterator itr = colNames.begin(); 
	itr != colNames.end(); itr++ ) {
    if ( itr->find("channel") == 0 ) { 
      dataColumns.push_back( table->getFieldIndex(*itr) );     
    } else {
      continue;
    }
  }
  
  int ncol = dataColumns.size();
  tip::Index_t nrow = table->getNumRecords();
  m_model.clear();
  m_model.reserve(ncol*nrow);

  double readVal(0.);
  for ( std::vector<tip::FieldIndex_t>::const_iterator itrData = dataColumns.begin();
	itrData != dataColumns.end(); itrData++ ) {    
    const tip::IColumn* col = table->getColumn(*itrData);
    for ( tip::Index_t irow(0); irow < nrow; irow++ ) {
      col->get(irow,readVal);
      m_model.push_back(readVal);
    }
  }
}


} // Likelihood
