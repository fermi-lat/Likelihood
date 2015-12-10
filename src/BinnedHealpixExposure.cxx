/**
 * @file BinnedHealpixExposure.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies. Healpix version
 * @author E. Charles, J. Cohen-Tanugi
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/BinnedHealpixExposure.cxx,v 1.5 2015/12/02 00:53:06 echarles Exp $
 */

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

//#include "astro/SkyProj.h"

#include "st_facilities/Util.h"

#include "astro/HealpixProj.h"

#include "Likelihood/BinnedHealpixExposure.h"
#include "Likelihood/Observation.h"


namespace Likelihood {

BinnedHealpixExposure::BinnedHealpixExposure(const evtbin::HealpixMap & cmap,
					     const Observation & observation,
					     bool useEbounds,
					     const st_app::AppParGroup * pars)
  : BinnedExposureBase(observation,useEbounds,pars) {
  setMapGeometry(cmap);
  if (!useEbounds) {
    for (size_t k(0); k < m_energies.size()-1; k++) {
      m_energies[k] = std::sqrt(m_energies[k]*m_energies[k+1]);
    }
    m_energies.pop_back();
  }
  if (pars) {
    setCosThetaBounds(*pars);
  }
  computeHealpixMap(cmap);
}
  
BinnedHealpixExposure::BinnedHealpixExposure(const std::string& filename)
  : BinnedExposureBase(filename),
    m_healpixProj(new astro::HealpixProj(filename,std::string("HPXEXPOSURES"))){
  setMapGeometry();

  std::auto_ptr<const tip::Table> 
    table(tip::IFileSvc::instance().readTable(filename, std::string("HPXEXPOSURES")));
  
  // This is a bit tricky, basically all the data we care about
  // are in the columns called "ENERGYx"
  // Note also that tip work in lowercase
  std::vector<tip::FieldIndex_t> dataColumns;
  const tip::Table::FieldCont& colNames = table->getValidFields();
  for ( tip::Table::FieldCont::const_iterator itr = colNames.begin(); 
        itr != colNames.end(); itr++ ) {
    if ( itr->find("energy") == 0 ) { 
      dataColumns.push_back( table->getFieldIndex(*itr) );     
    } else {
      continue;
    }
  }

  int ncol = dataColumns.size();
  tip::Index_t nrow = table->getNumRecords();

  m_exposureMap.resize(ncol);

  std::vector<float> oneColumn(nrow,0);
  // This keeps track of the energy slice we are filling.  
  // Starts with slice 0
  int iPlane(0);    
  for ( std::vector<tip::FieldIndex_t>::const_iterator itrData = dataColumns.begin();
        itrData != dataColumns.end(); itrData++, iPlane++ ) {    
    ExposurePlane_t& expPlane = m_exposureMap[iPlane];
    expPlane.SetNside(m_healpixProj->healpix().Nside(),m_healpixProj->healpix().Scheme());
    const tip::IColumn* col = table->getColumn(*itrData);
    for ( tip::Index_t irow(0); irow < nrow; irow++ ) {
      col->get(irow,expPlane[irow]);
    }    
  }
}
  
BinnedHealpixExposure::~BinnedHealpixExposure() {
}

double BinnedHealpixExposure::operator()(double energy, double ra, double dec) const {
  std::vector<double>::const_iterator ie = BinnedExposureBase::findNearest(m_energies, energy);
  unsigned int k = ie - m_energies.begin();
  std::pair<double, double> pixel;
  st_facilities::Util::skyDir2pixel(*m_proj, astro::SkyDir(ra, dec),
				    pixel.first, pixel.second);

  int i = static_cast<int>(pixel.first);
  bool within_bounds = (i >= 0 && i < m_healpixProj->healpix().Npix() &&
			k >= 0 && k < m_energies.size());

  if (m_enforce_boundaries && !within_bounds) {
    throw std::runtime_error("Request for exposure at a sky position that "
			     "is outside of the map boundaries.");
  }
  try {
    return m_exposureMap[k][i];
  } catch (std::out_of_range &) {
    // Range check performed already, so do nothing and return 0.
  }
  return 0.;
}


void BinnedHealpixExposure::setMapGeometry(const evtbin::HealpixMap & cmap) {
  m_energies = cmap.energies();
  m_isGalactic = cmap.isGalactic();  
  m_healpixProj = new astro::HealpixProj(cmap.nside(),cmap.scheme(),SET_NSIDE,m_isGalactic);
  m_proj = m_healpixProj;  
  m_allSky = true;
}

void BinnedHealpixExposure::setMapGeometry() {
  m_proj = m_healpixProj;  
  m_isGalactic = m_proj->isGalactic();
  m_allSky = true;
}


void BinnedHealpixExposure::computeHealpixMap(const evtbin::HealpixMap & cmap) {
  m_exposureMap.resize(m_energies.size());
  st_stream::StreamFormatter formatter("BinnedHealpixExposure", "computeHealpixMap", 2);
  formatter.warn() << "Computing Healpix binned exposure map";
  // Storing Aeff objects for use within loops over sky pixels.
  std::map<std::pair<unsigned int, int>, Aeff *> aeffs;
  for (unsigned int k(0); k < m_energies.size(); k++) {
    m_exposureMap[k].SetNside(m_healpixProj->healpix().Nside(),m_healpixProj->healpix().Scheme());
    std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
      resp = m_observation->respFuncs().begin();
    for (; resp != m_observation->respFuncs().end(); ++resp) {
      int evtType = resp->second->irfID();
      aeffs[std::make_pair(k, evtType)] =  
	new Aeff(m_energies[k], evtType, *m_observation, m_costhmin,
		 m_costhmax);
    }
  }
  int npix = m_healpixProj->healpix().Npix();
  for (int i = 0; i < npix; i++ ) {
    if (npix > 20 && (i % (npix/20)) == 0) {
      formatter.warn() << ".";
    }
    astro::SkyDir dir(i,0.,*m_healpixProj);
    for (unsigned int k = 0; k < m_energies.size(); k++) {
      m_exposureMap[k][i] = 0.;      
      std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
	resp = m_observation->respFuncs().begin();
      
      for (; resp != m_observation->respFuncs().end(); ++resp) {
	int evtType = resp->second->irfID();
	Aeff & aeff = *aeffs[std::make_pair(k, evtType)];
	double val = m_observation->expCube().value(dir, aeff, m_energies.at(k));
	m_exposureMap[k][i] += val;
      }
    }
  }
  // Release memory of Aeff map.
  for (unsigned int k(0); k < m_energies.size(); k++) {
    std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
      resp = m_observation->respFuncs().begin();
    for (; resp != m_observation->respFuncs().end(); ++resp) {
      int evtType = resp->second->irfID();
      delete aeffs[std::make_pair(k, evtType)];
    }
  }  
  formatter.warn() << "!" << std::endl;
}

void BinnedHealpixExposure::writeOutput(const std::string & filename) const {
  std::remove(filename.c_str());

  std::string ext="HPXEXPOSURES";
  tip::IFileSvc::instance().appendTable(filename, ext);
  tip::Table * image = tip::IFileSvc::instance().editTable(filename, ext);
  tip::Header & header(image->getHeader());

  m_healpixProj->setKeywords(header);
  for (long e_index = 0; e_index != m_exposureMap.size(); e_index++ ) {
    std::ostringstream e_channel;
    e_channel<<"ENERGY"<<e_index+1;
    //create new column
    image->appendField(e_channel.str(), std::string("D"));
    //loop over rows and fill healpix values in
    tip::Table::Iterator table_itor = image->begin();
    for(long hpx_index = 0;hpx_index != m_healpixProj->healpix().Npix();++hpx_index,++table_itor) {
      (*table_itor)[e_channel.str()].set(m_exposureMap[e_index][hpx_index]);
    }
  }
  //these should perhaps stay in the Primary header
  header["TELESCOP"].set("GLAST");
  header["INSTRUME"].set("LAT");
  astro::JulianDate current_time = st_facilities::Util::currentTime();
  header["DATE"].set(current_time.getGregorianDate());
  header["DATE-OBS"].set("");
  header["DATE-END"].set("");
  if (!m_isGalactic) {
    header["EQUINOX"].set(2000.0);
    header["RADECSYS"].set("FK5");
  }
  delete image;

  ext = "ENERGIES";
  tip::IFileSvc::instance().appendTable(filename, ext);
  tip::Table * table = tip::IFileSvc::instance().editTable(filename, ext);
  table->appendField("Energy", "1D");
  table->setNumRecords(m_energies.size());
    
  tip::Table::Iterator row = table->begin();
  tip::Table::Record & record = *row;
    
  std::vector<double>::const_iterator energy = m_energies.begin();
  for ( ; energy != m_energies.end(); ++energy, ++row) {
    record["Energy"].set(*energy);
  }
    
  delete table;

}

} // namespace Likelihood
