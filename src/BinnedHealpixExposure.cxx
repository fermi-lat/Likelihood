/**
 * @file BinnedHealpixExposure.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies. Healpix version
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedHealpixExposure.cxx,v 1.48 2013/01/10 08:56:12 sfegan Exp $
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

#include "astro/SkyProj.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedHealpixExposure.h"
#include "Likelihood/Observation.h"

#include "healpix/Healpix.h"

namespace Likelihood {

BinnedHealpixExposure::BinnedHealpixExposure(const evtbin::HealpixMap & cmap,
					     const Observation & observation,
					     bool useEbounds,
					     const st_app::AppParGroup * pars)
  : BinnedExposure(), m_observation(&observation) {
  setMapGeometry(cmap);
  if (!useEbounds) {
    for (size_t k(0); k < m_energies.size()-1; k++) {
      m_energies[k] = std::sqrt(m_energies[k]*m_energies[k+1]);
    }
    m_energies.pop_back();
    m_naxes[1] -= 1;
  }
  if (pars) {
    setCosThetaBounds(*pars);
  }
  std::string filename=(*pars)["outfile"];
  computeHealpixMap(cmap,filename);
}
  
BinnedHealpixExposure::~BinnedHealpixExposure() {
//   delete m_proj;
}

  void BinnedHealpixExposure::setMapGeometry(const evtbin::HealpixMap & cmap) {
   m_naxes.resize(2, 0);
   m_naxes[0] = cmap.hpx_binner()->getNumBins();
   m_naxes[1] = cmap.energies().size();
   m_energies = cmap.energies();
   m_isGalactic = cmap.isGalactic();

}

void BinnedHealpixExposure::computeHealpixMap(const evtbin::HealpixMap & cmap, const std::string filename) {
   m_exposureMap.resize(m_naxes.at(0)*m_energies.size(), 0);
   long iter(0);
   st_stream::StreamFormatter formatter("BinnedHealpixExposure", "computeHealpixMap", 2);
   formatter.warn() << "Computing Healpix binned exposure map";

   // Storing Aeff objects for use within loops over sky pixels.
   std::map<std::pair<unsigned int, int>, Aeff *> aeffs;
   for (unsigned int k(0); k < m_energies.size(); k++) {
      std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
         resp = m_observation->respFuncs().begin();
      for (; resp != m_observation->respFuncs().end(); ++resp) {
         int evtType = resp->second->irfID();
         aeffs[std::make_pair(k, evtType)] =  
            new Aeff(m_energies[k], evtType, *m_observation, m_costhmin,
                     m_costhmax);
      }
   }

   int npix(m_naxes[0]);
   astro::SkyDir::CoordSystem coordsys=cmap.isGalactic()?astro::SkyDir::GALACTIC:astro::SkyDir::EQUATORIAL;
   healpix::Healpix::Ordering order=cmap.ordering()=="RING"?healpix::Healpix::RING:healpix::Healpix::NEST;
   healpix::Healpix hp(cmap.nside(),order,coordsys);
   for (int i = 0; i < npix; i++, iter++) {
     if (npix > 20 && (iter % (npix/20)) == 0) {
       formatter.warn() << ".";
     }
     healpix::Healpix::Pixel pix(i,hp);
     astro::SkyDir dir=pix();

     for (unsigned int k = 0; k < m_energies.size(); k++) {
       unsigned int indx = k*m_naxes.at(0) + i;
       std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
	 resp = m_observation->respFuncs().begin();
       for (; resp != m_observation->respFuncs().end(); ++resp) {
	 int evtType = resp->second->irfID();
	 Aeff & aeff = *aeffs[std::make_pair(k, evtType)];
	 m_exposureMap.at(indx)
	   +=m_observation->expCube().value(dir, aeff, m_energies.at(k));
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

   //I need to start writing the outfile here to get the healpix info in the header
   std::string ext="HPXEXPOSURES";
   std::remove(filename.c_str());
   tip::IFileSvc::instance().appendTable(filename, ext);
   tip::Table * image = tip::IFileSvc::instance().editTable(filename, ext);
   tip::Header & header(image->getHeader());
   header["PIXTYPE"].set("HEALPIX");
   header["ORDERING"].set(cmap.ordering());
   header["ORDER"].set(cmap.order());
   header["NSIDE"].set(cmap.nside());
   header["FIRSTPIX"].set(0);
   header["LASTPIX"].set(12*cmap.nside()*cmap.nside()-1);
   header["NBRBINS"].set((m_naxes.at(1)==0?1:m_naxes.at(1)));
   delete image;
}

  void BinnedHealpixExposure::writeOutput(const std::string & filename) const {
    std::string ext="";
    ext="HPXEXPOSURES";
    //tip::IFileSvc::instance().appendTable(filename, ext);
    tip::Table * image = tip::IFileSvc::instance().editTable(filename, ext);
    for (long e_index = 0; e_index != (m_naxes.at(1)==0?1:m_naxes.at(1)); ++e_index) {
      std::ostringstream e_channel;
      e_channel<<"ENERGY"<<e_index+1;
      //create new column
      image->appendField(e_channel.str(), std::string("D"));
     //loop over rows and fill healpix values in
      tip::Table::Iterator table_itor = image->begin();
      for(long hpx_index = 0;hpx_index != m_naxes.at(0);++hpx_index,++table_itor) {
	(*table_itor)[e_channel.str()].set(m_exposureMap[hpx_index+e_index*m_naxes.at(0)]);
      }
    }
    //these should perhaps stay in the Primary header
    tip::Header & header(image->getHeader());
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
