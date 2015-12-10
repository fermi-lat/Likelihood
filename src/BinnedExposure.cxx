/**
 * @file BinnedExposure.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/BinnedExposure.cxx,v 1.5 2015/12/02 00:53:06 echarles Exp $
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

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/Observation.h"

namespace Likelihood {

BinnedExposure::BinnedExposure() : BinnedExposureBase() {}

BinnedExposure::BinnedExposure(const CountsMap & cmap,
                               const Observation & observation,
                               bool useEbounds,
                               const st_app::AppParGroup * pars)
  :  BinnedExposureBase(observation,useEbounds,pars){
   setMapGeometry(cmap);
   if (!useEbounds) {
      for (size_t k(0); k < m_energies.size()-1; k++) {
         m_energies[k] = std::sqrt(m_energies[k]*m_energies[k+1]);
      }
      m_energies.pop_back();
      m_naxes[2] -= 1;
   }
   computeMap();
}

BinnedExposure::BinnedExposure(const std::vector<double> & energies,
                               const Observation & observation,
                               const st_app::AppParGroup * pars) 
  :  BinnedExposureBase(energies,observation,pars){
   if (pars) {
      setMapGeometry(*pars);
   } else {
      setMapGeometry();
   }
   computeMap();
}

BinnedExposure::BinnedExposure(const std::string & filename) 
  :  BinnedExposureBase(filename){

   m_proj = new astro::SkyProj(filename);

   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(filename, ""));

   m_exposureMap.clear();
   image->get(m_exposureMap);

   image->getHeader().getKeyword("CDELT1",m_cdelt[0]);
   image->getHeader().getKeyword("CDELT2",m_cdelt[1]);
   image->getHeader().getKeyword("CRVAL1",m_crval[0]);
   image->getHeader().getKeyword("CRVAL2",m_crval[1]);
   image->getHeader().getKeyword("CRPIX1",m_crpix[0]);
   image->getHeader().getKeyword("CRPIX2",m_crpix[1]);
   m_naxes = image->getImageDimensions();

   // EAC, if we have an all-sky projection it saves
   // us the trouble of figuring out if stuff is actually 
   // in the map or not
   if ( std::fabs( m_cdelt[0]*m_naxes[0] ) >= 360. &&
	std::fabs( m_cdelt[1]*m_naxes[1] ) >= 180. ) {
     m_allSky = true;
   }

}

BinnedExposure::~BinnedExposure() {
}

double BinnedExposure::operator()(double energy, double ra, double dec) const {

   if ( m_proj == 0 ) {
     throw std::runtime_error("BinnedExposure::operator() "
			      "No projection loaded");
   }
  
   std::vector<double>::const_iterator ie = findNearest(m_energies, energy);
   unsigned int k = ie - m_energies.begin();

   std::pair<double, double> pixel;
   st_facilities::Util::skyDir2pixel(*m_proj, astro::SkyDir(ra, dec),
                                     pixel.first, pixel.second);

   int i = static_cast<int>(pixel.first - 1);
   int j = static_cast<int>(pixel.second - 1);

   unsigned int indx = (k*m_naxes.at(1) + j)*m_naxes.at(0) + i;

   bool within_bounds = (i >= 0 && i < m_naxes[0] &&
                         j >= 0 && j < m_naxes[1] &&
                         k >= 0 && k < m_energies.size());

   if (m_enforce_boundaries && !within_bounds) {
      throw std::runtime_error("Request for exposure at a sky position that "
                               "is outside of the map boundaries.");
   }

   try {
      return m_exposureMap.at(indx);
   } catch (std::out_of_range &) {
      // Range check performed already, so do nothing and return 0.
   }
   return 0;
}

void BinnedExposure::setMapGeometry(const CountsMap & cmap) {
   m_proj_name = cmap.proj_name();
   m_crpix[0] = cmap.crpix1();
   m_crpix[1] = cmap.crpix2();
   m_crval[0] = cmap.crval1();
   m_crval[1] = cmap.crval2();
   m_cdelt[0] = cmap.cdelt1();
   m_cdelt[1] = cmap.cdelt2();
   m_crota2 = cmap.crota2();
   m_naxes.resize(3, 0);
   m_naxes[0] = cmap.naxis1();
   m_naxes[1] = cmap.naxis2();
   m_naxes[2] = cmap.energies().size();
   m_energies = cmap.energies();
   m_isGalactic = cmap.isGalactic();
   if ( std::fabs( m_cdelt[0]*m_naxes[0] ) >= 360. &&
	std::fabs( m_cdelt[1]*m_naxes[1] ) >= 180. ) {
     m_allSky = true;
   }
}

void BinnedExposure::setMapGeometry(const st_app::AppParGroup & pars) {
   m_naxes.resize(3, 0);
   m_naxes[0] = pars["nxpix"];
   m_naxes[1] = pars["nypix"];
   m_naxes[2] = m_energies.size();
   double binsz = pars["binsz"];
   m_cdelt[0] = -binsz;
   m_cdelt[1] = binsz;
   m_crval[0] = pars["xref"];
   m_crval[1] = pars["yref"];
   m_crota2 = pars["axisrot"];
   std::string proj_name = pars["proj"];
   m_proj_name = proj_name;
   m_crpix[0] = m_naxes[0]/2. + 0.5;
   m_crpix[1] = m_naxes[1]/2. + 0.5;
   std::string coordsys = pars["coordsys"];
   m_isGalactic = (coordsys == "GAL");
   if ( std::fabs( m_cdelt[0]*m_naxes[0] ) >= 360. &&
	std::fabs( m_cdelt[1]*m_naxes[1] ) >= 180. ) {
     m_allSky = true;
   }
}

void BinnedExposure::setMapGeometry() {
   m_proj_name = "CAR";
   m_isGalactic = false;
   m_naxes.resize(3, 0);
   m_naxes[0] = 360;
   m_naxes[1] = 180;
   m_naxes[2] = m_energies.size();
   m_crpix[0] = m_naxes[0]/2. + 0.5;
   m_crpix[1] = m_naxes[1]/2. + 0.5;
   m_crval[0] = 180.;
   m_crval[1] = 0;
   m_cdelt[0] = -1;
   m_cdelt[1] = 1;
   m_crota2 = 0;
   m_allSky = true;
}

void BinnedExposure::computeMap() {
   m_proj = new astro::SkyProj(m_proj_name, &m_crpix[0], &m_crval[0],
                               &m_cdelt[0], m_crota2, m_isGalactic);

   m_exposureMap.resize(m_naxes.at(0)*m_naxes.at(1)*m_energies.size(), 0);
   long iter(0);
   st_stream::StreamFormatter formatter("BinnedExposure", "computeMap", 2);
   formatter.warn() << "Computing binned exposure map";

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

   long npix(m_naxes[0]*m_naxes[1]);
   for (int j = 0; j < m_naxes.at(1); j++) {
      for (int i = 0; i < m_naxes.at(0); i++, iter++) {
         if (npix > 20 && (iter % (npix/20)) == 0) {
            formatter.warn() << ".";
         }
         // std::pair<double, double> coord = m_proj->pix2sph(i + 1, j + 1);
         // astro::SkyDir dir(coord.first, coord.second, coordsys);
         astro::SkyDir dir;
         try {
            st_facilities::Util::pixel2SkyDir(*m_proj, i + 1, j + 1, dir);
         } catch (...) {
            // The astro::SkyProj class throws a SkyProjException
            // here, but SkyProjException is annoyingly defined in
            // SkyProj.cxx
            // http://www-glast.stanford.edu/cgi-bin/viewcvs/astro/src/SkyProj.cxx?revision=1.27&view=markup
            // so that client code cannot catch it directly. Amazing.
            continue;
         }
                                           
         for (unsigned int k = 0; k < m_energies.size(); k++) {
            unsigned int indx = (k*m_naxes.at(1) + j)*m_naxes.at(0) + i;
            std::map<unsigned int, irfInterface::Irfs *>::const_iterator 
               resp = m_observation->respFuncs().begin();
            for (; resp != m_observation->respFuncs().end(); ++resp) {
               int evtType = resp->second->irfID();
               // Aeff aeff(m_energies[k], evtType, *m_observation, m_costhmin,
               //           m_costhmax);
               Aeff & aeff = *aeffs[std::make_pair(k, evtType)];
               m_exposureMap.at(indx)
                  +=m_observation->expCube().value(dir, aeff, m_energies.at(k));
            }
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

void BinnedExposure::writeOutput(const std::string & filename) const {
   std::remove(filename.c_str());

   std::string ext("PRIMARY");
   tip::IFileSvc::instance().appendImage(filename, ext, m_naxes);
   tip::Image * image = tip::IFileSvc::instance().editImage(filename, ext);

   image->set(m_exposureMap);

   tip::Header & header(image->getHeader());

   header["TELESCOP"].set("GLAST");
   header["INSTRUME"].set("LAT");
   astro::JulianDate current_time = st_facilities::Util::currentTime();
   header["DATE"].set(current_time.getGregorianDate());
   header["DATE-OBS"].set("");
   header["DATE-END"].set("");

   header["CRVAL1"].set(m_crval[0]);
   double crpix1 = static_cast<double>(m_naxes[0] + 1)/2.;
   header["CRPIX1"].set(crpix1);
   header["CDELT1"].set(m_cdelt[0]);

   header["CRVAL2"].set(m_crval[1]);
   double crpix2 = static_cast<double>(m_naxes[1] + 1)/2.;
   header["CRPIX2"].set(crpix2);
   header["CDELT2"].set(m_cdelt[1]);

   if (m_isGalactic) {
      header["CTYPE1"].set("GLON-" + m_proj_name);
      header["CTYPE2"].set("GLAT-" + m_proj_name);
   } else {
      header["CTYPE1"].set("RA---" + m_proj_name);
      header["CTYPE2"].set("DEC--" + m_proj_name);
   }

   int nee = m_energies.size();
   header["CRVAL3"].set(log(m_energies.at(0)));
   header["CRPIX3"].set(1);
   if (nee == 1) {
      header["CDELT3"].set(0);
   } else {
      header["CDELT3"].set(log(m_energies.at(nee-1)/m_energies.at(0))/(nee-1));
   }
   header["CTYPE3"].set("log_Energy");

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
