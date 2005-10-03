/**
 * @file BinnedExposure.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedExposure.cxx,v 1.10 2005/07/13 17:17:43 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/EquinoxRotation.h"
#include "Likelihood/Observation.h"

#include "Verbosity.h"

namespace Likelihood {

BinnedExposure::BinnedExposure() : m_observation(0) {}

BinnedExposure::BinnedExposure(const std::vector<double> & energies,
                               const Observation & observation) 
   : m_observation(&observation), m_energies(energies) {
   computeMap();
}

BinnedExposure::BinnedExposure(const std::string & filename) {
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(filename, ""));

   m_exposureMap.clear();
   image->get(m_exposureMap);

   const tip::Header & header = image->getHeader();

   /// @bug tip doesn't allow direct access to NAXIS[1-3] keywords
   /// through the Header, so we are forced to use
   /// getImageDimensions() method.
   std::vector<tip::PixOrd_t> dims = image->getImageDimensions();

   unsigned int nra = dims[0];
   double ra0;
   double rastep;
   double refpix1;
   header["CRVAL1"].get(ra0);
   header["CDELT1"].get(rastep);
   header["CRPIX1"].get(refpix1);
   m_ras.resize(nra);
   for (unsigned int i = 0; i < nra; i++) {
      m_ras.at(i) = (i - refpix1)*rastep + ra0;
   }
   
   unsigned int ndec = dims[1];
   double dec0;
   double decstep;
   double refpix2;
   header["CRVAL2"].get(dec0);
   header["CDELT2"].get(decstep);
   header["CRPIX2"].get(refpix2);
   m_decs.resize(ndec);
   for (unsigned int i = 0; i < ndec; i++) {
      m_decs.at(i) = (i - refpix2)*decstep + dec0;
   }

   std::auto_ptr<const tip::Table>
      energies(tip::IFileSvc::instance().readTable(filename, "Energies"));

   m_energies.clear();
   tip::Table::ConstIterator it = energies->begin();
   tip::ConstTableRecord & row = *it;
   for ( ; it != energies->end(); ++it) {
      double value;
      row["Energy"].get(value);
      m_energies.push_back(value);
   }
}

double BinnedExposure::operator()(double energy, double ra, double dec) const {
   std::vector<double>::const_iterator ie = std::find(m_energies.begin(),
                                                      m_energies.end(),
                                                      energy);
   if (ie == m_energies.end()) {
      std::ostringstream what;
      what << "BinnedExposure::operator(): The energy " << energy 
           << " is not available.\nHere are the relevant energies:\n";
      for (unsigned int i = 0; i < m_energies.size(); i++) {
         what << m_energies.at(i) << "\n";
      }
      throw std::runtime_error(what.str());
   }
   unsigned int k = ie - m_energies.begin();
   if (ra < 0) ra += 360.; 
   if (ra > 360.) ra = fmod(ra, 360.);
   unsigned int i = findIndex(m_ras.begin(), m_ras.end(), ra);
   unsigned int j = findIndex(m_decs.begin(), m_decs.end(), dec);
   unsigned int indx = (k*m_decs.size() + j)*m_ras.size() + i;
   return m_exposureMap.at(indx);
}

void BinnedExposure::
getRotatedImage(double energy, 
                const std::vector<double> & lons, 
                const std::vector<double> & lats,
                const EquinoxRotation & rot,
                std::vector< std::vector<double> > & image) const {
   image.clear();
   image.reserve(lons.size());
   for (unsigned int i = 0; i < lons.size(); i++) {
      std::vector<double> row;
      row.reserve(lats.size());
      for (unsigned int j = 0; j < lats.size(); j++) {
         astro::SkyDir mapDir(lons.at(i), lats.at(j));
         astro::SkyDir trueDir;
         rot.do_rotation(mapDir, trueDir, false);
         row.push_back((*this)(energy, trueDir.ra(), trueDir.dec()));
      }
      image.push_back(row);
   }
}

unsigned int 
BinnedExposure::findIndex(std::vector<double>::const_iterator begin,
                          std::vector<double>::const_iterator end,
                          double value) const {
   std::vector<double>::const_iterator it 
      = std::upper_bound(begin, end, value);
   if (it == end) {
      throw std::range_error("BinnedExposure::findIndex: out of range value.");
   }
   int indx = it - begin - 1;
   return indx;
}

void BinnedExposure::computeMap() {
   linearArray(0., 360., 360, m_ras);
   linearArray(-90., 90., 180, m_decs);

   m_exposureMap.resize(m_ras.size()*m_decs.size()*m_energies.size(), 0);
   int iter(0);
   if (print_output()) {
      std::cerr << "Computing binned exposure map";
   }
   for (unsigned int j = 0; j < m_decs.size(); j++) {
      for (unsigned int i = 0; i < m_ras.size(); i++) {
         if ( print_output() && 
              (iter % ((m_decs.size()*m_ras.size())/20)) == 0 ) {
            std::cerr << ".";
         }
         astro::SkyDir dir(m_ras[i], m_decs[j], astro::SkyDir::EQUATORIAL);
         for (unsigned int k = 0; k < m_energies.size(); k++) {
            unsigned int indx 
               = k*m_ras.size()*m_decs.size() + j*m_ras.size() + i;
            for (int evtType = 0; evtType < 2; evtType++) {
               Aeff aeff(m_energies[k], evtType, *m_observation);
               m_exposureMap.at(indx)
                  += m_observation->expCube().value(dir, aeff);
            }
         }
         iter++;
      }
   }
   if (print_output()) {
      std::cerr << "!" << std::endl;
   }
}

void BinnedExposure::linearArray(double xmin, double xmax, unsigned int npts,
                                 std::vector<double> &xx) const {
   double xstep = (xmax - xmin)/(npts - 1.);
   for (unsigned int i = 0; i < npts; i++) {
      xx.push_back(i*xstep + xmin);
   }
}

double BinnedExposure::Aeff::s_phi(0);

double BinnedExposure::Aeff::operator()(double cosTheta) const {
   double inclination = acos(cosTheta)*180./M_PI;
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = m_observation.respFuncs().begin();
   for ( ; respIt != m_observation.respFuncs().end(); ++respIt) {
      if (respIt->second->irfID() == m_evtType) {
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double aeff_val = aeff->value(m_energy, inclination, s_phi);
         return aeff_val;
      }
   }
   return 0;
}

void BinnedExposure::writeOutput(const std::string & filename) const {
   std::remove(filename.c_str());

   std::vector<long> dims(3);
   dims.at(0) = m_ras.size();
   dims.at(1) = m_decs.size();
   dims.at(2) = m_energies.size();

   std::string ext("PRIMARY");
   tip::IFileSvc::instance().appendImage(filename, ext, dims);
   tip::Image * image = tip::IFileSvc::instance().editImage(filename, ext);

   image->set(m_exposureMap);

   tip::Header & header(image->getHeader());

   header["TELESCOP"].set("GLAST");
   header["INSTRUME"].set("LAT SIMULATION");
   header["DATE-OBS"].set("");
   header["DATE-END"].set("");

   header["CRVAL1"].set(m_ras[0]);
   header["CRPIX1"].set(0);
   header["CDELT1"].set(m_ras[1] - m_ras[0]);
   header["CTYPE1"].set("RA---CAR");

   header["CRVAL2"].set(m_decs[0]);
   header["CRPIX2"].set(0);
   header["CDELT2"].set(m_decs[1] - m_decs[0]);
   header["CTYPE2"].set("DEC--CAR");

   int nee = m_energies.size();
   header["CRVAL3"].set(log(m_energies.at(0)));
   header["CRPIX3"].set(1);
   header["CDELT3"].set(log(m_energies.at(nee-1)/m_energies.at(0))/(nee-1));
   header["CTYPE3"].set("log_Energy");

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
