/**
 * @file MapCubeFunction.cxx
 * @brief Encapsulation of 3D FITS image of a diffuse source with 
 * position-dependent spectral variation.
 * @author jchiang
 *
 * $Header$
 */

#include <stdexcept>

#include "Likelihood/FitsImage.h"
#include "Likelihood/MapCubeFunction.h"
#include "Likelihood/SkyDirArg.h"

namespace {
   bool reverse_cmp(double x, double y) {
      return x > y;
   }
}

namespace Likelihood {

MapCubeFunction::MapCubeFunction(const std::string & fitsFile) {
   FitsImage fitsImage(fitsFile);

   std::vector<std::string> axisNames;
   fitsImage.getAxisNames(axisNames);
   if (axisNames[0].find_first_of("RA") != std::string::npos) {
      m_coordSys = "Equatorial";
   } else if (axisNames[0].find_first_of("GLON") != std::string::npos) {
      m_coordSys = "Galactic";
   } else {
      std::ostringstream message;
      message << "Likelihood::MapCubeFunction\n"
              << "Unrecognized coordinate system in " << fitsFile << ".\n"
              << "Axis names: ";
      for (unsigned int i = 0; i < axisNames.size(); i++) {
         message << axisNames.at(i) << "  ";
      }
      throw std::runtime_error(message.str());
   }

   fitsImage.getAxisVector(0, m_lon);
   fitsImage.getAxisVector(1, m_lat);
   if (m_lon.front() < m_lon.back()) {
      m_lonMin = m_lon.front();
      m_lonMax = m_lon.back();
   } else {
      m_lonMin = m_lon.back();
      m_lonMax = m_lon.front();
   }
   if (m_lat.front() < m_lat.back()) {
      m_latMin = m_lat.front();
      m_latMax = m_lat.back();
   } else {
      m_latMin = m_lat.back();
      m_latMax = m_lat.front();
   }

   fitsImage.getImageData(m_image);

   readEnergyVector(fitsFile);

   setMaxNumParams(1);
   addParam("Normalization", 1, true);

   m_funcType = Addend;
   m_argType = "";
   m_genericName = "MapCubeFunction";
}

double MapCubeFunction::value(optimizers::Arg & x) const {
   SkyDirArg & dir = dynamic_cast<SkyDirArg&>(x);
   double lonValue;
   double latValue;
   double energy = dir.energy();
   if (m_coordSys == "Equatorial") {
      lonValue = dir().ra();
      latValue = dir().dec();
   } else {
      lonValue = dir().l();
      latValue = dir().b();
   }
   if (lonValue < m_lonMin || lonValue > m_lonMax ||
       latValue < m_latMin || latValue > m_latMax || 
       energy < m_energies.front() || energy > m_energies.back()) {
      return 0;
   }

   int i = findIndex(m_lon, lonValue);
   int j = findIndex(m_lat, latValue);
   int k = findIndex(m_energies, energy);

   int indx = k*m_lon.size()*m_lat.size() + j*m_lon.size() + i;
   return m_image.at(indx)*getParamValue("Normalization");
}

int MapCubeFunction::findIndex(std::vector<double> xx, double x) {
   std::vector<double>::iterator it;
   if (xx.front() < xx.back()) {
      return std::upper_bound(xx.begin(), xx.end(), x) - xx.begin() - 1;
   }
   return std::upper_bound(xx.begin(), xx.end(), x, ::reverse_cmp) 
      - xx.begin() - 1;
}

void MapCubeFunction::readEnergyVector(const std::string & fitsFile) {

   std::string routineName("readEnergyVector");

   int hdu = FitsImage::findHdu(fitsFile, "ENERGIES");
   
   int status(0);
   fitsfile * fptr = 0;

   fits_open_file(&fptr, fitsFile.c_str(), READONLY, &status);
   FitsImage::fitsReportError(status, routineName);

   int hdutype(0);
   fits_movabs_hdu(fptr, hdu, &hdutype, &status);
   FitsImage::fitsReportError(status, routineName);

   long nrows(0);
   fits_get_num_rows(fptr, &nrows, &status);
   FitsImage::fitsReportError(status, routineName);

   FitsImage::readColumn(fptr, "Energy", m_energies);

   fits_close_file(fptr, &status);
   FitsImage::fitsReportError(status, routineName);
}

}
