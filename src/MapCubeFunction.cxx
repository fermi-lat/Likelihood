/**
 * @file MapCubeFunction.cxx
 * @brief Encapsulation of 3D FITS image of a diffuse source with 
 * position-dependent spectral variation.
 * @author jchiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MapCubeFunction.cxx,v 1.4 2005/02/15 07:04:41 jchiang Exp $
 */

#include <algorithm>
#include <stdexcept>

#include "Likelihood/FitsImage.h"
#include "Likelihood/MapCubeFunction.h"
#include "Likelihood/SkyDirArg.h"

namespace {
   bool reverse_cmp(double x, double y) {
      return x > y;
   }
   double interpolatePowerLaw(double x, double x1, double x2,
                              double y1, double y2) {
      if (x1 <= 0 || x2 <= 0 || y1 <= 0 || y2 <= 0) {
         std::ostringstream message;
         message << "MapCubeFunction::interpolatePowerLaw:\n"
                 << "abscissa or ordinate values found that are <= 0: "
                 << "x1 = " << x1 << ", "
                 << "x2 = " << x2 << ", "
                 << "y1 = " << y1 << ", "
                 << "y2 = " << y2 << std::endl;
         throw std::runtime_error(message.str());
      }
      double gamma = std::log(y2/y1)/std::log(x2/x1);
      double n0 = y1/std::pow(x1, gamma);
      return n0*std::pow(x, gamma);
   }
}

namespace Likelihood {

double MapCubeFunction::value(optimizers::Arg & x) const {
   if (m_image.size() == 0) {
      return 0;
   }

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

   int i = findIndex(m_lon, lonValue) - 1;
   int j = findIndex(m_lat, latValue) - 1;
   unsigned int k = findIndex(m_energies, energy);
   k = std::min(k, m_energies.size() - 2);

   int nlon = m_lon.size() - 1;
   int nlat = m_lat.size() - 1;

   int indx = k*nlon*nlat + j*nlon + i;
   double y1 = m_image.at(indx);
   indx = (k+1)*nlon*nlat + j*nlon + i;
   double y2 = m_image.at(indx);
   double value = ::interpolatePowerLaw(energy, m_energies.at(k),
                                        m_energies.at(k+1), y1, y2);
   return value*getParam("Normalization").getTrueValue();
}

void MapCubeFunction::init() {
   setMaxNumParams(1);
// Leave this parameter fixed, modifying the overall normalization
// via a ConstantValue function as the spectral component.
   addParam("Normalization", 1, false);

   m_funcType = Addend;
   m_argType = "";
   m_genericName = "MapCubeFunction";
}

void MapCubeFunction::readFitsFile(const std::string & fitsFile) {
   m_fitsFile = fitsFile;
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

   fitsImage.getPixelBounds(0, m_lon);
   fitsImage.getPixelBounds(1, m_lat);
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
}

int MapCubeFunction::findIndex(std::vector<double> xx, double x) {
   if (xx.front() < xx.back()) {
      return std::upper_bound(xx.begin(), xx.end(), x) - xx.begin();
   }
   return std::upper_bound(xx.begin(), xx.end(), x, ::reverse_cmp) 
      - xx.begin();
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

double MapCubeFunction::mapIntegral() const {
   FitsImage fitsImage(m_fitsFile);
   std::vector<double> solidAngles;

   fitsImage.getSolidAngles(solidAngles);

   unsigned int nlon(m_lon.size()-1);
   unsigned int nlat(m_lat.size()-1);

   double map_integral(0);
   for (unsigned int j = 0; j < nlat; j++) {
      for (unsigned int i = 0; i < nlon; i++) {
         for (unsigned int k = 1; k < m_energies.size(); k++) {
            unsigned int indx = k*nlon*nlat + j*nlon + i;
            map_integral += solidAngles.at(j*nlon + i)*
               powerLawIntegral(m_energies.at(k-1), m_energies.at(k),
                                m_image.at(indx-nlon*nlat), m_image.at(indx));
         }
      }
   }
   return map_integral;
}

double MapCubeFunction::powerLawIntegral(double x1, double x2,
                                         double y1, double y2) const {
   double gamma = std::log(y2/y1)/std::log(x2/x1);
   double n0 = y1/std::pow(x1, gamma);
   double integral;
   if (gamma != 1.) {
      double gp1 = gamma + 1.;
      integral = n0/gp1*(std::pow(x2, gp1) - std::pow(x1, gp1));
   } else {
      integral = n0*std::log(x2/x1);
   }
   return integral;
}

}
