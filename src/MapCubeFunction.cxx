/**
 * @file MapCubeFunction.cxx
 * @brief Encapsulation of 3D FITS image of a diffuse source with 
 * position-dependent spectral variation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/MapCubeFunction.cxx,v 1.34 2011/02/09 23:11:50 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "facilities/Util.h"

#include "tip/Header.h"
#include "tip/Image.h"
#include "tip/IFileSvc.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/FitsImage.h"
#include "st_facilities/Util.h"

#include "Likelihood/ExposureMap.h"
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

   double my_round(double x) {
      int xint = static_cast<int>(x);
      if (x - xint >= 0.5) {
         return xint + 1.;
      }
      return xint;
   }
}

namespace Likelihood {

MapCubeFunction::MapCubeFunction() 
   : optimizers::Function(), MapBase(), 
     m_proj(0), m_nlon(0), m_nlat(0), m_isPeriodic(false) {
   init();
}

MapCubeFunction::MapCubeFunction(const std::string & fitsFile) 
   : optimizers::Function(), MapBase(fitsFile), m_proj(0), m_nlon(0), 
     m_nlat(0), m_isPeriodic(false)  {
   init();
   readFitsFile(fitsFile);
}

MapCubeFunction::MapCubeFunction(const MapCubeFunction & rhs)
   : optimizers::Function(rhs), MapBase(rhs),
     m_nlon(rhs.m_nlon), m_nlat(rhs.m_nlat), m_isPeriodic(rhs.m_isPeriodic) {
// astro::SkyProj copy constructor is not implemented properly so we
// must share this pointer, ensure it is not deleted in the destructor,
// and live with the resulting memory leak when this object is deleted.
//   m_proj = new astro::SkyProj(*(rhs.m_proj));
   m_proj = rhs.m_proj;
   m_energies = rhs.m_energies;
   m_image = rhs.m_image;
   m_mapIntegrals = rhs.m_mapIntegrals;
}

MapCubeFunction & MapCubeFunction::operator=(const MapCubeFunction &rhs) {
   if (this != &rhs) {
// astro::SkyProj copy constructor is not implemented properly so we
// must share this pointer, ensure it is not deleted in the destructor,
// and live with the resulting memory leak when this object is deleted.
//       delete m_proj;
//       m_proj = new astro::SkyProj(*(rhs.m_proj));
      m_proj = rhs.m_proj;
      optimizers::Function::operator=(rhs);
      MapBase::operator=(rhs);
      m_nlon = rhs.m_nlon;
      m_nlat = rhs.m_nlat;
      m_energies = rhs.m_energies;
      m_image = rhs.m_image;
      m_mapIntegrals = rhs.m_mapIntegrals;
   }
   return *this;
}

MapCubeFunction::~MapCubeFunction() {
// astro::SkyProj copy constructor is not implemented properly so we
// must ensure this pointer is not deleted here and live with the
// resulting memory leak when this object is deleted.
//   delete m_proj;
}

double MapCubeFunction::value(optimizers::Arg & xarg) const {
   if (m_image.size() == 0) {
      return 0;
   }

   SkyDirArg & dir = dynamic_cast<SkyDirArg &>(xarg);
   double energy = dir.energy();

   size_t k = findIndex(m_energies, energy) - 1;
   k = std::min(k, m_energies.size() - 2);

   std::pair<double, double> pixel = dir().project(*m_proj);

   double x(pixel.first);
   double y(pixel.second);

   if (m_isPeriodic) {
      x = std::fmod(x, m_nlon);
   }

   if ((!m_isPeriodic && (x < 0.5 || x > m_nlon + 0.5)) ||
       y < 0.5 || y > m_nlat + 0.5) {
      // Sky location is outside of map, so do not extrapolate and return 0.
      return 0;
   }

// This code tries to do a bilinear interpolation on the pixel values.
   int ix(static_cast<int>(x));
   int iy(static_cast<int>(y));

// Points within half a pixel of the edges of the map need to be
// extrapolated in the context of the bilinear scheme, even though
// they are formally inside the map.
   if (!m_isPeriodic) {
      if (ix < 1) {
         ix = 1;
      }
      if (ix >= m_nlon) {
         ix = m_nlon - 1;
      }
   }
   if (iy < 1) {
      iy = 1;
   }
   if (iy >= m_nlat) {
      iy = m_nlat - 1;
   }
   
// Perform bilinear interpolation for the image planes that bracket
// the desired energy.
   double y1 = bilinear_interpolation(k, ix, iy, x, y);
   double y2 = bilinear_interpolation(k + 1, ix, iy, x, y);
   if (y1 == 0 || y2 == 0) {
      return 0;
   }
   double value = ::interpolatePowerLaw(energy, m_energies.at(k),
                                        m_energies.at(k+1), y1, y2);
   return value*getParam("Normalization").getTrueValue();
}

double MapCubeFunction::
bilinear_interpolation(size_t k, int ix, int iy, double x, double y) const {
   double tt(x - ix);
   double uu(y - iy);

   double y1, y4;

   if (m_isPeriodic && ix == 0) {
      y1 = m_image.at((k*m_nlat + iy - 1)*m_nlon + m_nlon - 1);
      y4 = m_image.at((k*m_nlat + iy)*m_nlon + m_nlon - 1);
   } else {
      y1 = m_image.at((k*m_nlat + iy - 1)*m_nlon + ix - 1);
      y4 = m_image.at((k*m_nlat + iy)*m_nlon + ix - 1);
   }
   double y2(m_image.at((k*m_nlat + iy - 1)*m_nlon + ix));
   double y3(m_image.at((k*m_nlat + iy)*m_nlon + ix));
   
   double value((1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
                + tt*uu*y3 + (1. - tt)*uu*y4);
   return value;
}                                               

void MapCubeFunction::init() {
   setMaxNumParams(1);
// Leave this parameter fixed, modifying the overall normalization
// via a ConstantValue function as the spectral component.
   addParam("Normalization", 1, false);

   m_funcType = Addend;
   m_argType = "";
   m_genericName = "MapCubeFunction";
   m_normParName = "Normalization";
}

void MapCubeFunction::readFitsFile(const std::string & fitsFile, 
                                   const std::string & extension,
                                   bool loadMap) {
   MapBase::readFitsFile(fitsFile, extension, loadMap);

   std::string expandedFileName(m_fitsFile);
   facilities::Util::expandEnvVar(&expandedFileName);
   ExposureMap::readEnergyExtension(expandedFileName, m_energies);
}

void MapCubeFunction::readFitsFile() {
   MapBase::readFitsFile();

   std::string expandedFileName(m_fitsFile);
   facilities::Util::expandEnvVar(&expandedFileName);

   m_proj = new astro::SkyProj(expandedFileName);

   st_facilities::FitsImage fitsImage(expandedFileName);
   m_image = fitsImage.imageData();
   
   std::vector<int> naxes;
   fitsImage.getAxisDims(naxes);
   m_nlon = naxes.at(0);
   m_nlat = naxes.at(1);
   
   double cdelt1;
   const tip::Image * image =
      tip::IFileSvc::instance().readImage(expandedFileName, "");
   const tip::Header & header = image->getHeader();
   header["CDELT1"].get(cdelt1);
   if (::my_round(m_nlon*cdelt1) == 360.) {
      m_isPeriodic = true;
   }
   delete image;
}

void MapCubeFunction::deleteMap() {
   MapBase::deleteMap();
   m_image.clear();
}

int MapCubeFunction::
findIndex(const std::vector<double> & xx, double x) const {
   if (xx.front() < xx.back()) {
      return std::upper_bound(xx.begin(), xx.end(), x) - xx.begin();
   }
   return std::upper_bound(xx.begin(), xx.end(), x, ::reverse_cmp) 
      - xx.begin();
}

double MapCubeFunction::mapIntegral() const {
   const std::vector< std::vector<double> > & solidAngles 
      = wcsmap().solidAngles();
   double map_integral(0);
   for (int j = 0; j < m_nlat; j++) {
      for (int i = 0; i < m_nlon; i++) {
         for (unsigned int k = 1; k < m_energies.size(); k++) {
            unsigned int indx = (k*m_nlat + j)*m_nlon + i;
            map_integral += solidAngles.at(i).at(j)*
               powerLawIntegral(m_energies.at(k-1), m_energies.at(k),
                                m_image.at(indx-m_nlon*m_nlat),
                                m_image.at(indx));
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

double MapCubeFunction::mapIntegral(double energy) const {
   if (energy < m_energies.front() || energy > m_energies.back()) {
      return 0;
   }
   size_t k = std::upper_bound(m_energies.begin(), m_energies.end(),
                               energy) - m_energies.begin();
   if (m_mapIntegrals.empty()) {
      const_cast<MapCubeFunction *>(this)->computeMapIntegrals();
   }
   if (energy == m_energies.at(k)) {
      return m_mapIntegrals.at(k);
   }
   double value = (m_mapIntegrals.at(k-1)*
                   std::exp((std::log(energy/m_energies.at(k-1)))
                            /(std::log(m_energies.at(k)/m_energies.at(k-1)))
                            *std::log(m_mapIntegrals.at(k)/m_mapIntegrals.at(k-1))));
   return value;
}

void MapCubeFunction::computeMapIntegrals() {
   const std::vector< std::vector<double> > & solidAngles
      = wcsmap().solidAngles();
   m_mapIntegrals.clear();
   for (size_t k(0); k < m_energies.size(); k++) {
      m_mapIntegrals.push_back(0);
      for (int j(0); j < m_nlat; j++) {
         for (int i(0); i < m_nlon; i++) {
            size_t indx((k*m_nlat + j)*m_nlon + i);
            m_mapIntegrals.back() += solidAngles.at(i).at(j)*m_image.at(indx);
         }
      }
   }
}

}
