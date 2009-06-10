/**
 * @file MapCubeFunction.cxx
 * @brief Encapsulation of 3D FITS image of a diffuse source with 
 * position-dependent spectral variation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MapCubeFunction.cxx,v 1.30 2009/06/10 17:57:10 jchiang Exp $
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

double MapCubeFunction::value(optimizers::Arg & x) const {
   if (m_image.size() == 0) {
      return 0;
   }

   SkyDirArg & dir = dynamic_cast<SkyDirArg &>(x);
   double energy = dir.energy();

   size_t k = findIndex(m_energies, energy) - 1;
   k = std::min(k, m_energies.size() - 2);

   std::pair<double, double> pixel = dir().project(*m_proj);

// NB: wcslib (through astro::SkyProj) starts indexing pixels with
// 1, not 0, so apply correction here to avoid off-by-one error.
   int i = static_cast<int>(::my_round(pixel.first)) - 1;
   int j = static_cast<int>(::my_round(pixel.second)) - 1;

   if (m_isPeriodic && i >= m_nlon) {
      i -= m_nlon;
   }

   if ((!m_isPeriodic && (i < 0 || i >= m_nlon)) 
       || j < 0 || j >= m_nlat) {
      return 0;
   }
   if (i < 0 && i >= -1) {
      i = 0;
   }

   int indx = (k*m_nlat + j)*m_nlon + i;
   try {
      double y1 = m_image.at(indx);
      double y2 = m_image.at(indx + m_nlon*m_nlat);

      double value = ::interpolatePowerLaw(energy, m_energies.at(k),
                                           m_energies.at(k+1), y1, y2);
      return value*getParam("Normalization").getTrueValue();
   } catch (std::exception &) {
      return 0;
   }
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
                                   const std::string & extension) {
   MapBase::readFitsFile(fitsFile, extension);

   std::string expandedFileName(fitsFile);
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

   ExposureMap::readEnergyExtension(expandedFileName, m_energies);
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
      = m_wcsmap->solidAngles();
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
      = m_wcsmap->solidAngles();
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
