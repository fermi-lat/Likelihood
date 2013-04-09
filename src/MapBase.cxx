/**
 * @file MapBase.cxx
 * @brief Base class to provide the functionality needed to find the limits
 * for diffuse response integrals
 *
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/MapBase.cxx,v 1.14 2012/06/19 04:02:38 jchiang Exp $
 */

#include <cmath>

#include <stdexcept>
#include <sstream>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Util.h"

#include "Likelihood/DiffRespIntegrand.h"
#include "Likelihood/EquinoxRotation.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/WcsMapLibrary.h"

namespace Likelihood {

MapBase::MapBase() : m_wcsmap(0), m_fitsFile(""), 
                     m_expandedFileName(""), m_extension("") {}

MapBase::MapBase(const std::string & fitsFile, const std::string & extension) 
   : m_wcsmap(0), m_fitsFile(fitsFile), m_extension(extension) {
/// Comment out so that fits file is not read in by default.  Intention is
/// to have fits file read in only when it is first needed, i.e., when
/// wcsmap() is called from subclasses.
//   readFitsFile();
}

MapBase::~MapBase() {
   WcsMapLibrary::instance()->remove_observer(this);
}

MapBase::MapBase(const MapBase & other) 
   : m_wcsmap(other.m_wcsmap), 
     m_fitsFile(other.m_fitsFile),
     m_expandedFileName(other.m_expandedFileName),
     m_extension(other.m_extension) {
}

MapBase & MapBase::operator=(const MapBase & rhs) {
   if (this != &rhs) {
      m_wcsmap = rhs.m_wcsmap;
      m_fitsFile = rhs.m_fitsFile;
      m_expandedFileName = rhs.m_expandedFileName;
      m_extension = rhs.m_extension;
   }
   return *this;
}

void MapBase::readFitsFile(const std::string & fitsFile,
                           const std::string & extension,
                           bool loadMap) {
   m_fitsFile = fitsFile;
   m_extension = extension;
   if (loadMap) {
      readFitsFile();
   }
}

void MapBase::readFitsFile() {
   m_expandedFileName = m_fitsFile;

   facilities::Util::expandEnvVar(&m_expandedFileName);

// The following to StreamFormatter is necessary since Xerces seems to
// corrupt the exception handling when this method is called from
// SourceFactory::readXml and the program simply aborts.
   st_stream::StreamFormatter formatter("MapBase", "readFitsFile", 2);
   if (!st_facilities::Util::fileExists(m_expandedFileName)) {
      formatter.err() << "File not found: " << m_expandedFileName << std::endl;
      throw std::runtime_error("File not found: " + m_expandedFileName);
   }

   formatter.info(4) << "MapBase::readFitsFile: creating WcsMap2 object" 
                     << std::endl;
   m_wcsmap = WcsMapLibrary::instance()->wcsmap(m_expandedFileName,
                                                m_extension);
   WcsMapLibrary::instance()->add_observer(this);
}

WcsMap2 & MapBase::wcsmap() {
   if (m_wcsmap == 0) {
      readFitsFile();
   }
   return *m_wcsmap;
}

void MapBase::deleteMap() {
   st_stream::StreamFormatter formatter("MapBase", "deleteMap", 2);
   formatter.info(4) << "MapBased::deleteMap: " << m_expandedFileName
                     << std::endl;
   WcsMapLibrary::instance()->delete_map(m_expandedFileName, m_extension);
   m_wcsmap = 0;
}

void MapBase::update() {
   if (!WcsMapLibrary::instance()->has_map(m_expandedFileName, m_extension)) {
      m_wcsmap = 0;
   }
}

bool MapBase::insideMap(const astro::SkyDir & dir) const {
   return wcsmap().insideMap(dir);
}

void MapBase::getDiffRespLimits(const astro::SkyDir & dir, 
                                double & mumin, double & mumax,
                                double & phimin, double & phimax) const {
   EquinoxRotation eqRot(dir.ra(), dir.dec());
   mumin = -1;
   mumax = 1;
   phimin = 0;
   phimax = 2.*M_PI;
   astro::SkyDir closest, farthest;
   getMinMaxDistPixels(dir, closest, farthest);
   mumin = farthest().dot(dir());
   if (!insideMap(dir)) {
      mumax = closest().dot(dir());
      std::vector<astro::SkyDir> corners;
      getCorners(corners);
      double mu(corners.front()().dot(dir()));
// Search for minimum and maximum phi-values, initializing both to the value
// associated with the first corner
      phimin = DiffRespIntegrand::phiValue(mu, corners.front(), eqRot);
      phimax = phimin;
      for (size_t k(1); k < corners.size(); k++) {
         mu = corners.at(k)().dot(dir());
         double phitest = DiffRespIntegrand::phiValue(mu, corners.at(k), eqRot);
         if (phitest < phimin) {
            phimin = phitest;
         }
         if (phitest > phimax) {
            phimax = phitest;
         }
      }
      double dphi(phimax - phimin);
// The opening angle of any source map from a point outside of that
// map cannot subtend more than 180 deg (is this true for any
// projection?), so remap (by 360 deg) and reorder phi values.
      if (dphi > M_PI) {
         if (phimin < 0) {
            double tmp = phimax;
            phimax = 2*M_PI + phimin;
            phimin = tmp;
         } else if (phimax < 0) {
            double tmp = phimin;
            phimin = phimax - 2*M_PI;
            phimax = tmp;
         }
      }
   }
}

void MapBase::rebin(unsigned int factor, bool average) {
   if (!m_wcsmap) {
      // Do nothing for now, though will need an implementation for 
      // MapCubeFunctions.
      return;
   }
   WcsMap2 * tmp = m_wcsmap->rebin(factor, average);
   deleteMap();
   m_wcsmap = tmp;
}

void MapBase::getMinMaxDistPixels(const astro::SkyDir & dir,
                                  astro::SkyDir & closestPixel, 
                                  astro::SkyDir & farthestPixel) const {
   std::pair<astro::SkyDir, astro::SkyDir> pixels =
      wcsmap().minMaxDistPixels(dir);
   closestPixel = pixels.first;
   farthestPixel = pixels.second;
}

void MapBase::getCorners(std::vector<astro::SkyDir> & corners) const {
   wcsmap().getCorners(corners);
}

double MapBase::interpolatePowerLaw(double x, double x1, double x2,
                                    double y1, double y2) {
   if (y1 == 0 && y2 == 0) {
      return 0;
   }
   if (x1 <= 0 || x2 <= 0 || y1 <= 0 || y2 <= 0) {
      std::ostringstream message;
      message << "MapBase::interpolatePowerLaw:\n"
              << "abscissa or ordinate values found that are <= 0: "
              << "x1 = " << x1 << ", "
              << "x2 = " << x2 << ", "
              << "y1 = " << y1 << ", "
              << "y2 = " << y2 << std::endl;
      throw std::runtime_error(message.str());
   }
   double gamma = std::log(y2/y1)/std::log(x2/x1);
   // double n0 = y1/std::pow(x1, gamma);
   // return n0*std::pow(x, gamma);
   return y1*std::pow(x/x1, gamma);
}

} // namespace Likelihood
