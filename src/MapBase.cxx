/**
 * @file MapBase.cxx
 * @brief Base class to provide the functionality needed to find the limits
 * for diffuse response integrals
 *
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/MapBase.cxx,v 1.7 2010/07/07 01:05:30 jchiang Exp $
 */

#include <cmath>

#include <stdexcept>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Util.h"

#include "Likelihood/DiffRespIntegrand.h"
#include "Likelihood/EquinoxRotation.h"
#include "Likelihood/MapBase.h"

namespace Likelihood {

MapBase::MapBase() : m_wcsmap(0), m_fitsFile(""), m_extension("") {}

MapBase::MapBase(const std::string & fitsFile, const std::string & extension) 
   : m_wcsmap(0), m_fitsFile(fitsFile), m_extension(extension) {
   readFitsFile();
}

MapBase::~MapBase() {
   delete m_wcsmap;
}

MapBase::MapBase(const MapBase & other) 
   : m_wcsmap(0), m_fitsFile(other.m_fitsFile),
     m_extension(other.m_extension) {
   if (other.m_wcsmap) {
      m_wcsmap = new WcsMap(*(other.m_wcsmap));
   }
}

MapBase & MapBase::operator=(const MapBase & rhs) {
   if (this != &rhs) {
      delete m_wcsmap;
      m_wcsmap = 0;
      if (rhs.m_wcsmap) {
         m_wcsmap = new WcsMap(*(rhs.m_wcsmap));
      }
      m_fitsFile = rhs.m_fitsFile;
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
   std::string expandedFileName(m_fitsFile);

   facilities::Util::expandEnvVar(&expandedFileName);

   if (!st_facilities::Util::fileExists(expandedFileName)) {
// The following to StreamFormatter is necessary since Xerces seems to
// corrupt the exception handling when this method is called from
// SourceFactory::readXml and the program simply aborts.
      st_stream::StreamFormatter formatter("MapBase", "readFitsFile", 2);
      formatter.err() << "File not found: " << expandedFileName << std::endl;
      throw std::runtime_error("File not found: " + expandedFileName);
   }

   delete m_wcsmap;
   m_wcsmap = new WcsMap(expandedFileName, m_extension, true);
}

WcsMap & MapBase::wcsmap() {
   if (m_wcsmap == 0) {
      readFitsFile();
   }
   return *m_wcsmap;
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
   WcsMap * tmp = m_wcsmap->rebin(factor, average);
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

} // namespace Likelihood
