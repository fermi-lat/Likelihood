/**
 * @file MapBase.cxx
 * @brief Base class to provide the functionality needed to find the limits
 * for diffuse response integrals
 *
 * @author J. Chiang
 * 
 * $Header$
 */

#include "Likelihood/MapBase.h"

//namespace Likelihood {

MapBase::MapBase() : m_wcsmap(0), m_fitsFile(""), m_extension("") {}

MapBase::MapBase(const std::string & fitsFile, const std::string & extension) 
   : m_wcsmap(new WcsMap(fitsFile, extension)),
     m_fitsFile(fitsFile), m_extension(extension) {}

MapBase::~MapBase() {
   delete m_wcsmap;
}

MapBase::MapBase(const MapBase & other) 
   : m_wcsmap(new WcsMap(*(other.m_wcsmap))), m_fitsFile(other.m_fitsFile),
     m_extension(other.m_extension) {}

MapBase & MapBase::operator=(const MapBase & rhs) {
   if (this != &rhs) {
      delete m_wcsmap;
      if (rhs.m_wcsmap) {
         m_wcsmap = new WcsMap(*(rhs.m_wcsmap));
      }
      m_fitsFile = rhs.m_fitsFile;
      m_extension = rhs.m_extension;
   }
   return *this;
}

bool MapBase::insideMap(const astro::SkyDir & dir) const {
   return m_wcsmap->insideMap(dir);
}

void MapBase::getDiffRespLimits(double & mumin, double & mumax,
                                double & phimin, double & phimax) const {

}

void MapBase::getMinMaxDistPixels(const astro::SkyDir &,
                                  astro::SkyDir & closestPixel, 
                                  astro::SkyDir & farthestPixel) const {

}

void MapBase::getCorners(std::vector<astro::SkyDir & corners) const {

}

} // namespace Likelihood

