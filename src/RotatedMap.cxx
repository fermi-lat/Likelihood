/**
 * @file RotatedMap.cxx
 * @brief Map-like representation of a diffuse source that is
 * rotated to an equinox-centered coordinate system.
 * @author J. Chiang
 *
 * $Header$
 */

#include "st_facilities/Util.h"

#include "Likelihood/Convolve.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/RotatedMap.h"

namespace Likelihood {

RotatedMap::RotatedMap(const DiffuseSource & diffuseSource,
                       double ra, double dec, double radius, int npts) {
   m_rot = EquinoxRotation(ra, dec);
   linearArray(-radius, radius, npts, m_lons);
   linearArray(-radius, radius, npts, m_lats);
   m_image.reserve(npts);
   astro::SkyDir trueDir(0, 0);
   for (unsigned int i = 0; i < m_lons.size(); i++) {
      std::vector<double> row;
      for (unsigned int j = 0; j < m_lats.size(); j++) {
         astro::SkyDir localDir(m_lons.at(i), m_lats.at(j));
         m_rot.do_rotation(localDir, trueDir);
         row.push_back(diffuseSource.spatialDist(trueDir));
      }
      m_image.push_back(row);
   }
}

double RotatedMap::operator()(const astro::SkyDir & dir) const {
   astro::SkyDir mapDir(0, 0);
   m_rot.do_rotation(dir, mapDir, true);
   double lon(mapDir.ra());
   if (lon > 180.) {
      lon -= 360.;
   }
   double lat(mapDir.dec());
   std::vector<double>::const_iterator lon_it
      = std::upper_bound(m_lons.begin(), m_lons.end(), lon);
   if (lon_it == m_lons.begin() || lon_it == m_lons.end()) {
      return 0;
   }
   std::vector<double>::const_iterator lat_it
      = std::upper_bound(m_lats.begin(), m_lats.end(), lat);
   if (lat_it == m_lats.begin() || lat_it == m_lats.end()) {
      return 0;
   }
   return st_facilities::Util::bilinear(m_lons, lon, m_lats, lat, m_image);
}

RotatedMap RotatedMap::convolve(const MeanPsf & psf, double energy) const {
   std::vector< std::vector<double> > psf_image;
   psf.getImage(energy, 0, 0, m_lons, m_lats, psf_image);
   return RotatedMap(Convolve::convolve2d(m_image, psf_image),
                     m_lons, m_lats, m_rot);
}

void RotatedMap::getUnrotatedMap(std::vector< std::vector<double> > & map,
                                 const std::vector<double> & ras,
                                 const std::vector<double> & decs,
                                 astro::SkyDir::CoordSystem coordSys) const {
   map.clear();
   map.reserve(ras.size());
   for (unsigned int i = 0; i < ras.size(); i++) {
      std::vector<double> row;
      for (unsigned int j = 0; j < decs.size(); j++) {
         astro::SkyDir dir(ras.at(i), decs.at(j), coordSys);
         row.push_back((*this)(dir));
      }
      map.push_back(row);
   }
}

void RotatedMap::linearArray(double xmin, double xmax, int npts,
                             std::vector<double> & x) const {
   double xstep = (xmax - xmin)/(npts - 1);
   x.clear();
   x.reserve(npts);
   for (unsigned int i = 0; i < x.size(); i++) {
      x.push_back(xstep*i + xmin);
   }
}

} // namespace Likelihood
