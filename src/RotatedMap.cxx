/**
 * @file RotatedMap.cxx
 * @brief Map-like representation of a diffuse source that is
 * rotated to an equinox-centered coordinate system.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RotatedMap.cxx,v 1.3 2005/05/26 00:23:25 jchiang Exp $
 */

#include <algorithm>

#include "st_facilities/Util.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/Convolve.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/RotatedMap.h"
#include "Likelihood/SkyDirArg.h"

namespace {
   class Image : public std::vector< std::vector<double> > {
   public:
      Image() {}
      void multiply(const std::vector< std::vector<double> > & other) {
         if (this->size() != other.size() ||
             this->at(0).size() != other.at(0).size()) {
            throw std::runtime_error("Image::multiply: array sizes do "
                                     "not match.");
         }
         for (unsigned int i = 0; i < this->size(); i++) {
            if (this->at(i).size() != other.at(i).size()) {
               throw std::runtime_error("Image::multiply: inconsistent "
                                        "row lengths.");
            }
            for (unsigned int j = 0; j < this->at(i).size(); j++) {
               this->at(i).at(j) *= other.at(i).at(j);
            }
         }
      }
      void normalize() {
         double total(0);
         for (unsigned int i = 0; i < this->size(); i++) {
            for (unsigned int j = 0; j < this->at(i).size(); j++) {
               total += this->at(i).at(j);
            }
         }
         for (unsigned int i = 0; i < this->size(); i++) {
            for (unsigned int j = 0; j < this->at(i).size(); j++) {
               this->at(i).at(j) /= total;
            }
         }
      }
   };
} // unnamed namespace

namespace Likelihood {

RotatedMap::RotatedMap(const DiffuseSource & diffuseSource,
                       double ra, double dec, double radius, int npts,
                       double energy) {
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
         SkyDirArg my_dir(trueDir, energy);
         row.push_back(diffuseSource.spatialDist(my_dir));
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

RotatedMap RotatedMap::convolve(double energy, const MeanPsf & psf,
                                const BinnedExposure & exposure) const {
   ::Image counts;
   exposure.getRotatedImage(energy, m_lons, m_lats, m_rot, counts);
   counts.multiply(m_image);
   ::Image psf_image;
   psf.getImage(energy, 0, 0, m_lons, m_lats, psf_image);
   psf_image.normalize();
   return RotatedMap(Convolve::convolve2d(counts, psf_image),
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
   x.resize(npts);
   for (unsigned int i = 0; i < x.size(); i++) {
      x.at(i) = xstep*i + xmin;
   }
}

} // namespace Likelihood
