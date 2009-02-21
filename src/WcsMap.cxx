/**
 * @file WcsMap.cxx
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/WcsMap.cxx,v 1.33 2009/02/18 18:13:38 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/Convolve.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/WcsMap.h"

namespace {
   class Image : public std::vector< std::vector<double> > {
   public:
      Image() {}
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

   double my_round(double x) {
      int xint = static_cast<int>(x);
      if (x - xint >= 0.5) {
         return xint + 1.;
      }
      return xint;
   }
} // unnamed namespace

namespace Likelihood {

WcsMap::WcsMap() : m_proj(0), m_interpolate(false), m_mapIntegral(0) {}


WcsMap::WcsMap(const std::string & filename,
               const std::string & extension,
               bool interpolate) 
   : m_proj(0), m_interpolate(interpolate), m_isPeriodic(false),
     m_mapIntegral(0) {
   
   m_proj = new astro::SkyProj(filename, extension);
   
   const tip::Image * image = 
      tip::IFileSvc::instance().readImage(filename, extension);
   
   std::vector<float> my_image;
   image->get(my_image);
   
   const tip::Header & header = image->getHeader();
   
   header["NAXIS1"].get(m_naxis1);
   header["NAXIS2"].get(m_naxis2);

   double ix, iy;
   header["CRPIX1"].get(ix);
   header["CRPIX2"].get(iy);

   double cdelt1;
   header["CDELT1"].get(cdelt1);
   if (::my_round(m_naxis1*cdelt1) == 360.) {
      m_isPeriodic = true;
   }

   if (m_proj->isGalactic()) {
      m_coordSys = astro::SkyDir::GALACTIC;
   } else {
      m_coordSys = astro::SkyDir::EQUATORIAL;
   }

   std::pair<double, double> coord = m_proj->pix2sph(ix, iy);
   m_refDir = astro::SkyDir(coord.first, coord.second, m_coordSys);

   delete image;

   m_image.reserve(m_naxis2);
   std::vector<float>::const_iterator pixelValue = my_image.begin();
   for (int j = 0; j < m_naxis2; j++) {
      std::vector<double> row(m_naxis1);
      for (int i = 0; i < m_naxis1 && pixelValue != my_image.end();
           i++, ++pixelValue) {
         row.at(i) = *pixelValue;
      }
      m_image.push_back(row);
   }

   computeMapIntegral();
}

WcsMap::WcsMap(const DiffuseSource & diffuseSource,
               double ra, double dec, double pix_size, int npts,
               double energy, const std::string & proj_name, bool use_lb,
               bool interpolate) 
   : m_refDir(ra, dec), m_interpolate(interpolate), m_mapIntegral() {
   if (use_lb) { // convert to l, b
      ra = m_refDir.l();
      dec = m_refDir.b();
   }
   double crpix[] = {npts/2., npts/2.};
   double crval[] = {ra, dec};
   double cdelt[] = {-pix_size, pix_size};

   m_proj = new astro::SkyProj(proj_name, crpix, crval, cdelt, 0, use_lb);

   if (use_lb) {
      m_coordSys = astro::SkyDir::GALACTIC;
   } else {
      m_coordSys = astro::SkyDir::EQUATORIAL;
   }

   m_image.reserve(npts);
   double ix, iy;
   for (int j = 0; j < npts; j++) {
      iy = j + 1.;
      std::vector<double> row(npts, 0);
      for (int i = 0; i < npts; i++) {
         ix = i + 1.;
         if (m_proj->testpix2sph(ix, iy) == 0) {
            std::pair<double, double> coord = m_proj->pix2sph(ix, iy);
            astro::SkyDir dir(coord.first, coord.second, m_coordSys);
            SkyDirArg my_dir(dir, energy);
            row.at(i) = diffuseSource.spatialDist(my_dir);
         } else {
            row.at(i) = 0;
         }
      }
      m_image.push_back(row);
   }
   m_naxis1 = npts;
   m_naxis2 = npts;

   computeMapIntegral();
}

WcsMap::~WcsMap() {
// // astro::SkyProj copy constructor is not implemented properly so we
// // must ensure this pointer is not deleted here and live with the
// // resulting memory leak when this object is deleted.
//   delete m_proj;
}

WcsMap::WcsMap(const WcsMap & rhs) 
   : m_refDir(rhs.m_refDir), m_image(rhs.m_image), 
     m_naxis1(rhs.m_naxis1), m_naxis2(rhs.m_naxis2), 
     m_interpolate(rhs.m_interpolate),
     m_isPeriodic(rhs.m_isPeriodic),
     m_coordSys(rhs.m_coordSys),
     m_mapIntegral(rhs.m_mapIntegral) {
// astro::SkyProj copy constructor is not implemented properly so we
// must share this pointer, ensure it is not deleted in the destructor,
// and live with the resulting memory leak when this object is deleted.
//   m_proj = new astro::SkyProj(*(rhs.m_proj));
   m_proj = rhs.m_proj;
}

WcsMap & WcsMap::operator=(const WcsMap & rhs) {
   if (this != &rhs) {
// astro::SkyProj copy constructor is not implemented properly so we
// must share this pointer, ensure it is not deleted in the destructor,
// and live with the resulting memory leak when this object is deleted.
//      delete m_proj;
//      m_proj = new astro::SkyProj(*(rhs.m_proj));
      m_proj = rhs.m_proj;
      m_refDir = rhs.m_refDir;
      m_image = rhs.m_image;
      m_naxis1 = rhs.m_naxis1;
      m_naxis2 = rhs.m_naxis2;
      m_interpolate = rhs.m_interpolate;
      m_mapIntegral = rhs.m_mapIntegral;
   }
   return *this;
}

double WcsMap::operator()(const astro::SkyDir & dir) const {
// NB: wcslib starts indexing pixels with 1, not 0.
   std::pair<double, double> pixel = dir.project(*m_proj);

   double x(pixel.first);
   double y(pixel.second);

   if (m_isPeriodic) {
      x = std::fmod(x, m_naxis1);
   }

   if (!m_interpolate) {
      return pixelValue(x, y);
   }
// This code tries to do a bilinear interpolation on the pixel values.

   int ix(static_cast<int>(x));
   int iy(static_cast<int>(y));

   if ((!m_isPeriodic && (ix < 1 || ix >= m_naxis1))
       || iy < 1 || iy >= m_naxis2) {
      return 0;
   }
   
   double tt(x - ix);
   double uu(y - iy);

   double y1, y4;

   if (m_isPeriodic && ix == 0) {
      y1 = m_image.at(iy-1).back();
      y4 = m_image.at(iy).back();
   } else {
      y1 = m_image.at(iy-1).at(ix-1);
      y4 = m_image.at(iy).at(ix-1);
   }
   double y2(m_image.at(iy-1).at(ix));
   double y3(m_image.at(iy).at(ix));
   
   double value((1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
                + tt*uu*y3 + (1. - tt)*uu*y4);

//    if (value < 0) {
//       value = 0;
//    }
   
   return value;
}

WcsMap WcsMap::convolve(double energy, const MeanPsf & psf,
                        const BinnedExposure & exposure,
                        bool performConvolution) const {
// Compute unconvolved counts map by multiplying intensity image by exposure.
   ::Image counts;
   counts.resize(m_naxis2);

   for (int j = 0; j < m_naxis2; j++) {
      counts.at(j).resize(m_naxis1, 0);
      for (int i = 0; i < m_naxis1; i++) {
         if (m_proj->testpix2sph(i+1, j+1) == 0) {
            std::pair<double, double> coord = m_proj->pix2sph(i+1, j+1);
            astro::SkyDir dir(coord.first, coord.second, m_coordSys);
            counts.at(j).at(i) = 
               m_image.at(j).at(i)*exposure(energy, dir.ra(), dir.dec());
         }
      }
   }

// Fill a square array with an image of the Psf at the same binning
// resolution as the source image.  Use the smaller of the image map
// dimensions to determine the psf image size.
// @todo handle cases where m_naxis1 or m_naxis2 are odd.
   int npix, imin, imax, jmin, jmax;
   if (m_naxis1 <= m_naxis2) {
      npix = m_naxis1;
      imin = 0;
      imax = npix;
      jmin = (m_naxis2 - m_naxis1)/2;
      jmax = (m_naxis2 + m_naxis1)/2;
   } else {
      npix = m_naxis2;
      imin = (m_naxis1 - m_naxis2)/2;
      imax = (m_naxis1 + m_naxis2)/2;
      jmin = 0;
      jmax = npix;
   }

   ::Image psf_image;
   psf_image.resize(npix);

   for (int j = jmin; j < jmax; j++) {
      psf_image.at(j).resize(npix, 0);
      for (int i = imin; i < imax; i++) {
         if (m_proj->testpix2sph(i+1, j+1) == 0) {
            std::pair<double, double> coord = m_proj->pix2sph(i+1, j+1);
            astro::SkyDir dir(coord.first, coord.second, m_coordSys);
            double theta = m_refDir.difference(dir)*180./M_PI;
            psf_image.at(j-jmin).at(i-imin) = psf(energy, theta, 0);
         }
      }
   }

   psf_image.normalize();

   WcsMap my_image(*this);
   if (performConvolution) {
      my_image.m_image = Convolve::convolve2d(counts, psf_image);
   } else {
      my_image.m_image = counts;
   }

   return my_image;
}

double WcsMap::solidAngle(const astro::SkyProj & proj, 
                          double ilon, double ilat) {
   std::pair<double, double> center(proj.pix2sph(ilon, ilat));
   std::pair<double, double> left(proj.pix2sph(ilon - 0.5, ilat));
   std::pair<double, double> right(proj.pix2sph(ilon + 0.5, ilat));
   std::pair<double, double> bottom(proj.pix2sph(ilon, ilat - 0.5));
   std::pair<double, double> top(proj.pix2sph(ilon, ilat + 0.5));

   double cos_lat = std::fabs(std::cos(center.second*M_PI/180.));
   double delta_lon = std::fabs(std::fmod((right.first - left.first + 360.),
                                          360.))*M_PI/180.;
   double delta_lat = (top.second - bottom.second)*M_PI/180.;

   double dOmega = std::fabs(delta_lon*delta_lat*cos_lat);

   return dOmega;
}

double WcsMap::solidAngle(double ilon, double ilat) const {
   return solidAngle(*m_proj, ilon, ilat);
}

double WcsMap::pixelValue(double ilon, double ilat) const {
// Find the pixel in which the sky location lives and returns the value.
   int ix(static_cast<int>(::my_round(ilon)) - 1);
   int iy(static_cast<int>(::my_round(ilat)) - 1);

   if ((!m_isPeriodic && (ix < 0 || ix >= m_naxis1)) 
       || iy < 0 || iy >= m_naxis2) {
      return 0;
   }
   if (ix < 0 && ix >= -1) {
      ix = 0;
   }
   return m_image.at(iy).at(ix);
}

astro::SkyDir WcsMap::skyDir(double ilon, double ilat) const {
   std::pair<double, double> coords(m_proj->pix2sph(ilon, ilat));
   return astro::SkyDir(coords.first, coords.second, m_coordSys);
}

bool WcsMap::insideMap(const astro::SkyDir & dir) const {
   std::pair<double, double> pixel = dir.project(*m_proj);

   double x(pixel.first);
   double y(pixel.second);

   int ix(static_cast<int>(x));
   int iy(static_cast<int>(y));

   if ((!m_isPeriodic && (ix < 1 || ix >= m_naxis1))
       || iy < 1 || iy >= m_naxis2) {
      return false;
   }
   return true;
}

std::pair<astro::SkyDir, astro::SkyDir> 
WcsMap::minMaxDistPixels(const astro::SkyDir & dir) const {
   astro::SkyDir closest(skyDir(1, 1));
   double min_dist = dir.difference(closest);
   astro::SkyDir farthest(skyDir(1, 1));
   double max_dist = dir.difference(closest);
   int i(2);
   int j(1);
   for ( ; i < m_naxis1 + 1; i++) { // i = 2, m_naxis1; j = 1
      astro::SkyDir current(skyDir(i, j));
      double dist(dir.difference(current));
      if (dist < min_dist) {
         min_dist = dist;
         closest = current;
      }
      if (dist > max_dist) {
         max_dist = dist;
         farthest = current;
      }
   }
   for (j = 2 ; j < m_naxis2 + 1; j++) { // i = m_naxis1; j = 2, m_naxis2
      astro::SkyDir current(skyDir(i, j));
      double dist(dir.difference(current));
      if (dist < min_dist) {
         min_dist = dist;
         closest = current;
      }
      if (dist > max_dist) {
         max_dist = dist;
         farthest = current;
      }
   }
   for (i = 1; i < m_naxis1; i++) { // i = 1, m_naxis1-1; j = m_naxis2
      astro::SkyDir current(skyDir(i, j));
      double dist(dir.difference(current));
      if (dist < min_dist) {
         min_dist = dist;
         closest = current;
      }
      if (dist > max_dist) {
         max_dist = dist;
         farthest = current;
      }
   }
   i = 1;
   for (j = 2; j < m_naxis2; j++) { // i = 1; j = 2, m_naxis2-1
      astro::SkyDir current(skyDir(i, j));
      double dist(dir.difference(current));
      if (dist < min_dist) {
         min_dist = dist;
         closest = current;
      }
      if (dist > max_dist) {
         max_dist = dist;
         farthest = current;
      }
   }
   return std::make_pair(closest, farthest);
}

void WcsMap::getCorners(std::vector<astro::SkyDir> & corners) const {
   corners.clear();
   corners.push_back(skyDir(1, 1));
   corners.push_back(skyDir(1, m_naxis2));
   corners.push_back(skyDir(m_naxis1, m_naxis2));
   corners.push_back(skyDir(m_naxis1, 1));
}

double WcsMap::mapIntegral() const {
   return m_mapIntegral;
}

void WcsMap::computeMapIntegral() {
   m_mapIntegral = 0;
   for (int i(0); i < m_naxis1; i++) {
      for (int j(0); j < m_naxis2; j++) {
         m_mapIntegral += solidAngle(i, j)*m_image.at(j).at(i);
      }
   }
}

} // namespace Likelihood
