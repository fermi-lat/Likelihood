/**
 * @file WcsMap.cxx
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/WcsMap.cxx,v 1.17 2006/06/13 03:30:03 jchiang Exp $
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
} // unnamed namespace

namespace Likelihood {

WcsMap::WcsMap(const std::string & filename,
               const std::string & extension) : m_proj(0) {

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

   astro::SkyDir::CoordSystem coordSys;
   if (m_proj->isGalactic()) {
      coordSys = astro::SkyDir::GALACTIC;
   } else {
      coordSys = astro::SkyDir::EQUATORIAL;
   }

   std::pair<double, double> coord = m_proj->pix2sph(ix, iy);
   m_refDir = astro::SkyDir(coord.first, coord.second, coordSys);

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
}

WcsMap::WcsMap(const DiffuseSource & diffuseSource,
               double ra, double dec, double radius, int npts,
               double energy, const std::string & proj_name, bool use_lb) 
   : m_refDir(ra, dec) {
   if (use_lb) { // convert to l, b
      ra = m_refDir.l();
      dec = m_refDir.b();
   }
   double crpix[] = {npts/2. + 0.5, npts/2. + 0.5};
   double crval[] = {ra, dec};
   double cdelt[] = {2.*radius/npts, 2.*radius/npts};

   m_proj = new astro::SkyProj(proj_name, crpix, crval, cdelt, 0, use_lb);

   astro::SkyDir::CoordSystem coordSys;
   if (use_lb) {
      coordSys = astro::SkyDir::GALACTIC;
   } else {
      coordSys = astro::SkyDir::EQUATORIAL;
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
            astro::SkyDir dir(coord.first, coord.second, coordSys);
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
}

WcsMap::~WcsMap() {
// // astro::SkyProj copy constructor is not implemented properly so we
// // must ensure this pointer is not deleted here and live with the
// // resulting memory leak when this object is deleted.
//   delete m_proj;
}

WcsMap::WcsMap(const WcsMap & rhs) 
   : m_refDir(rhs.m_refDir), m_image(rhs.m_image), 
     m_naxis1(rhs.m_naxis1), m_naxis2(rhs.m_naxis2) {
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
   }
   return *this;
}

double WcsMap::operator()(const astro::SkyDir & dir) const {
// NB: wcslib starts indexing pixels with 1, not 0.
   std::pair<double, double> pixel = dir.project(*m_proj);

   double x(pixel.first);
   double y(pixel.second);

   int ix(static_cast<int>(std::floor(x - 0.5)));
   int iy(static_cast<int>(std::floor(y - 0.5)));

   if (ix < 0 || ix >= m_naxis1 || iy < 0 || iy >= m_naxis2) {
      return 0;
   }

   return m_image.at(iy).at(ix);

//    ix = std::min(std::max(1, ix), m_naxis1 - 1);
//    iy = std::min(std::max(1, iy), m_naxis2 - 1);

//    double uu(x - ix);
//    double tt(y - iy);

//    double y1(m_image.at(iy-1).at(ix-1));
//    double y2(m_image.at(iy).at(ix-1));
//    double y3(m_image.at(iy).at(ix));
//    double y4(m_image.at(iy-1).at(ix));

//    double value((1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
//                 + tt*uu*y3 + (1. - tt)*uu*y4);

//    if (value < 0) {
//       throw std::runtime_error("WcsMap::operator(): value < 0");
//    }

//    return value;
}

WcsMap WcsMap::convolve(double energy, const MeanPsf & psf,
                        const BinnedExposure & exposure,
                        bool performConvolution) const {
   ::Image counts;
   counts.resize(m_naxis2);
   ::Image psf_image;
   psf_image.resize(m_naxis2);

   astro::SkyDir::CoordSystem coordSys;
   if (m_proj->isGalactic()) {
      coordSys = astro::SkyDir::GALACTIC;
   } else {
      coordSys = astro::SkyDir::EQUATORIAL;
   }

   for (int j = 0; j < m_naxis2; j++) {
      counts.at(j).resize(m_naxis1, 0);
      psf_image.at(j).resize(m_naxis1);
      for (int i = 0; i < m_naxis1; i++) {
         if (m_proj->testpix2sph(i+1, j+1) == 0) {
            std::pair<double, double> coord = m_proj->pix2sph(i+1, j+1);
            astro::SkyDir dir(coord.first, coord.second, coordSys);
            counts.at(j).at(i) = 
               m_image.at(j).at(i)*exposure(energy, dir.ra(), dir.dec());
            double theta = m_refDir.difference(dir)*180./M_PI;
            psf_image.at(j).at(i) = psf(energy, theta, 0);
         } else {
            psf_image.at(j).at(i) = 0;
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

} // namespace Likelihood
