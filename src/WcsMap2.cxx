/**
 * @file WcsMap2.cxx
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/WcsMap2.cxx,v 1.20 2015/10/17 17:19:17 mdwood Exp $
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/Convolve.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/SpatialFunction.h"
#include "Likelihood/WcsMap2.h"

namespace {
//   class Image : public std::vector< std::vector<double> > {
   class Image : public std::vector< std::vector<float> > {
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

WcsMap2::WcsMap2() : m_refDir(0, 0), m_proj(0), m_interpolate(false),
                     m_mapIntegral(0), m_enforceEnergyRange(true),
                     m_extrapolated(0) {}

WcsMap2::WcsMap2(const std::string & filename,
                 const std::string & extension,
                 bool interpolate,
                 bool enforceEnergyRange) 
   : m_proj(0), m_naxis3(1), m_interpolate(interpolate),
     m_isPeriodic(false), m_mapIntegral(0),
     m_enforceEnergyRange(enforceEnergyRange),
     m_extrapolated(0), m_filename(filename) {
   
   m_proj = new astro::SkyProj(filename, extension);

   const tip::Image * image = 
      tip::IFileSvc::instance().readImage(filename, extension);
   
   std::vector<float> my_image;
   image->get(my_image);

   const tip::Header & header = image->getHeader();

   header["NAXIS"].get(m_naxes);
   header["NAXIS1"].get(m_naxis1);
   header["NAXIS2"].get(m_naxis2);

   if (m_naxes == 3) {
      header["NAXIS3"].get(m_naxis3);
   }
   if (m_naxis3 > 1) {
      ExposureMap::readEnergyExtension(filename, m_energies);
      if (m_naxis3 != m_energies.size()) {
         throw std::runtime_error("NAXIS3 does not match the number of rows "
                                  "in the ENERGIES extension.");
      }
   } else {
      m_energies.push_back(100.);
   }

   header["CDELT1"].get(m_cdelt1);
   if (::my_round(m_naxis1*m_cdelt1) == 360.) {
      m_isPeriodic = true;
   }
   header["CDELT2"].get(m_cdelt2);

   header["CRPIX1"].get(m_crpix1);
   header["CRPIX2"].get(m_crpix2);

   header["CRVAL1"].get(m_crval1);
   header["CRVAL2"].get(m_crval2);

   m_crota2 = 0;

   if (m_proj->isGalactic()) {
      m_coordSys = astro::SkyDir::GALACTIC;
   } else {
      m_coordSys = astro::SkyDir::EQUATORIAL;
   }

   std::pair<double, double> coord = m_proj->pix2sph(m_crpix1, m_crpix2);
   m_refDir = astro::SkyDir(coord.first, coord.second, m_coordSys);

   try {
      std::vector<astro::SkyDir> corners;
      getCorners(corners);
      m_mapRadius = corners[0].difference(m_refDir);
      for (size_t i(1); i < corners.size(); i++) {
         double sep(corners[i].difference(m_refDir));
         if (sep > m_mapRadius) {
            m_mapRadius = sep;
         }
      }
   } catch(std::exception & eObj) {
      m_mapRadius = 180.;
   }

   delete image;

   m_image.reserve(m_naxis3);
   std::vector<float>::const_iterator pixelValue = my_image.begin();
   for (int k = 0; k < m_naxis3; k++) {
      ImagePlane_t image_plane;
      image_plane.reserve(m_naxis2);
      for (int j = 0; j < m_naxis2; j++) {
//          std::vector<double> row(m_naxis1);
         std::vector<float> row(m_naxis1);
         for (int i = 0; i < m_naxis1 && pixelValue != my_image.end();
              i++, ++pixelValue) {
            row[i] = *pixelValue;
         }
         image_plane.push_back(row);
      }
      m_image.push_back(image_plane);
   }

   computeMapIntegrals();
}

WcsMap2::WcsMap2(const DiffuseSource & diffuseSource,
                 double ra, double dec, double pix_size, int npts,
                 double energy, const std::string & proj_name, bool use_lb,
                 bool interpolate, bool enforceEnergyRange) 
   : m_refDir(ra, dec), m_interpolate(interpolate), m_mapIntegral(0),
     m_enforceEnergyRange(enforceEnergyRange),
     m_extrapolated(0), m_filename("") {
   if (use_lb) { // convert to l, b
      ra = m_refDir.l();
      dec = m_refDir.b();
   }
   double refpix = static_cast<double>(npts + 1)/2.;
   double crpix[] = {refpix, refpix};
   m_crpix1 = crpix[0];
   m_crpix2 = crpix[1];
   double crval[] = {ra, dec};
   m_crval1 = crval[0];
   m_crval2 = crval[1];
   double cdelt[] = {-pix_size, pix_size};
   m_cdelt1 = -pix_size;
   m_cdelt2 = pix_size;
   m_crota2 = 0;

   m_proj = new astro::SkyProj(proj_name, crpix, crval, cdelt, 0, use_lb);

   if (use_lb) {
      m_coordSys = astro::SkyDir::GALACTIC;
   } else {
      m_coordSys = astro::SkyDir::EQUATORIAL;
   }

   // Create a single image plane (at the specified energy).
   ImagePlane_t image_plane;
   image_plane.reserve(npts);
   double ix, iy;
   for (int j = 0; j < npts; j++) {
      iy = j + 1.;
//      std::vector<double> row(npts, 0);
      std::vector<float> row(npts, 0);
      for (int i = 0; i < npts; i++) {
         ix = i + 1.;
         if (m_proj->testpix2sph(ix, iy) == 0) {
            std::pair<double, double> coord = m_proj->pix2sph(ix, iy);
            astro::SkyDir dir(coord.first, coord.second, m_coordSys);
            SkyDirArg my_dir(dir, energy);
            row[i] = diffuseSource.spatialDist(my_dir);
         } else {
            row[i] = 0;
         }
      }
      image_plane.push_back(row);
   }
   m_image.clear();
   m_image.push_back(image_plane);
   m_naxis1 = npts;
   m_naxis2 = npts;
   m_naxis3 = 1;
   m_energies.push_back(energy);

   computeMapIntegrals();
}

WcsMap2::WcsMap2(const DiffuseSource & diffuseSource,
                 double ra, double dec, 
                 double crpix1, double crpix2, 
                 double cdelt1, double cdelt2,
                 int naxis1, int naxis2,
                 double energy, const std::string & proj_name, bool use_lb,
                 bool interpolate, bool enforceEnergyRange) 
   : m_refDir(ra, dec), 
     m_naxis1(naxis1), 
     m_naxis2(naxis2), 
     m_naxis3(1),
     m_crpix1(crpix1),
     m_crpix2(crpix2),
     m_cdelt1(cdelt1),
     m_cdelt2(cdelt2), 
     m_crota2(0),
     m_interpolate(interpolate), 
     m_mapIntegral(0),
     m_enforceEnergyRange(enforceEnergyRange),
     m_extrapolated(0),
     m_filename("") {
   if (::my_round(m_naxis1*m_cdelt1) == 360.) {
      m_isPeriodic = true;
   }
   if (use_lb) { // convert to l, b
      ra = m_refDir.l();
      dec = m_refDir.b();
   }
   double crpix[] = {crpix1, crpix2};
   double crval[] = {ra, dec};
   double cdelt[] = {cdelt1, cdelt2};

   m_proj = new astro::SkyProj(proj_name, crpix, crval, cdelt, 0, use_lb);

   if (use_lb) {
      m_coordSys = astro::SkyDir::GALACTIC;
   } else {
      m_coordSys = astro::SkyDir::EQUATORIAL;
   }

   ImagePlane_t image_plane;
   image_plane.reserve(naxis2);
   double ix, iy;
   for (int j = 0; j < naxis2; j++) {
      iy = j + 1.;
//      std::vector<double> row(naxis1, 0);
      std::vector<float> row(naxis1, 0);
      for (int i = 0; i < naxis1; i++) {
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
      image_plane.push_back(row);
   }
   check_negative_pixels(image_plane);
   m_image.clear();
   m_image.push_back(image_plane);
   m_energies.push_back(energy);
   computeMapIntegrals();
}

WcsMap2::~WcsMap2() {
// // astro::SkyProj copy constructor is not implemented properly so we
// // must ensure this pointer is not deleted here and live with the
// // resulting memory leak when this object is deleted.
//   delete m_proj;
   st_stream::StreamFormatter formatter("WcsMap2", "", 2);
   if (m_extrapolated > 0) {
      formatter.info(3) << "WcsMap2: extrapolated beyond the maximum "
                        << "energy for map cube file " << m_filename << " "
                        << m_extrapolated << " times.\n";
   }
}

WcsMap2::WcsMap2(const WcsMap2 & rhs) 
   : m_refDir(rhs.m_refDir),
     m_image(rhs.m_image), 
     m_solidAngles(rhs.m_solidAngles),
     m_naxes(rhs.m_naxes),
     m_naxis1(rhs.m_naxis1),
     m_naxis2(rhs.m_naxis2),
     m_naxis3(rhs.m_naxis3),
     m_crpix1(rhs.m_crpix1),
     m_crpix2(rhs.m_crpix2),
     m_crval1(rhs.m_crval1),
     m_crval2(rhs.m_crval2),
     m_cdelt1(rhs.m_cdelt1),
     m_cdelt2(rhs.m_cdelt2),
     m_crota2(rhs.m_crota2),
     m_energies(rhs.m_energies),
     m_interpolate(rhs.m_interpolate),
     m_isPeriodic(rhs.m_isPeriodic),
     m_coordSys(rhs.m_coordSys),
     m_mapIntegral(rhs.m_mapIntegral),
     m_enforceEnergyRange(rhs.m_enforceEnergyRange),
     m_extrapolated(rhs.m_extrapolated),
     m_filename(rhs.m_filename) {
// astro::SkyProj copy constructor is not implemented properly so we
// must share this pointer, ensure it is not deleted in the
// destructor, and live with the resulting memory leak when this
// object is deleted.
//   m_proj = new astro::SkyProj(*(rhs.m_proj));
   m_proj = rhs.m_proj;
}

WcsMap2 & WcsMap2::operator=(const WcsMap2 & rhs) {
   if (this != &rhs) {
// astro::SkyProj copy constructor is not implemented properly so we
// must share this pointer, ensure it is not deleted in the destructor,
// and live with the resulting memory leak when this object is deleted.
//      delete m_proj;
//      m_proj = new astro::SkyProj(*(rhs.m_proj));
      m_proj = rhs.m_proj;
      m_refDir = rhs.m_refDir;
      m_image = rhs.m_image;
      m_solidAngles = rhs.m_solidAngles;
      m_naxes = rhs.m_naxes;
      m_naxis1 = rhs.m_naxis1;
      m_naxis2 = rhs.m_naxis2;
      m_naxis3 = rhs.m_naxis3;
      m_crpix1 = rhs.m_crpix1;
      m_crpix2 = rhs.m_crpix2;
      m_crval1 = rhs.m_crval1;
      m_crval2 = rhs.m_crval2;
      m_cdelt1 = rhs.m_cdelt1;
      m_cdelt2 = rhs.m_cdelt2;
      m_crota2 = rhs.m_crota2;
      m_energies = rhs.m_energies;
      m_interpolate = rhs.m_interpolate;
      m_isPeriodic = rhs.m_isPeriodic;
      m_coordSys = rhs.m_coordSys;
      m_mapIntegral = rhs.m_mapIntegral;
      m_enforceEnergyRange = rhs.m_enforceEnergyRange;
      m_extrapolated = rhs.m_extrapolated;
      m_filename = rhs.m_filename;
   }
   return *this;
}

double WcsMap2::operator()(const astro::SkyDir & dir, int k) const {
   check_energy_index(k);
// NB: wcslib starts indexing pixels with 1, not 0.
   std::pair<double, double> pixel;
   try {
      pixel = dir.project(*m_proj);
   } catch (...) {
      // The annoying astro::SkyProj class throws a SkyProjException
      // but does not expose it! (It lives in the .cxx file in the
      // anonymous namespace.) So we have no choice but to catch
      // everything and assume the exception occurs because the
      // direction is outside the map.
      return 0;
   }

   double x(pixel.first);
   double y(pixel.second);

   if (m_isPeriodic) {
      x = std::fmod(x, m_naxis1);
   }

   if ((!m_isPeriodic && (x < 0.5 || x > m_naxis1 + 0.5)) ||
       y < 0.5 || y > m_naxis2 + 0.5) {
      // Sky location is outside of map, so do not extrapolate and return 0.
      return 0;
   }

   if (!m_interpolate) {
      return pixelValue(x, y, k);
   }

// This code tries to do a bilinear interpolation on the pixel values.
   int ix(static_cast<int>(x));
   int iy(static_cast<int>(y));

// For points within half a pixel of the edges of the map,
// extrapolation would be required in the context of the bilinear
// scheme.  However, this could result in unphysical negative values,
// so just return the un-interpolated value of the pixel in which the
// point lies.
   if (!m_isPeriodic) {
      if (ix < 1 || ix >= m_naxis1) {
         return pixelValue(x, y, k);
      }
   }
   if (iy < 1 || iy >= m_naxis2) {
      return pixelValue(x, y, k);
   }
   
   double tt(x - ix);
   double uu(y - iy);

   double y1, y4;

   if (m_isPeriodic && ix == 0) {
      y1 = m_image[k].at(iy-1).back();
      y4 = m_image[k].at(iy).back();
   } else {
      y1 = m_image[k].at(iy-1).at(ix-1);
      y4 = m_image[k].at(iy).at(ix-1);
   }
   double y2(m_image[k].at(iy-1).at(ix));
   double y3(m_image[k].at(iy).at(ix));
   
   double value((1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
                + tt*uu*y3 + (1. - tt)*uu*y4);

   return value;
}

double WcsMap2::operator()(const astro::SkyDir & dir, double energy) const {
   if (energy < 0) {
      energy = m_energies.front();
   }
   check_energy(energy);

   int k(0);
   if (m_naxes == 3 && m_energies.size() > 1) {
      k = std::upper_bound(m_energies.begin(), m_energies.end(), energy)
         - m_energies.begin() - 1;
      /// Extrapolate beyond highest energy.  This will only occur if
      /// m_enforceEnergyRange == false.
      if (k > static_cast<int>(m_energies.size() - 2)) {
         k = m_energies.size() - 2;
         m_extrapolated += 1;
      }
   }
   double y1 = operator()(dir, k);
   if (energy == m_energies[k]) { 
      return y1;
   }
   double y2 = operator()(dir, k+1);

   double value = interpolatePowerLaw(energy, m_energies[k],
                                      m_energies[k+1], y1, y2);
   return value;
}

WcsMap2 WcsMap2::convolve(double energy, const MeanPsf & psf,
                          const BinnedExposure & exposure,
			  const SpatialFunction& fn,
                          int k) const {

   // Convolve for a single image plane.
   check_energy_index(k);

   std::pair<astro::SkyDir, astro::SkyDir> minMaxDir;
   minMaxDir = minMaxDistPixels(fn.dir());

   const double cdelt_max = std::max(std::abs(cdelt1()),std::abs(cdelt2()));
   const double cdelt_min = std::min(std::abs(cdelt1()),std::abs(cdelt2()));
   const double theta_max = minMaxDir.second.difference(fn.dir())*180./M_PI + 2.0*cdelt_max;
   const double theta_min = cdelt_min*0.5;
   const int theta_nstep = ::round(theta_max/cdelt_min)*2;
   const double theta_step = (theta_max-theta_min)/double(theta_nstep-1);

   // Build a vector for fast interpolation
   std::vector<double> theta(theta_nstep+1);
   std::vector<double> value(theta_nstep+1);

   theta[0] = 0.0;
   value[0] = fn.spatialResponse(0.0,energy,psf);
   for(int i = 1; i < theta.size(); i++) {
       theta[i] = theta_min+(i-1)*theta_step;
       value[i] = fn.spatialResponse(theta[i],energy,psf);
   }       

   ::Image counts;
   counts.resize(m_naxis2);

   astro::SkyDir fndir =  fn.dir();
   double exp = exposure(energy, fndir.ra(), fndir.dec());

   for (int j = 0; j < m_naxis2; j++) {
      counts.at(j).resize(m_naxis1, 0);
      for (int i = 0; i < m_naxis1; i++) {
         if (m_proj->testpix2sph(i+1, j+1) == 0) {
            std::pair<double, double> coord = m_proj->pix2sph(i+1, j+1);
            astro::SkyDir dir(coord.first, coord.second, m_coordSys);

	    double dtheta = fndir.difference(dir)*180./M_PI;
	    double v = st_facilities::Util::interpolate(theta, value, 
							dtheta);
            counts[j][i] = v*exp;
         }
      }
   }

   WcsMap2 my_image(*this);
   my_image.m_image.clear();
   my_image.m_naxis3 = 1;

   my_image.m_image.push_back(counts);
   return my_image;
}

WcsMap2 WcsMap2::convolve(double energy, const MeanPsf & psf,
                          const BinnedExposure & exposure,
                          bool performConvolution,
                          int k) const {
// Convolve for a single image plane.
   check_energy_index(k);

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
               m_image[k].at(j).at(i)*exposure(energy, dir.ra(), dir.dec());
         }
      }
   }

   WcsMap2 my_image(*this);
   my_image.m_image.clear();
   my_image.m_naxis3 = 1;

   if (!performConvolution) {
      my_image.m_image.push_back(counts);
      return my_image;
   }
      
// Fill a square array with an image of the Psf at the same binning
// resolution as the source image.  Use the smaller of the image map
// dimensions to determine the psf image size.
   int npix;
   if (m_naxis1 <= m_naxis2) {
      npix = m_naxis1;
   } else {
      npix = m_naxis2;
   }

   ::Image psf_image;

   // Ensure the psf array size is odd in each dimension, so that the
   // center pixel corresponds to the center of the PSF.
   if (npix % 2 == 0) {
      npix -= 1;
   }
   double refpix = static_cast<double>(npix + 1)/2.;  // refpix at center
   double crpix[] = {refpix, refpix};
   double crval[] = {m_refDir.ra(), m_refDir.dec()}; // actually arbitrary
   double cdelt[] = {m_cdelt1, m_cdelt2}; // ensure same resolution as input map
   astro::SkyProj my_proj(m_proj->projType(), crpix, crval, cdelt);
      
   psf_image.resize(npix);
   for (int j(0); j < npix; j++) {
      psf_image.at(j).resize(npix, 0);
      for (int i(0); i < npix; i++) {
         if (my_proj.testpix2sph(i+1, j+1) == 0) {
            std::pair<double, double> coord = my_proj.pix2sph(i+1, j+1);
            astro::SkyDir dir(coord.first, coord.second, 
                              astro::SkyDir::EQUATORIAL);
            double theta = m_refDir.difference(dir)*180./M_PI;
            psf_image.at(j).at(i) = psf(energy, theta, 0);
         }
      }
   }

   psf_image.normalize();

   check_negative_pixels(counts);
   check_negative_pixels(psf_image);
   my_image.m_image.push_back(Convolve::convolve2d(counts, psf_image));

   return my_image;
}

double WcsMap2::solidAngle(const astro::SkyProj & proj, 
                           double ilon, double ilat) {
   bool atpole(false);
   std::pair<double, double> center(proj.pix2sph(ilon, ilat));
   std::pair<double, double> left(proj.pix2sph(ilon - 0.5, ilat));
   std::pair<double, double> right(proj.pix2sph(ilon + 0.5, ilat));
   std::pair<double, double> bottom(proj.pix2sph(ilon, ilat - 0.5));
   std::pair<double, double> top(proj.pix2sph(ilon, ilat + 0.5));
   if (std::fabs(center.second) == 90.) {
      // The center of this pixel is at the pole of the map, so need
      // to shift center by half a pixel in latitude and compute half
      // of the implied "bow-tie" shape for the pixel and will
      // multiply dOmega by two on output.
      left = proj.pix2sph(ilon - 0.5, ilat - 0.25);
      right = proj.pix2sph(ilon + 0.5, ilat - 0.25);
      top = proj.pix2sph(ilon, ilat);
      left = proj.pix2sph(ilon, ilat - 0.5);
      atpole = true;
   }

   astro::SkyDir rightDir(right.first, right.second);
   astro::SkyDir leftDir(left.first, left.second);
   double delta_lon = std::acos(leftDir().dot(rightDir()));

   double delta_lat = (top.second - bottom.second)*M_PI/180.;

   double dOmega = std::fabs(delta_lon*delta_lat);
   if (atpole) {
      return dOmega*2;
   }

   return dOmega;
}

double WcsMap2::solidAngle(double ilon, double ilat) const {
   try {
      return solidAngle(*m_proj, ilon, ilat);
   } catch (std::exception & eObj) {
      if (!st_facilities::Util::
          expectedException(eObj, "SkyProj wcslib error")) {
         throw;
      }
   }
   return 0;
}

//const std::vector< std::vector<double> > & WcsMap2::solidAngles() const {
const std::vector< std::vector<float> > & WcsMap2::solidAngles() const {
   if (m_solidAngles.empty()) {
//      m_solidAngles.resize(m_naxis1, std::vector<double>(m_naxis2, 0));
      m_solidAngles.resize(m_naxis1, std::vector<float>(m_naxis2, 0));
      for (int i(0); i < m_naxis1; i++) {
         for (int j(0); j < m_naxis2; j++) {
            m_solidAngles[i][j] = solidAngle(i, j);
         }
      }
   }
   return m_solidAngles;
}

double WcsMap2::pixelValue(double ilon, double ilat, int k) const {
   check_energy_index(k);

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
   return m_image[k].at(iy).at(ix);
}

astro::SkyDir WcsMap2::skyDir(double ilon, double ilat) const {
   std::pair<double, double> coords(m_proj->pix2sph(ilon, ilat));
   return astro::SkyDir(coords.first, coords.second, m_coordSys);
}

bool WcsMap2::withinMapRadius(const astro::SkyDir & dir) const {
   if (dir.difference(m_refDir) <= m_mapRadius) {
      return true;
   }
   return false;
}

bool WcsMap2::insideMap(const astro::SkyDir & dir) const {
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
WcsMap2::minMaxDistPixels(const astro::SkyDir & dir) const {
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
   i = m_naxis1;
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
   j = m_naxis2;
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

void WcsMap2::getCorners(std::vector<astro::SkyDir> & corners) const {
   corners.clear();
   corners.push_back(skyDir(1, 1));
   corners.push_back(skyDir(1, m_naxis2));
   corners.push_back(skyDir(m_naxis1, m_naxis2));
   corners.push_back(skyDir(m_naxis1, 1));
}

double WcsMap2::mapIntegral() const {
   return m_mapIntegral;
}

double WcsMap2::mapIntegral(double energy) const {
   check_energy(energy);

   int k(0);
   if (m_naxes == 3) {
      k = std::upper_bound(m_energies.begin(), m_energies.end(), energy)
         - m_energies.begin() - 1;
   }
   if (energy == m_energies.at(k)) {
      return m_mapIntegrals.at(k);
   }
   if (k == m_energies.size() - 1) {
      k = m_energies.size() - 2;
   }

   double value = (m_mapIntegrals[k]*
                   std::exp((std::log(energy/m_energies[k]))
                            /(std::log(m_energies[k+1]/m_energies[k]))
                            *std::log(m_mapIntegrals[k+1]/m_mapIntegrals[k])));
   return value;
}

void WcsMap2::computeMapIntegrals() {
   m_mapIntegrals.clear();
   for (int k(0); k < m_naxis3; k++) {
      m_mapIntegrals.push_back(0);
      for (int j(0); j < m_naxis2; j++) {
         for (int i(0); i < m_naxis1; i++) {
            // NB: Indexing for solidAngles() is reversed from usual
            // convention.
            m_mapIntegrals.back() += solidAngles()[i][j]*m_image[k][j][i];
         }
      }
   }

// For a 2D map, the map integral is just the angle-integrated map.
   if (m_naxis3 == 1) {
      m_mapIntegral = m_mapIntegrals.at(0);
      return;
   }

// 3D maps: integrate the angle-integrate maps over energy.
   m_mapIntegral = 0;
   for (size_t k(0); k < m_naxis3-1; k++) {
      m_mapIntegral += ( (m_mapIntegrals[k+1]*m_energies[k+1] +
                          m_mapIntegrals[k]*m_energies[k])/2.
                         *std::log(m_energies[k+1]/m_energies[k]) );
   }
}

WcsMap2 * WcsMap2::rebin(unsigned int factor, bool average) {
   WcsMap2 * my_map = new WcsMap2(*this);
   int dnxp = factor - (m_naxis1 % factor);
   if (dnxp == factor) {
      dnxp = 0;
   }
   int dnyp = factor - (m_naxis2 % factor);
   if (dnyp == factor) {
      dnyp = 0;
   }
   my_map->m_naxis1 = (m_naxis1 + dnxp)/factor;
   my_map->m_naxis2 = (m_naxis2 + dnyp)/factor;

   // Set reference pixel, keeping same reference direction.
   my_map->m_crpix1 = (m_crpix1 - 0.5)/factor + 0.5;
   my_map->m_crpix2 = (m_crpix2 - 0.5)/factor + 0.5;

   // apply the rebinning factor to the pixel size at the reference direction.
   my_map->m_cdelt1 = m_cdelt1*factor;
   my_map->m_cdelt2 = m_cdelt2*factor;

   // Set the projection (assuming the astro::SkyProj destructor is
   // still not implemented correctly so that we cannot delete the
   // pointer to it).
   double crpix[2] = {my_map->m_crpix1, my_map->m_crpix2};
   double cdelt[2] = {my_map->m_cdelt1, my_map->m_cdelt2};
   double crval[2];
   if (m_proj->isGalactic()) {
      crval[0] = m_refDir.l();
      crval[1] = m_refDir.b();
   } else {
      crval[0] = m_refDir.ra();
      crval[1] = m_refDir.dec();
   }
   my_map->m_proj = new astro::SkyProj(m_proj->projType(), crpix, crval, cdelt,
                                       m_crota2, m_proj->isGalactic());

   st_stream::StreamFormatter formatter("WcsMap2", "", 2);

   formatter.info(4) << "naxis1: " << my_map->m_naxis1 << "\n"
                     << "naxis2: " << my_map->m_naxis2 << "\n";

   formatter.info(4) << "crpix1: " << crpix[0] << "\n"
                     << "crpix2: " << crpix[1] << "\n";

   formatter.info(4) << "cdelt1: " << cdelt[0] << "\n"
                     << "cdelt2: " << cdelt[1] << "\n";

   formatter.info(4) << "crval1: " << crval[0] << "\n"
                     << "crval2: " << crval[1] << "\n";

   /// Fill array with solid angle sums used for averaging.
//   std::vector< std::vector<double> > my_solidAngles;
   std::vector< std::vector<float> > my_solidAngles;
   if (average) {
      for (int j(0); j < my_map->m_naxis2; j++) {
//         my_solidAngles.push_back(std::vector<double>(my_map->m_naxis1, 0));
         my_solidAngles.push_back(std::vector<float>(my_map->m_naxis1, 0));
      }
      for (size_t i(0); i < m_naxis1; i++) {
         unsigned int ii = i/factor;
         for (size_t j(0); j < m_naxis2; j++) {
            unsigned int jj = j/factor;
            // Note that solidAngles() has opposite ordering of indexes.
            my_solidAngles[jj][ii] += solidAngles()[i][j];
         }
      }
   }

   for (int k(0); k < m_naxis3; k++) {
      my_map->m_image[k].clear();
      my_map->m_image[k].resize(my_map->m_naxis2);
      for (size_t j(0); j < my_map->m_naxis2; j++) {
         my_map->m_image[k][j].resize(my_map->m_naxis1, 0);
      }

      for (size_t i(0); i < m_naxis1; i++) {
         unsigned int ii = i/factor;
         for (size_t j(0); j < m_naxis2; j++) {
            unsigned int jj = j/factor;
            if (average) {
               my_map->m_image[k][jj][ii] += 
                  m_image[k][j][i]*solidAngles()[i][j];
            } else {
               my_map->m_image[k][jj][ii] += m_image[k][j][i];
            }
         }
      }
      if (average) {
         for (size_t ii(0); ii < my_map->m_naxis1; ii++) {
            for (size_t jj(0); jj < my_map->m_naxis2; jj++) {
               my_map->m_image[k][jj][ii] /= my_solidAngles[jj][ii];
            }
         }
      }
   }
   my_map->m_solidAngles.clear();

   my_map->computeMapIntegrals();

   return my_map;
}

void WcsMap2::check_energy_index(int k) const {
   if (m_naxes == 3 && (k < 0 || k > m_naxis3)) {
      throw std::runtime_error("WcsMap2: Requested energy index is "
                               "out-of-range.");
   }
}

void WcsMap2::check_energy(double energy) const {
   if (m_enforceEnergyRange && 
       (energy < m_energies.front() || energy > m_energies.back())) {
      throw std::runtime_error("WcsMap2: Requested energy is out-of-range.");
   }
}

void WcsMap2::setExtrapolation(bool enforceEnergyRange) {
   m_enforceEnergyRange = enforceEnergyRange;
}

bool WcsMap2::enforceEnergyRange() const {
   return m_enforceEnergyRange;
}

void WcsMap2::check_negative_pixels(const ImagePlane_t & image) const {
   for (size_t j(0); j < image.size(); j++) {
      for (size_t i(0); i < image[j].size(); i++) {
         if (image[j][i] < 0) {
            throw std::runtime_error("Image pixel value less than zero.");
         }
      }
   }
}

double WcsMap2::interpolatePowerLaw(double x, double x1, double x2,
                                    double y1, double y2) {
   if (y1 == 0 && y2 == 0) {
      return 0;
   }
   if (y1 == 0 || y2 == 0) {
      if ( (x < x1 && y1 == 0) || (x > x2 && y2 == 0) ) {
         // Linear extrapolation in these cases would produce
         // negative values.
         std::ostringstream message;
         message << "WcsMap2::interpolatePowerLaw: "
                 << "linear extrapolation selected for "
                 << "zero-valued ordinates.\n"
                 << "x = " << x << ", "
                 << "x1 = " << x1 << ", "
                 << "x2 = " << x2 << ", "
                 << "y1 = " << y1 << ", "
                 << "y2 = " << y2 << std::endl;
         throw std::runtime_error(message.str());
      }
      // Use linear interpolation
      double value((x - x1)/(x2 - x1)*(y2 - y1) + y1);
      return value;
   }
   if (x1 <= 0 || x2 <= 0 || y1 <= 0 || y2 <= 0) {
      std::ostringstream message;
      message << "WcsMap2::interpolatePowerLaw:\n"
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
