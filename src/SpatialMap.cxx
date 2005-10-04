/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.16 2005/10/03 15:02:43 jchiang Exp $
 *
 */

#include <memory>
#include <stdexcept>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Header.h"

#include "astro/SkyProj.h"

#include "st_facilities/Util.h"
#include "st_facilities/FitsImage.h"

#include "Likelihood/SkyDirArg.h"
#include "Likelihood/SpatialMap.h"

namespace Likelihood {

SpatialMap::SpatialMap(const SpatialMap & rhs) 
   : optimizers::Function(rhs) {
   m_proj = new astro::SkyProj(*(rhs.m_proj));
   m_fitsFile = rhs.m_fitsFile;
   m_extension = rhs.m_extension;
   m_naxis1 = rhs.m_naxis1;
   m_naxis2 = rhs.m_naxis2;
   m_naxis3 = rhs.m_naxis3;
   m_image = rhs.m_image;
}

SpatialMap & SpatialMap::operator=(const SpatialMap & rhs) {
   if (this != &rhs) {
      m_proj = new astro::SkyProj(*(rhs.m_proj));
      optimizers::Function::operator=(rhs);
      m_fitsFile = rhs.m_fitsFile;
      m_extension = rhs.m_extension;
      m_naxis1 = rhs.m_naxis1;
      m_naxis2 = rhs.m_naxis2;
      m_naxis3 = rhs.m_naxis3;
      m_image = rhs.m_image;
      delete m_proj;
   }
   return *this;
}

SpatialMap::~SpatialMap() {
   delete m_proj;
}

void SpatialMap::init() {
// This Function has one Parameter, an overall normalization, 
// but set it to be unit constant.
   int nParams = 1;
   setMaxNumParams(nParams);
   m_genericName = "SpatialMap";
   addParam("Prefactor", 1, false);
   setParamAlwaysFixed("Prefactor");
}

void SpatialMap::readFitsFile(const std::string & fitsFile,
                              const std::string & extension) {
   m_fitsFile = fitsFile;
   m_extension = extension;

   facilities::Util::expandEnvVar(&m_fitsFile);

   m_proj = st_facilities::FitsImage::skyProjCreate(m_fitsFile, extension);

   const tip::Image * image = 
      tip::IFileSvc::instance().readImage(m_fitsFile, extension);
   
   image->get(m_image);
   
   const tip::Header & header = image->getHeader();
   
   header["NAXIS1"].get(m_naxis1);
   header["NAXIS2"].get(m_naxis2);
   m_naxis3 = 1;
   try {
      header["NAXIS3"].get(m_naxis3);
   } catch (...) {
   }

   delete image;
}

double SpatialMap::value(optimizers::Arg & arg) const {
   astro::SkyDir dir;
   dynamic_cast<SkyDirArg &>(arg).fetchValue(dir);
   
   std::pair<double, double> pixels = dir.project(*m_proj);
   if (pixels.first < 0 || pixels.first > m_naxis1 ||
       pixels.second < 0 || pixels.second > m_naxis2) {
      return 0;
   }

   double x, y;
   if (m_proj->isGalactic()) {
      x = dir.l();
      y = dir.b();
   } else {
      x = dir.ra();
      y = dir.dec();
   }

   size_t ix = static_cast<size_t>(pixels.first);
   size_t iy = static_cast<size_t>(pixels.second);

   std::pair<double, double> lower_left = m_proj->pix2sph(ix, iy);
   std::pair<double, double> upper_right = m_proj->pix2sph(ix+1, iy+1);

   double uu((x - lower_left.first)/(upper_right.first - lower_left.first));
   double tt((y - lower_left.second)/(upper_right.second - lower_left.second));

// NB: wcslib starts indexing pixels with 1, not 0; so we need to 
// apply correction for the off-by-one error here.
   double y1(m_image.at(m_naxis1*(iy-1) + (ix-1)));
   double y2(m_image.at(m_naxis1*iy + ix - 1));
   double y3(m_image.at(m_naxis1*iy + ix));
   double y4(m_image.at(m_naxis1*(iy-1) + ix));

   double value((1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
                + tt*uu*y3 + (1. - tt)*uu*y4);

   return value;
}

} // namespace Likelihood
