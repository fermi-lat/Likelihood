/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.5 2003/08/13 18:01:16 jchiang Exp $
 *
 */

#include "Likelihood/FitsImage.h" 
#include "Likelihood/Response.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/SpatialMap.h"

namespace Likelihood {

void SpatialMap::init() {
// This Function has one Parameter, an overall normalization, 
// but set it to be unit constant.
   int nParams = 1;
   setMaxNumParams(nParams);
   m_genericName = "SpatialMap";
   addParam("Prefactor", 1, false);
}

void SpatialMap::readFitsFile(const std::string &fitsFile) {
   m_fitsFile = fitsFile;

   FitsImage fitsImage(fitsFile);

// Assume 0th and 1st axes are RA and DEC.
   fitsImage.fetchAxisVector(0, m_ra);
   fitsImage.fetchAxisVector(1, m_dec);

   fitsImage.fetchImageData(m_image);
}

double SpatialMap::value(optimizers::Arg& arg) const {
   astro::SkyDir dir;
   dynamic_cast<SkyDirArg &>(arg).fetchValue(dir);

   double ra = dir.ra();

// wrap to +/-180
   if (ra > 180) ra = ra - 360;

   double my_value 
      = Response::bilinear(m_dec, dir.dec(), m_ra, ra, m_image);

   return m_parameter[0].getTrueValue()*my_value;
}

} // namespace Likelihood
