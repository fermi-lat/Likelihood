/** @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header$
 *
 */

#include "Likelihood/Response.h"
#include "Likelihood/SpatialMap.h"

namespace Likelihood {

SpatialMap::SpatialMap(std::string fitsfile) : FitsImage(fitsfile) {

// Assume 0th and 1st axes are RA and DEC.
   fetchAxisVector(0, m_ra);
   fetchAxisVector(1, m_dec);

// This Function has one Parameter, an overall normalization, 
// but set it to be unit constant.
   int nParams = 1;
   setMaxNumParams(nParams);
   addParam("Prefactor", 1, false);
}

double SpatialMap::value(Arg& arg) const {
   astro::SkyDir dir;
   dynamic_cast<SkyDirArg &>(arg).fetchValue(dir);

   double ra = dir.ra();

// wrap to +/-180
   if (ra > 180) ra = ra - 360;

   double my_value 
      = Response::bilinear(m_dec, dir.dec(), m_ra, ra, m_image);

   return my_value;
}

} // namespace Likelihood
