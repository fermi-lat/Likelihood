/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.12 2004/09/28 04:32:25 jchiang Exp $
 *
 */

#include <algorithm>
#include <numeric>

#include "facilities/Util.h"

#include "st_facilities/Util.h"
#include "st_facilities/FitsImage.h"

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
   setParamAlwaysFixed("Prefactor");
}

void SpatialMap::readFitsFile(const std::string &fitsFile) {
   m_fitsFile = fitsFile;

   std::string inFile(m_fitsFile);
   facilities::Util::expandEnvVar(&inFile);

   FitsImage fitsImage(inFile);

// Assume 0th and 1st axes are RA and DEC.
   fitsImage.getAxisVector(0, m_ra);
   fitsImage.getAxisVector(1, m_dec);

   fitsImage.getImageData(m_image);
}

double SpatialMap::value(optimizers::Arg& arg) const {
   astro::SkyDir dir;
   dynamic_cast<SkyDirArg &>(arg).fetchValue(dir);

   double ra = dir.ra();

// wrap to +/-180
   if (ra > 180) ra = ra - 360;

//    if (dir.dec() < m_dec.front() || dir.dec() > m_dec.back()
//        || dir.ra() < m_ra.front() || dir.ra() > m_ra.back()) {
//       return 0;
//    }
   double my_value = 
      st_facilities::Util::bilinear(m_dec, dir.dec(), m_ra, ra, m_image);
      
   return m_parameter[0].getTrueValue()*my_value;
}

} // namespace Likelihood
