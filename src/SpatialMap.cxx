/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.14 2005/02/17 23:22:32 jchiang Exp $
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

   m_coordSys = fitsImage.coordSys();

// Assume 0th and 1st axes are RA and DEC.
   fitsImage.getAxisVector(0, m_ra);
// wrap to +/- 180
   for (unsigned int i = 0; i < m_ra.size(); i++) {
      if (m_ra.at(i) > 180) {
         m_ra.at(i) -= 360;
      }
   }
   if (m_ra.front() < m_ra.back()) {
      m_raMin = m_ra.front();
      m_raMax = m_ra.back();
   } else {
      m_raMin = m_ra.back();
      m_raMax = m_ra.front();
   }
   fitsImage.getAxisVector(1, m_dec);
   if (m_dec.front() < m_dec.back()) {
      m_decMin = m_dec.front();
      m_decMax = m_dec.back();
   } else {
      m_decMin = m_dec.back();
      m_decMax = m_dec.front();
   }
   fitsImage.getImageData(m_image);
}

double SpatialMap::value(optimizers::Arg& arg) const {
   astro::SkyDir dir;
   dynamic_cast<SkyDirArg &>(arg).fetchValue(dir);

   double ra, dec;
   if (m_coordSys == "Equatorial") {
      ra = dir.ra();
      dec = dir.dec();
   } else {
      ra = dir.l();
      dec = dir.b();
   }

// wrap to +/-180
   if (ra > 180) ra -= 360;

   if (dec < m_decMin || dec > m_decMax || ra < m_raMin || ra > m_raMax) {
      return 0;
   }
   double my_value = 
      st_facilities::Util::bilinear(m_dec, dec, m_ra, ra, m_image);
      
   return m_parameter[0].getTrueValue()*my_value;
}

} // namespace Likelihood
