/** 
 * @file SpatialMap.cxx
 * @brief Implementation of Function object class that returns interpolated
 * image values of a FITS image file.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SpatialMap.cxx,v 1.9 2003/11/30 23:14:40 jchiang Exp $
 *
 */

#include <algorithm>
#include <numeric>

#include "facilities/Util.h"

#include "Likelihood/FitsImage.h" 
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/SpatialMap.h"

namespace {
   double bilinear(const std::vector<double> &xx, double x, 
                   const std::vector<double> &yy, double y, 
                   const std::valarray<double> &z) 
      throw(Likelihood::Exception) {

      std::vector<double>::const_iterator ix;
      if (x < *(xx.begin())) {
         ix = xx.begin() + 1;
      } else if (x >= *(xx.end()-1)) {
         ix = xx.end() - 1;
      } else {
         ix = std::upper_bound(xx.begin(), xx.end(), x);
      }
      int i = ix - xx.begin();

      std::vector<double>::const_iterator iy;
      if (y < *(yy.begin())) {
         iy = yy.begin() + 1;
      } else if (y >= *(yy.end()-1)) {
         iy = yy.end() - 1;
      } else {
         iy = std::upper_bound(yy.begin(), yy.end(), y);
      }
      int j = iy - yy.begin();

      double tt = (x - *(ix-1))/(*(ix) - *(ix-1));
      double uu = (y - *(iy-1))/(*(iy) - *(iy-1));

      double y1 = z[yy.size()*(i-1) + (j-1)];
      double y2 = z[yy.size()*(i) + (j-1)];
      double y3 = z[yy.size()*(i) + (j)];
      double y4 = z[yy.size()*(i-1) + (j)];

      double value = (1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
         + tt*uu*y3 + (1. - tt)*uu*y4; 
      if (value < 0.) {
         std::cerr << "bilinear: value = " << value << " < 0\n";
         std::cerr << xx[i-1] << "  " << *(ix-1) << "  " 
                   << x << "  " << *ix << "\n";
         std::cerr << yy[j-1] << "  " << *(iy-1) << "  " 
                   << y << "  " << *iy << "\n";
         std::cerr << tt << "  " << uu << "  " 
                   << y1 << "  " << y2 << "  "
                   << y3 << "  " << y4 << "  " << std::endl;
         throw Likelihood::Exception("SpatialMap::bilinear: \n bailing...");
      }
      return value;
   }
} // unnamed namespace

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

   std::string inFile(m_fitsFile);
   facilities::Util::expandEnvVar(&inFile);

   FitsImage fitsImage(inFile);

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

   double my_value = ::bilinear(m_dec, dir.dec(), m_ra, ra, m_image);
      
   return m_parameter[0].getTrueValue()*my_value;
}

} // namespace Likelihood
