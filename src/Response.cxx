/** 
 * @file Response.cxx
 * @brief Implementation for Response base class.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Response.cxx,v 1.15 2003/05/29 20:10:46 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "Likelihood/Response.h"
#include "Likelihood/Table.h"
#include "LikelihoodException.h"

namespace Likelihood {

// default maximum inclination in degrees for strawman response files
double Response::s_incMax = 70.;

Response::Response() throw(LikelihoodException) {
   ScData *scData = ScData::instance();
   if (!scData) {
      throw LikelihoodException("*scData not properly instantiated!");
   }
}

double Response::bilinear(const std::vector<double> &xx, double x, 
                          const std::vector<double> &yy, double y, 
                          const std::valarray<double> &z) 
   throw(LikelihoodException) {

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
      throw LikelihoodException("Response::bilinear: \n bailing...");
   }
   return value;
}

} // namespace Likelihood
