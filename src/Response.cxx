/** @file Response.cxx
 * @brief Implementation for Response base class.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Response.cxx,v 1.9 2003/03/17 00:53:44 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>

#include "Likelihood/Response.h"
#include "Likelihood/Table.h"

namespace Likelihood {

// default maximum inclination in degrees for strawman response files
double Response::s_incMax = 70.;

Response::Response() {
   ScData *scData = ScData::instance();
   if (!scData) {
      std::cerr << "*scData not properly instantiated!" << std::endl;
      exit(0);
   }
}

double Response::m_bilinear(int nx, double *xx, int i, double x, 
                            int ny, double *yy, int j, double y, double *z) {

/* be sure to pass xx as xx-1 and yy as yy-1 to account for NR
   unit offset kludge */

   double tt, uu, y1, y2, y3, y4, value;

   tt = (x - xx[i])/(xx[i+1] - xx[i]);
   uu = (y - yy[j])/(yy[j+1] - yy[j]);

   y1 = z[ny*(i-1) + (j-1)];
   y2 = z[ny*i + (j-1)];
   y3 = z[ny*i + j];
   y4 = z[ny*(i-1) + j];

   value = (1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
      + tt*uu*y3 + (1. - tt)*uu*y4; 
   if (value < 0.) {
      std::cerr << "bilinear: value = " << value << " <0\n";
      std::cerr << i << "  " << xx[i] << "  " << x << "  " << xx[i+1] << "\n";
      std::cerr << j << "  " << yy[j] << "  " << y << "  " << yy[j+1] << "\n";
      std::cerr << tt << "  " << uu << "  " 
                << y1 << "  " << y2 << "  "
                << y3 << "  " << y4 << "  " << std::endl;
   }
   return value;
}

} // namespace Likelihood
