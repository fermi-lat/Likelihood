#include <vector>
#include <string>
#include <cmath>

#include "../Likelihood/Response.h"
#include "../Likelihood/Table.h"

namespace Likelihood {

Response::Response() {
   ScData *scData = ScData::instance();
}

void Response::m_hunt(double *xx, int n, double x, int *jlo) {
   int jm,jhi,inc;
   int ascnd;
   
   ascnd=(xx[n] > xx[1]);
   if (*jlo <= 0 || *jlo > n) {
      *jlo=0;
      jhi=n+1;
   } else {
      inc=1;
      if (x >= xx[*jlo] == ascnd) {
	 if (*jlo == n) return;
	 jhi=(*jlo)+1;
	 while (x >= xx[jhi] == ascnd) {
	    *jlo=jhi;
	    inc += inc;
	    jhi=(*jlo)+inc;
	    if (jhi > n) {
	       jhi=n+1;
	       break;
	    }
	 }
      } else {
	 if (*jlo == 1) {
	    *jlo=0;
	    return;
	 }
	 jhi=(*jlo)--;
	 while (x < xx[*jlo] == ascnd) {
	    jhi=(*jlo);
	    inc <<= 1;
	    if (inc >= jhi) {
	       *jlo=0;
	       break;
	    }
	    else *jlo=jhi-inc;
	 }
      }
   }
   while (jhi-(*jlo) != 1) {
      jm=(jhi+(*jlo)) >> 1;
      if (x > xx[jm] == ascnd)
	 *jlo=jm;
      else
	 jhi=jm;
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
      std::cout << "bilinear: value = " << value << " <0\n";
      std::cout << i << "  " << xx[i] << "  " << x << "  " << xx[i+1] << "\n";
      std::cout << j << "  " << yy[j] << "  " << y << "  " << yy[j+1] << "\n";
      std::cout << tt << "  " << uu << "  " 
                << y1 << "  " << y2 << "  "
                << y3 << "  " << y4 << "  " << std::endl;
   }
   return value;
}


} // namespace Likelihood
