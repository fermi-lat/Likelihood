/** @file TrapQuad.h
 * @brief Declaration of the TrapQuad class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef TrapQuad_h
#define TrapQuad_h

#include "Likelihood/Function.h"

namespace Likelihood {

/** 
 * @class TrapQuad
 *
 * @brief This class performs simple 1D trapezoidal quadrature.  It is
 * instantiated either using (a) two vectors of abscissa and ordinate
 * values or (b) a Function object.  For the latter, the abscissa
 * values are specified in the call to the integral() member function
 * either by a vector or by minimum and maximum values and the number
 * of abscissa points, which will be uniformly spaced.
 *
 * @authors J. Chiang
 *    
 * $Header$
 */

class TrapQuad {
    
public:
   
   TrapQuad(const std::vector<double> &x, const std::vector<double> &y) :
      m_x(x), m_y(y) {m_haveFunc = false;}
   TrapQuad(Function *func) {
      m_func = func;
      m_haveFunc = true;
   }
   ~TrapQuad() {}

   double integral();
   double integral(double xmin, double xmax, int npts = 100);
   double integral(std::vector<double> &xvals);

private:

   std::vector<double> m_x;
   std::vector<double> m_y;

   bool m_haveFunc;
   Function *m_func;

   double compute_integral();

};

} // namespace Likelihood

#endif // TrapQuad_h
