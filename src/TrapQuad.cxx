/** 
 * @file TrapQuad.cxx
 * @brief Implementation of the TrapQuad class, which performs simple 1D 
 * trapezoidal quadrature.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/TrapQuad.cxx,v 1.9 2003/07/21 22:14:58 jchiang Exp $
 */

#include "Likelihood/TrapQuad.h"
#include "optimizers/dArg.h"
#include "optimizers/Exception.h"

#include <sstream>

namespace Likelihood {

double TrapQuad::integral() throw(optimizers::Exception) {
   if (m_haveFunc) {
      std::ostringstream errorMessage;
      errorMessage << "TrapQuad::integral:\n"
                   << "This operation requires "
                   << "instantiation with abscissa and ordinate vectors.\n";
      throw optimizers::Exception(errorMessage.str());
   }
   return compute_integral();
}

double TrapQuad::integral(double xmin, double xmax, int npts) 
   throw(optimizers::Exception) {
   if (!m_haveFunc) {
      std::ostringstream errorMessage;
      errorMessage << "TrapQuad::integral:\n"
                   << "This operation requires "
                   << "instantiation with a pointer to a Function object.\n";
      throw optimizers::Exception(errorMessage.str());
   }
   m_x.clear();
   m_y.clear();
   m_x.reserve(npts);
   m_y.reserve(npts);
   double xstep = (xmax - xmin)/(npts - 1);
   for (int i = 0; i < npts; i++) {
      m_x.push_back(i*xstep + xmin);
      optimizers::dArg xarg(m_x[i]);
      m_y.push_back((*m_func)(xarg));
   }
   return compute_integral();
}

double TrapQuad::integral(std::vector<double> &xvals) 
   throw(optimizers::Exception) {
   if (!m_haveFunc) {
      std::ostringstream errorMessage;
      errorMessage << "TrapQuad::integral:\n"
                   << "This operation requires "
                   << "instantiation with a pointer to a Function object.\n";
      throw optimizers::Exception(errorMessage.str());
   }
   m_x = xvals;
   int npts = m_x.size();
   m_y.clear();
   m_y.reserve(npts);
   for (int i = 0; i < npts; i++) {
      optimizers::dArg xarg(m_x[i]);
      m_y.push_back((*m_func)(xarg));
   }
   return compute_integral();
}

double TrapQuad::compute_integral() {
   int n = m_y.size();
   double sum = m_y[0]*(m_x[1] - m_x[0])/2. 
      + m_y[n-1]*(m_x[n-1] - m_x[n-2])/2.;
   n--;
   while (--n > 0) sum += m_y[n]*(m_x[n+1] - m_x[n-1])/2.;
   return sum;
}

} // namespace Likelihood
