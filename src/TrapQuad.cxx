/** 
 * @file TrapQuad.cxx
 * @brief Implementation of the TrapQuad class, which performs simple 1D 
 * trapezoidal quadrature.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/TrapQuad.cxx,v 1.14 2009/03/27 19:24:32 jchiang Exp $
 */

#include "Likelihood/TrapQuad.h"
#include "optimizers/dArg.h"
#include "optimizers/Exception.h"

#include <sstream>

namespace Likelihood {

double TrapQuad::integral() {
   if (m_func) {
      std::ostringstream errorMessage;
      errorMessage << "TrapQuad::integral:\n"
                   << "This operation requires "
                   << "object creation with abscissa and ordinate vectors.\n";
      throw optimizers::Exception(errorMessage.str());
   }
   double value;
   if (m_useLog) {
      value = compute_log_integral();
   } else {
      value = compute_integral();
   }
   return value;
}

double TrapQuad::integral(double xmin, double xmax, int npts) {
   if (!m_func) {
      std::ostringstream errorMessage;
      errorMessage << "TrapQuad::integral:\n"
                   << "This operation requires "
                   << "object creation with a pointer to a Function object.\n";
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
   double value;
   if (m_useLog) {
      value = compute_log_integral();
   } else {
      value = compute_integral();
   }
   return value;
}

double TrapQuad::integral(const std::vector<double> &xvals) {
   if (!m_func) {
      std::ostringstream errorMessage;
      errorMessage << "TrapQuad::integral:\n"
                   << "This operation requires "
                   << "object creation with a pointer to a Function object.\n";
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
   double value;
   if (m_useLog) {
      value = compute_log_integral();
   } else {
      value = compute_integral();
   }
   return value;
}

double TrapQuad::compute_integral() {
   int n = m_y.size();
   double sum = m_y[0]*(m_x[1] - m_x[0])/2. 
      + m_y[n-1]*(m_x[n-1] - m_x[n-2])/2.;
   n--;
   while (--n > 0) sum += m_y[n]*(m_x[n+1] - m_x[n-1])/2.;
   return sum;
}

double TrapQuad::compute_log_integral() {
//    double sum(0);
//    for (unsigned int i = 0; i < m_x.size()-1; i++) {
//       if (m_y[i+1] > 0 && m_y[i] > 0) {
//          double b = log(m_y[i+1]/m_y[i])/log(m_x[i+1]/m_x[i]);
//          if (b != -1.) {
//             sum += m_y[i]/(b + 1.)*(m_x[i+1]*pow(m_x[i+1]/m_x[i], b) - m_x[i]);
//          } else {
//             double a = m_y[i]/pow(m_x[i], b);
//             sum += a*log(m_x[i+1]/m_x[i]);
//          }
//       } else {
//          sum += (m_y[i+1] + m_y[i])/2.*(m_x[i+1] - m_x[i]);
//       }
//    }
//    return sum;
   int n = m_y.size();
   double sum = m_y[0]*m_x[0]*std::log(m_x[1]/m_x[0])/2. 
      + m_y[n-1]*m_x[n-1]*std::log(m_x[n-1]/m_x[n-2])/2.;
   n--;
   while (--n > 0) sum += m_y[n]*m_x[n]*std::log(m_x[n+1]/m_x[n-1])/2.;
   return sum;
}

} // namespace Likelihood
