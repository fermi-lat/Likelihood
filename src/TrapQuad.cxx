#include "../Likelihood/TrapQuad.h"
#include "../Likelihood/dArg.h"

namespace Likelihood {

double TrapQuad::integral() {
   if (m_haveFunc) {
      std::cerr << "TrapQuad::integral:  This operation requires "
                << "instantiation with abscissa and ordinate vectors." 
                << std::endl;
      return 0;
   }
   return compute_integral();
}

double TrapQuad::integral(double xmin, double xmax, int npts) {
   if (!m_haveFunc) {
      std::cerr << "TrapQuad::integral:  This operation requires "
                << "instantiation with a pointer to a Function object."
                << std::endl;
      return 0;
   }
   m_x.reserve(npts);
   m_y.reserve(npts);
   double xstep = (xmax - xmin)/(npts - 1);
   for (int i = 0; i < npts; i++) {
      m_x[i] = i*xstep + xmin;
      dArg xarg(m_x[i]);
      m_y[i] = (*m_func)(xarg);
   }
   return compute_integral();
}

double TrapQuad::integral(std::vector<double> &xvals) {
   if (!m_haveFunc) {
      std::cerr << "TrapQuad::integral:  This operation requires "
                << "instantiation with a pointer to a Function object."
                << std::endl;
      return 0;
   }
   m_x = xvals;
   int npts = m_x.size();
   m_y.reserve(npts);
   for (int i = 0; i < npts; i++) {
      dArg xarg(m_x[i]);
      m_y[i] = (*m_func)(xarg);
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
