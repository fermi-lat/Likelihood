/** @file logLike_gauss.cxx
 * @brief logLike_gauss class implementation
 *
 * $Header:
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "logLike_gauss.h"
#include "../Likelihood/dArg.h"

namespace Likelihood {

//! objective function as a function of the free parameters
//! do EML (specialized to Gaussian functions for now)
double logLike_gauss::value(const std::vector<double> &paramVec) {
   setParamValues(paramVec);
   
   double my_value = 0;
   
   for (int j = 0; j < m_eventData[0].dim; j++) {
      double src_sum = 0.;
      for (unsigned int i = 0; i < getNumSrcs(); i++) {
// NB: Here evaluate_at(Arg) is inherited from SourceModel and
// evaluates as a function of the data variable. 
	 dArg xarg(m_eventData[0].val[j]);
         src_sum += evaluate_at(xarg);
      }
      my_value += log(src_sum);
   }
   for (unsigned int i = 0; i < getNumSrcs(); i++) {
      Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      dArg xmin(-1e3);
      dArg xmax(1e3);
      for (; func_it != srcFuncs.end(); func_it++) {
         my_value -= (*func_it).second->integral(xmin, xmax);
      }
   }
   
   return my_value;
}

} // namespace Likelihood
