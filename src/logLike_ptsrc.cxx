/** @file logLike_ptsrc.cxx
 * @brief logLike_ptsrc class implementation
 *
 * $Header:
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "logLike_ptsrc.h"

namespace Likelihood {

/* compute the EML log-likelihood for a single-point source */

double logLike_ptsrc::value(const std::vector<double> &paramVec) {
   setParamValues(paramVec);
   
   double my_value = 0;
   
   for (int j = 0; j < m_eventData[0].dim; j++) {
      double src_sum = 0.;
      for (unsigned int i = 0; i < getNumSrcs(); i++)
// NB: Here evaluate_at(double) is inherited from SourceModel and
// evaluates as a function of the data variable.  This will need
// generalization in SourceModel, i.e., an Event class object should
// be passed instead.
         src_sum += evaluate_at(m_eventData[0].val[j]);
      my_value += log(src_sum);
   }
   for (unsigned int i = 0; i < getNumSrcs(); i++) {
      Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++)
         my_value -= (*func_it).second->integral(-1e3, 1e3);
   }
   
   return my_value;
}

} // namespace Likelihood
