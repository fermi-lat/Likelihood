/** @file logSrcModel.cxx
 * @brief logSrcModel class implementation
 *
 * $Header:
 */

#include <cmath>
#include "../Likelihood/logSrcModel.h"

namespace Likelihood {

double logSrcModel::value(Arg &xarg) const {
   Event evt;
   dynamic_cast<EventArg &>(xarg).fetchValue(evt);

   double my_value = 0;
   for (unsigned int i = 0; i < getNumSrcs(); i++) {
      my_value += m_sources[i]->fluxDensity(evt);
   }
   return my_value;
}

void logSrcModel::fetchDerivs(Arg &xarg, std::vector<double> &derivs, 
                              bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   Event evt;
   dynamic_cast<EventArg &>(xarg).fetchValue(evt);
   double srcSum = value(xarg);

   for (unsigned int i = 0; i < m_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<std::string> paramNames;
         if (getFree) {
            (*func_it).second->getFreeParamNames(paramNames);
         } else {
            (*func_it).second->getParamNames(paramNames);
         }
         for (unsigned int j = 0; j < paramNames.size(); j++) {
            derivs.push_back(
               (*m_sources[i]).fluxDensityDeriv(evt, paramNames[j])
               /srcSum);
         }
      }
   }
}

} // namespace Likelihood
