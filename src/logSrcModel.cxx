/** @file logSrcModel.cxx
 * @brief logSrcModel class implementation
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include "Likelihood/logSrcModel.h"

namespace Likelihood {

double logSrcModel::value(Arg &xarg) const {
   Event evt;
   dynamic_cast<EventArg &>(xarg).fetchValue(evt);

   double my_value = 0;
   for (unsigned int i = 0; i < getNumSrcs(); i++) {
      my_value += s_sources[i]->fluxDensity(evt);
   }
   if (my_value > 0) {
      return log(my_value);
   } else {
//      std::cerr << "logSrcModel::value: my_value < 0" << std::endl;
      return 0;
   }
}

void logSrcModel::fetchDerivs(Arg &xarg, std::vector<double> &derivs, 
                              bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   Event evt;
   dynamic_cast<EventArg &>(xarg).fetchValue(evt);
   double srcSum = exp(value(xarg));

   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
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
               (*s_sources[i]).fluxDensityDeriv(evt, paramNames[j])/srcSum
               );
         }
      }
   }
}

} // namespace Likelihood
