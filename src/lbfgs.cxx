/** @file lbfgs.cxx
 * @brief lbfgs class implementation
 *
 * $Header:
 */

#include <vector>
#include <utility>

//  extern "C" int lbfgs_min__(long *n, double *x, double *xmin, double *xmax, 
//                             long (*func)(long *, double *, double *, double *), 
//                             long *iprint, long *ncall);

#include "Likelihood/Parameter.h"
#include "lbfgs.h"

namespace Likelihood {

Statistic *lbfgs::s_stat = 0;
int lbfgs::s_verbose = 0;

long lbfgs::statInterface(long *nparams, double *param_vals, 
                          double *func_val, double *derivs) {
   
   std::vector<double>paramValues;
   for (int i = 0; i < *nparams; i++) 
      paramValues.push_back(param_vals[i]);

   *func_val = -s_stat->value(paramValues);

   std::vector<double> derivsVec;
   s_stat->getFreeDerivs(derivsVec);

   for (int i = 0; i < *nparams; i++) {
      derivs[i] = -derivsVec[i];
   }

   if (s_verbose) {
      std::cout << *func_val << "  ";
      for (int i = 0; i < *nparams; i++)
         std::cout << param_vals[i] << "  "
                   << derivs[i] << " : ";
      std::cout << std::endl;
   }
   return 0;
}

void lbfgs::find_min(int verbose, double tol) {

   s_verbose = verbose;

   std::vector<Parameter> Params;
   s_stat->getFreeParams(Params);

   long nparams = Params.size();
   double *param_vals = new double[nparams];
   double *param_mins = new double[nparams];
   double *param_maxs = new double[nparams];
   long iprint = verbose - 2;
   long nullval = 0;

   for (int i = 0; i < nparams; i++) {
      param_vals[i] = Params[i].getValue();
      param_mins[i] = Params[i].getBounds().first;
      param_maxs[i] = Params[i].getBounds().second;
   }

//     lbfgs_min__(&nparams, param_vals, param_mins, param_maxs,
//                 &lbfgs::statInterface, &iprint, &nullval);

   delete [] param_vals;
   delete [] param_mins;
   delete [] param_maxs;
}

} // namespace Likelihood
