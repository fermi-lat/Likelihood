/** @file OptPP.cxx
 * @brief Implementation OptPP class that provides an interface to the
 * OPT++ package
 * @author J. Chiang
 *
 * $Header$
 */

#include <cassert>
#include <fstream>

#ifdef HAVE_OPTIMIZERS

#include "NLF.h"
#include "CompoundConstraint.h"
#include "BoundConstraint.h"
#include "OptBCQNewton.h"
//#include "OptBaQNewton.h"
//#include "OptBCEllipsoid.h"

#endif //HAVE_OPTIMIZERS

#include "Likelihood/Parameter.h"
#include "OptPP.h"

namespace Likelihood {

Statistic *OptPP::s_stat = 0;
int OptPP::s_verbose = 0;

#ifdef HAVE_OPTIMIZERS

void OptPP::statInterface(int mode, int ndim, const ColumnVector &x,
                          double &fx, ColumnVector &gx, int &result) {
   
   assert(ndim == static_cast<int>(s_stat->getNumFreeParams()));

   if (mode & NLPFunction) {
      std::vector<double> paramValues;
      if (s_verbose) std::cout << "parameter values: ";
      for (int i = 0; i < ndim; i++) {
         double param_val = x(i+1);
         paramValues.push_back(param_val);
         if (s_verbose) std::cout << param_val << "  ";
      }

      fx = -s_stat->value(paramValues);
      result = NLPFunction;

      if (s_verbose) std::cout << "f(x) = " << fx << std::endl;
   }

   if (mode & NLPGradient) {
      std::vector<double> derivsVec;
      s_stat->getFreeDerivs(derivsVec);
      if (s_verbose) std::cout << "gradients: ";
      for (int i = 0; i < ndim; i++) {
         gx(i+1) = -derivsVec[i];
         if (s_verbose) std::cout << gx(i+1) << "  ";
      }
      if (s_verbose) std::cout << std::endl;

      result = NLPGradient;
   }
}

void OptPP::statInit(int ndim, ColumnVector &x) {

   assert(ndim == static_cast<int>(s_stat->getNumFreeParams()));

   std::vector<double> paramValues;
   s_stat->getFreeParamValues(paramValues);

   for (int i = 0; i < ndim; i++)
      x(i+1) = paramValues[i];
}

void OptPP::find_min(int verbose, double tol) {

   s_verbose = verbose;

   int ndim = s_stat->getNumFreeParams();
   
// retrieve the current set of Parameters

   std::vector<Parameter> params;
   s_stat->getFreeParams(params);

// set the lower and upper bounds constraints

   ColumnVector lower(ndim), upper(ndim);
   for (int i = 0; i < ndim; i++) {
      std::pair<double, double> bounds = params[i].getBounds();
      lower(i+1) = bounds.first;
      upper(i+1) = bounds.second;
//        if (s_verbose) std::cout << lower(i+1) << "  "
//                                 << upper(i+1) << std::endl;
   }

   Constraint myBounds = new BoundConstraint(ndim, lower, upper);
   CompoundConstraint *myConstraints = new CompoundConstraint(myBounds);

// instantiate a non-linear function object with 1st derivatives
// using the appropriate static functions and constraints

   NLF1 nlp(ndim, OptPP::statInterface, OptPP::statInit, myConstraints);
   nlp.initFcn();

// Create an objective function for bound constrained quasi-Newton
// optimization, using a LineSearch strategy

   OptBCQNewton myObjFunc(&nlp, update_model);
   myObjFunc.setSearchStrategy(LineSearch);

//   OptBaQNewton myObjFunc(&nlp, update_model);
//   OptBCEllipsoid myObjFunc(&nlp, update_model);

   if (s_verbose) myObjFunc.setOutputFile(std::cout);

// do the business....

   myObjFunc.optimize();
}
#else

void OptPP::find_min(int verbose, double tol) {
// do nothing
}

#endif // HAVE_OPTIMIZERS


} // namespace Likelihood
