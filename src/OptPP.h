/** @file OptPP.h
 * @brief Declaration of OptPP class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef OptPP_h
#define OptPP_h

#include "Likelihood/Optimizer.h"
#include "Likelihood/Statistic.h"

#ifdef HAVE_OPTIMIZER
#include "Opt.h"
#endif //HAVE_OPTIMIZER

namespace Likelihood {

/** 
 * @class OptPP
 *
 * @brief Wrapper class for the OPT++ package
 * (http://csmr.ca.sandia.gov/projects/opt/).  Presently, we use their
 * bound constrained quasi-Newton optimizer OptBCQNewton with a
 * LineSearch strategy.  A more general interface that allows other
 * optimizer methods to be chosen will be implemented...eventually.
 *
 * @author J. Chiang
 *    
 * $Header$
 *
 */

class OptPP : public Optimizer {
    
public:
    
   OptPP(Statistic *stat) {s_stat = stat;}
   virtual ~OptPP() {}

   void find_min(int verbose = 0, double tol = 1e-5);
    
protected:

   static int s_verbose;

#ifdef HAVE_OPTIMIZERS

   //! interface to the objective function that OPT++ expects
   static void statInterface(int mode, int ndim, const ColumnVector &x,
                             double &fx, ColumnVector &gx, int &result);

   //! returns the initial parameter values to the OPT++ routines
   static void statInit(int ndim, ColumnVector &x);

   //! do-nothing helper function for use with OptBCQNewton
   static void update_model(int, int, ColumnVector) {}

#endif //HAVE_OPTIMIZERS

   static Statistic *s_stat;

};

} // namespace Likelihood

#endif // OptPP_h
