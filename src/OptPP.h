/** 
 * @file OptPP.h
 * @brief OptPP declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/OptPP.h,v 1.4 2003/05/02 19:02:15 jchiang Exp $
 */

#ifndef Likelihood_OptPP_h
#define Likelihood_OptPP_h

#include "Likelihood/Optimizer.h"
#include "Likelihood/Statistic.h"

#ifdef HAVE_OPT_PP
#include "Opt.h"
#endif

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/OptPP.h,v 1.4 2003/05/02 19:02:15 jchiang Exp $
 */

class OptPP : public Optimizer {
    
public:
    
   OptPP(Statistic &stat) {s_stat = &stat;}
   virtual ~OptPP() {}

   void find_min(int verbose = 0, double tol = 1e-5);
    
protected:

   static int s_verbose;

#ifdef HAVE_OPT_PP
   //! interface to the objective function that OPT++ expects
   static void statInterface(int mode, int ndim, const ColumnVector &x,
                             double &fx, ColumnVector &gx, int &result);

   //! returns the initial parameter values to the OPT++ routines
   static void statInit(int ndim, ColumnVector &x);

   //! do-nothing helper function for use with OptBCQNewton
   static void update_model(int, int, ColumnVector) {}
#endif

   static Statistic *s_stat;

};

} // namespace Likelihood

#endif // Likelihood_OptPP_h
