/** @file Optimizer.h
 * @brief Declaration of Optimizer base class
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Optimizer.h,v 1.9 2003/05/02 19:02:13 jchiang Exp $
 */

#ifndef Optimizer_h
#define Optimizer_h

#include <iostream>

//#define HAVE_OPT_PP
#define HAVE_OPT_LBFGS

namespace Likelihood {

/** 
 * @class Optimizer
 *
 * @brief Abstract base class for objective function optimizers.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Optimizer.h,v 1.9 2003/05/02 19:02:13 jchiang Exp $
 */

class Optimizer {
    
public:
    
   Optimizer() {}
   virtual ~Optimizer() {}

   virtual void find_min(int verbose, double tol) = 0;
    
protected:

};

} // namespace Likelhood

#endif // Optimizer_h
