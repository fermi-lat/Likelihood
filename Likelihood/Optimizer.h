/** 
 * @file Optimizer.h
 * @brief Declaration of Optimizer base class
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Optimizer.h,v 1.12 2003/06/10 19:19:55 jchiang Exp $
 */

#ifndef Likelihood_Optimizer_h
#define Likelihood_Optimizer_h

#include <iostream>

//#define HAVE_OPT_PP

namespace Likelihood {

/** 
 * @class Optimizer
 *
 * @brief Abstract base class for objective function optimizers.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Optimizer.h,v 1.12 2003/06/10 19:19:55 jchiang Exp $
 */

class Optimizer {
    
public:
    
   Optimizer() {}
   virtual ~Optimizer() {}

   virtual void find_min(int verbose, double tol) = 0;
    
protected:

};

} // namespace Likelhood

#endif // Likelihood_Optimizer_h
