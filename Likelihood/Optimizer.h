#ifndef Optimizer_h
#define Optimizer_h

#include "Likelihood/Statistic.h"

namespace Likelihood {

/** 
 * @class Optimizer
 *
 * @brief Abstract base class for objective function optimizers.
 *
 * @author J. Chiang
 *    
 * $Header:
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
