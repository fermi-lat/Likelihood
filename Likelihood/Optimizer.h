#ifndef SciTools_Optimizer_h
#define SciTools_Optimizer_h

#include <vector>
#include <string>

/** 
 * @class Optimizer
 *
 * @brief Base class for Science Tools Optimizers, e.g., variable metric
 * methods, Powell, Nelder-Mead, etc.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/Optimizer.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class Optimizer {
    
public:
    
   virtual Optimizer();
   virtual ~Optimizer();
    
private:
        
};

#endif // SciTools_Optimizer_h
