#ifndef Statistic_h
#define Statistic_h

#include "Function.h"

/** 
 * @class Statistic
 *
 * @brief Objective functions used for parameter estimation.  
 * The standard example is the negative log-likelihood.
 * 
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/Statistic.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class Statistic : Function {
    
public:
    
    Statistic();
    ~Statistic();
    
private:
        
};

#endif // Statistic_h
