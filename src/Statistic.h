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
 * $Header$
 */

class Statistic : Function {
    
public:
    
    Statistic();
    ~Statistic();
    
private:
        
};

#endif // Statistic_h
