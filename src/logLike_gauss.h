/** @file logLike_gauss.h
 * @brief Declaration of logLike_gauss class
 * $Header:
 */

#ifndef logLike_gauss_h
#define logLike_gauss_h

#include "../Likelihood/Statistic.h"

namespace Likelihood {

/** 
 * @class logLike_gauss
 *
 * @brief Objective function for the log(likelihood) of a sampled
 * 1D Gaussian function.
 *
 * @author J. Chiang
 *    
 * $Header: 
 */

class logLike_gauss : public Statistic {
    
public:

   logLike_gauss(){};
   virtual ~logLike_gauss(){};

   //! return the objective function value taking the free parameters 
   //! as the function argument
   double value(const std::vector<double> &paramVec);
   double operator()(const std::vector<double> &paramVec) 
      {return value(paramVec);};

};

} // namespace Likelihood

#endif // logLike_gauss_h
