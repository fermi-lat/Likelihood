/** @file logLike_gauss.h
 * @brief Declaration of logLike_gauss class
 * $Header:
 */

#ifndef logLike_gauss_h
#define logLike_gauss_h

#include "Likelihood/Statistic.h"
#include "Likelihood/Arg.h"

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

   logLike_gauss(){
      deleteAllSources();
   };
   virtual ~logLike_gauss(){};

   //! return the objective function value taking the free parameters 
   //! as the function argument
   double value(const std::vector<double> &paramVec);
   
   void getFreeDerivs(std::vector<double> &freeDerivs){};

};

} // namespace Likelihood

#endif // logLike_gauss_h
