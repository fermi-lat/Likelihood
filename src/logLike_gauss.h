/** 
 * @file logLike_gauss.h
 * @brief Declaration of logLike_gauss class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/logLike_gauss.h,v 1.9 2003/06/11 17:08:05 jchiang Exp $
 */

#ifndef Likelihood_logLike_gauss_h
#define Likelihood_logLike_gauss_h

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/logLike_gauss.h,v 1.9 2003/06/11 17:08:05 jchiang Exp $
 */

class logLike_gauss : public Statistic {
    
public:

   logLike_gauss() {deleteAllSources();}

   virtual ~logLike_gauss() {}

   //! return the objective function value taking the free parameters 
   //! as the function argument
   double value(const std::vector<double> &paramVec);
   
   void getFreeDerivs(std::vector<double> &){}

};

} // namespace Likelihood

#endif // Likelihood_logLike_gauss_h
