/** @file logLike_gauss.h
 * @brief Declaration of logLike_gauss class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/logLike_gauss.h,v 1.7 2003/03/22 01:22:51 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/logLike_gauss.h,v 1.7 2003/03/22 01:22:51 jchiang Exp $
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

#endif // logLike_gauss_h
