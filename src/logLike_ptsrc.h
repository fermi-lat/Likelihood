/** @file logLike_ptsrc.h
 * @brief Declaration of logLike_ptsrc class
 * $Header:
 */

#ifndef logLike_ptsrc_h
#define logLike_ptsrc_h

#include "../Likelihood/Statistic.h"

namespace Likelihood {

/** 
 * @class logLike_ptsrc
 *
 * @brief Objective function for the log(likelihood) of a 
 * single point source.
 *
 * @author J. Chiang
 *    
 * $Header: 
 */

class logLike_ptsrc : public Statistic {
    
public:

   logLike_ptsrc(){};
   virtual ~logLike_ptsrc(){};

   //! return the objective function value taking the free parameters 
   //! as the function argument
   double value(const std::vector<double> &paramVec);
   double operator()(const std::vector<double> &paramVec) 
      {return value(paramVec);};

};

} // namespace Likelihood

#endif // logLike_ptsrc_h
