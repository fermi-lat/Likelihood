/** 
 * @file LogNormal.h
 * @brief LogNormal class for spectral modeling
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LogNormal.h,v 1.2 2015/03/03 18:05:36 jchiang Exp $
 */

#ifndef Likelihood_LogNormal_h
#define Likelihood_LogNormal_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

/** 
 * @class LogNormal
 *
 */
    
class LogNormal : public optimizers::Function {

public:

   LogNormal(double prefactor=1, double log10_mean=3,
             double log10_sigma=2);

   virtual Function * clone() const {
      return new LogNormal(*this);
   }

protected:

   double value(const optimizers::Arg &) const;

   double derivByParamImp(const optimizers::Arg & x, 
                          const std::string & paramName) const;

};

} // namespace Likelihood

#endif // Likelihood_LogNormal_h
