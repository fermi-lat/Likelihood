/** 
 * @file LogGaussianDeriv.h
 * @brief Declaration for the LogGaussianDeriv Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/LogGaussianDeriv.h,v 1.3 2009/10/19 19:15:04 jchiang Exp $
 */

#ifndef Likelihood_LogGaussianDeriv_h
#define Likelihood_LogGaussianDeriv_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class LogGaussianDeriv
 *
 * @brief A power-law function that uses integrated flux and index
 * as free parameters and upper and lower bounds of integration as 
 * fixed parameters.
 *
 */
    
class LogGaussianDeriv : public optimizers::Function {

public:

   LogGaussianDeriv(double norm=1, double mean=0, double sigma=1);

   double value(optimizers::Arg & arg) const;

   double derivByParam(optimizers::Arg & x, 
                       const std::string & paramName) const;

   double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const;
   
   virtual Function * clone() const {
      return new LogGaussianDeriv(*this);
   }

private:

};

} // namespace Likelihood

#endif // Likelihood_LogGaussianDeriv_h
