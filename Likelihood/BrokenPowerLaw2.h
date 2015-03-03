/** 
 * @file BrokenPowerLaw2.h
 * @brief Declaration for the BrokenPowerLaw2 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BrokenPowerLaw2.h,v 1.2 2012/07/31 19:38:34 jchiang Exp $
 */

#ifndef Likelihood_BrokenPowerLaw2_h
#define Likelihood_BrokenPowerLaw2_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class BrokenPowerLaw2
 *
 * @brief A broken power-law function that uses integrated flux,
 * indices, and break value as free parameters and upper and lower
 * bounds of integration as fixed parameters.
 *
 */
    
class BrokenPowerLaw2 : public optimizers::Function {

public:

   BrokenPowerLaw2(double Integral=1., double Index1=-2., 
                   double Index2=-3., double BreakValue=1000.,
                   double LowerLimit=20., double UpperLimit=2e5);

   double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const;

   virtual Function * clone() const {
      return new BrokenPowerLaw2(*this);
   }

protected:

   double value(optimizers::Arg & x) const;

   double derivByParamImp(optimizers::Arg & x, 
                          const std::string & paramName) const;

private:

   double N0_value() const;

};

} // namespace Likelihood

#endif // Likelihood_BrokenPowerLaw2_h
