/** 
 * @file BrokenPowerLaw2.h
 * @brief Declaration for the BrokenPowerLaw2 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BrokenPowerLaw2.h,v 1.3 2015/03/03 18:05:36 jchiang Exp $
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

   double integral(const optimizers::Arg & xmin,
                   const optimizers::Arg & xmax) const;

   virtual Function * clone() const {
      return new BrokenPowerLaw2(*this);
   }

protected:

   double value(const optimizers::Arg & x) const;

   double derivByParamImp(const optimizers::Arg & x, 
                          const std::string & paramName) const;

private:

   double N0_value() const;

};

} // namespace Likelihood

#endif // Likelihood_BrokenPowerLaw2_h
