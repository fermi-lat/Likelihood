/** 
 * @file BrokenPowerLaw3.h
 * @brief Declaration for the BrokenPowerLaw3 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/BrokenPowerLaw3.h,v 1.1 2005/09/19 00:19:16 jchiang Exp $
 */

#ifndef Likelihood_BrokenPowerLaw3_h
#define Likelihood_BrokenPowerLaw3_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class BrokenPowerLaw3
 *
 * @brief A broken power-law function that uses integrated flux,
 * indices, and break value as free parameters and upper and lower
 * bounds of integration as fixed parameters.
 *
 */
    
class BrokenPowerLaw3 : public optimizers::Function {

public:

   BrokenPowerLaw3(double Integral1=1., double Index1=-2., 
                   double Integral2=1., double Index2=-3., 
                   double LowerLimit1=100., double UpperLimit1=1e4,
                   double LowerLimit2=2e4, double UpperLimit2=1e5);

   double value(optimizers::Arg & x) const;

   double derivByParam(optimizers::Arg & x, 
                       const std::string & paramName) const;

   double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const;

   virtual Function * clone() const {
      return new BrokenPowerLaw3(*this);
   }

private:

   double x0_value() const;
   double N0_value(double x0) const;

};

} // namespace Likelihood

#endif // Likelihood_BrokenPowerLaw3_h
