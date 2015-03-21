/**
 * @file BrokenPowerLawExpCutoff.h
 * @brief Declaration for the BrokenPowerLawExpCutoff Function class
 * @author Jennifer Carson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BrokenPowerLawExpCutoff.h,v 1.2 2015/03/03 18:05:36 jchiang Exp $
 */

#ifndef Likelihood_BrokenPowerLawExpCutoff_h
#define Likelihood_BrokenPowerLawExpCutoff_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

/**
 * @class BrokenPowerLawExpCutoff
 *
 * @brief Broken Power Law with Exponential Cutoff function
 *
 */

class BrokenPowerLawExpCutoff : public optimizers::Function {

public:

   BrokenPowerLawExpCutoff(double Prefactor=10 , double Index1=-2.1, 
			   double Index2=-2.1, double BreakValue=1000.,
                           double Eabs=10., double P1=150.);

   virtual optimizers::Function * clone() const {
      return new BrokenPowerLawExpCutoff(*this);
   }

protected:

   double value(const optimizers::Arg&) const;

   double derivByParamImp(const optimizers::Arg & x,
                          const std::string & paramName) const;
};

} // namespace Likelihood

#endif // Likelihood_BrokenPowerLawExpCutoff_h
