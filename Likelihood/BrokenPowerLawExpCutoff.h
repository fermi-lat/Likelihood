/**
 * @file BrokenPowerLawExpCutoff.h
 * @brief Declaration for the BrokenPowerLawExpCutoff Function class
 * @author Jennifer Carson
 *
 * $Header$
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
 * @author Jennifer Carson
 *
 */

class BrokenPowerLawExpCutoff : public optimizers::Function {

public:

   BrokenPowerLawExpCutoff(){init(10., -2.1, -2.1, 1000., 10., 150.);}
   BrokenPowerLawExpCutoff(double Prefactor, double Index1, 
			   double Index2, double BreakValue,
            double Eabs, double P1)
      {init(Prefactor, Index1, Index2, BreakValue, Eabs, P1);}

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &x, const std::string &paramName) const;

   //double integral(Arg &xmin, Arg &xmax) const;

   virtual optimizers::Function *clone() const {
      return new BrokenPowerLawExpCutoff(*this);
   }

private:

   void init(double Prefactor, double Index1, double Index2, 
	     double BreakValue, double Eabs, double P1);

};

} // namespace Likelihood

#endif // Likelihood_BrokenPowerLawExpCutoff_h
