/** 
 * @file BrokenPowerLaw.h
 * @brief Declaration for the BrokenPowerLaw Function class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_BrokenPowerLaw_h
#define Likelihood_BrokenPowerLaw_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

/** 
 * @class BrokenPowerLaw
 *
 * @brief A broken power-law function
 *
 * @author J. Chiang
 *    
 * $Header$
 */
    
class BrokenPowerLaw : public optimizers::Function {

public:

   BrokenPowerLaw(){init(0., -2., -2., 1.);}
   BrokenPowerLaw(double Prefactor, double Index1, double Index2, 
                  double BreakValue)
      {init(Prefactor, Index1, Index2, BreakValue);}

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &x, const std::string &paramName) const
      throw(optimizers::ParameterNotFound);

   virtual optimizers::Function *clone() const {
      return new BrokenPowerLaw(*this);
   }

private:

   void init(double Prefactor, double Index1, 
             double Index2, double BreakValue);

// Disable this method.
   double integral(optimizers::Arg &, optimizers::Arg &) const 
      {return 0;}

};

} // namespace Likelihood

#endif // Likelihood_BrokenPowerLaw_h
