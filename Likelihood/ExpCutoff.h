/**
 * @file ExpCutoff.h
 * @brief Declaration for the ExpCutoff Function class
 * @author Luis C. Reyes
 *
 * $Header: $
 */

#ifndef Likelihood_ExpCutoff_h
#define Likelihood_ExpCutoff_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

/**
 * @class ExpCutoff
 *
 * @brief Power Law with Exponential Cutoff function
 *
 */

class ExpCutoff : public optimizers::Function {

public:

   ExpCutoff(double Prefactor=10., double Index=-2.1, double Scale=100.,
             double Ebreak=10., double P1=150., double P2=0, double P3=0);

   virtual optimizers::Function * clone() const {
      return new ExpCutoff(*this);
   }

protected:

   double value(optimizers::Arg &) const;

   double derivByParamImp(optimizers::Arg & x, 
                          const std::string & paramName) const;

};

} // namespace Likelihood

#endif // Likelihood_ExpCutoff_h
