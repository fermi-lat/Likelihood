/**
 * @file ExpCutoff.h
 * @brief Declaration for the ExpCutoff Function class
 * @author Luis C. Reyes
 *
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
 * @author Luis C. Reyes
 *
 */

class ExpCutoff : public optimizers::Function {

public:

   ExpCutoff(){init(10., -2.1, 100., 10., 150., 0., 0.);}
   ExpCutoff(double Prefactor, double Index, double Scale,
            double Ebreak, double P1, double P2, double P3)
      {init(Prefactor, Index, Scale, Ebreak, P1, P2, P3);}

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &x, const std::string &paramName) const;

   //double integral(Arg &xmin, Arg &xmax) const;

   virtual optimizers::Function *clone() const {
      return new ExpCutoff(*this);
   }

private:

   void init(double Prefactor, double Index, double Scale,
             double Ebreak, double P1, double P2, double P3);

};

} // namespace Likelihood

#endif // Likelihood_ExpCutoff_h
