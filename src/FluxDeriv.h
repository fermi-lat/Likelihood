/**
 * @file FluxDeriv.h
 * @brief Functor class that wraps a Function to provide an interface
 * to that function's partial derivative wrt a named parameter.
 * 
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * 
 * $Header$
 */

#ifndef Likelihood_FluxDeriv
#define Likelihood_FluxDeriv

#include "optimizers/Function.h"

namespace Likelihood {

class FluxDeriv : public optimizers::Function {

public:

   FluxDeriv(const optimizers::Function & func, const std::string & parName) 
      : m_func(func), m_parName(parName) {}

   virtual double value(optimizers::Arg & x) const {
      return m_func.derivByParam(x, m_parName);
   }
   
   virtual double derivByParam(optimizers::Arg &, const std::string &) const {
      throw std::runtime_error("FluxDeriv::deriveByParam not implemented");
   }

protected:

   virtual Function * clone() const {
      return 0;
   }

private:

   const optimizers::Function & m_func;

   std::string m_parName;

};

} // namespace Likelihood

#endif // Likelihood_FluxDeriv
