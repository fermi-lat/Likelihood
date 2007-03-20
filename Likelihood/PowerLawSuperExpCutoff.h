/*
 * @file PowerLawSuperExpCutoff.h
 * @brief Declaration for the PowerLawSuperExpCutoff Function class
 * @author Damien Parent
 *
 * $Header$
 */

#ifndef Likelihood_PowerLawSuperExpCutoff_h
#define Likelihood_PowerLawSuperExpCutoff_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class PowerLawSuperExpCutoff
   *
   * @brief Power Law with Super Exponential Cutoff function
   *
   * @author Damien Parent
   *
   */
  
  class PowerLawSuperExpCutoff : public optimizers::Function {
    
  public:
    
    PowerLawSuperExpCutoff(){ init( 10., -2.1, 1000., 10000., 2. );}
    PowerLawSuperExpCutoff(
			   double Prefactor,
			   double Index1, 
			   double Scale,
			   double Cutoff,
			   double Index2)
      {init(Prefactor, Index1, Scale, Cutoff, Index2);}
    
    double value(optimizers::Arg&) const;
    
    double derivByParam(optimizers::Arg &x, const std::string &paramName) const;
    
    virtual optimizers::Function *clone() const {
      return new PowerLawSuperExpCutoff(*this);
    }
    
  private:
    
    void init(double Prefactor, double Index1, double Scale, double Cutoff, double Index2);
    
  };
  
} // namespace Likelihood

#endif // Likelihood_PowerLawSuperExpCutoff_h
