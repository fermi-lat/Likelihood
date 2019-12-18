/*
 * @file PowerLawSuperExpCutoff2.h
 * @brief Declaration for the PowerLawSuperExpCutoff2 Function class
 * @author Matthew Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PowerLawSuperExpCutoff2.h,v 1.1 2017/06/17 00:02:39 mdwood Exp $
 */

#ifndef Likelihood_PowerLawSuperExpCutoff2_h
#define Likelihood_PowerLawSuperExpCutoff2_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class PowerLawSuperExpCutoff2
   *
   * @brief Power Law with Super Exponential Cutoff function
   *
   * @author Matthew Wood
   *
   */
  
  class PowerLawSuperExpCutoff2 : public optimizers::Function {
    
  public:
    
    PowerLawSuperExpCutoff2(double Prefactor=10.,
			    double Index1=-2.1, 
			    double Scale=1000.,
			    double Expfactor=0,
			    double Index2=2.);
    
    virtual optimizers::Function * clone() const {
      return new PowerLawSuperExpCutoff2(*this);
    }
    
  protected:

    double value(const optimizers::Arg &) const;
    
    double derivByParamImp(const optimizers::Arg & x,
                           const std::string &paramName) const;
    
  };
  
} // namespace Likelihood

#endif // Likelihood_PowerLawSuperExpCutoff2_h
