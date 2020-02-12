/*
 * @file PowerLawSuperExpCutoff3.h
 * @brief Declaration for the PowerLawSuperExpCutoff3 Function class
 * @author Matthew Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PowerLawSuperExpCutoff3.h,v 1.1 2017/06/17 00:02:39 mdwood Exp $
 */

#ifndef Likelihood_PowerLawSuperExpCutoff3_h
#define Likelihood_PowerLawSuperExpCutoff3_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class PowerLawSuperExpCutoff3
   *
   * @brief Power Law with Super Exponential Cutoff function
   *
   * @author Jean Ballet, Philippe Bruel
   *
   */

  class PowerLawSuperExpCutoff3 : public optimizers::Function {
    
  public:
    
    PowerLawSuperExpCutoff3(double Prefactor=1.,
			    double IndexS=-2.0, 
			    double Scale=1000.,
			    double Expfactor2=1,
			    double Index2=0.5);
    
    virtual optimizers::Function * clone() const {
      return new PowerLawSuperExpCutoff3(*this);
    }
    
  protected:

    double value(const optimizers::Arg &) const;
    
    double derivByParamImp(const optimizers::Arg & x,
                           const std::string &paramName) const;
    
  };
  
} // namespace Likelihood

#endif // Likelihood_PowerLawSuperExpCutoff3_h
