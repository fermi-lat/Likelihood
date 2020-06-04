/*
 * @file PowerLawSuperExpCutoff4.h
 * @brief Declaration for the PowerLawSuperExpCutoff4 Function class
 * @author Jean Ballet, Philippe Bruel
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PowerLawSuperExpCutoff4.h,v 1.1 2017/06/17 00:02:39 mdwood Exp $
 */

#ifndef Likelihood_PowerLawSuperExpCutoff4_h
#define Likelihood_PowerLawSuperExpCutoff4_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class PowerLawSuperExpCutoff4
   *
   * @brief Power Law with Super Exponential Cutoff function
   *
   * @author Jean Ballet, Philippe Bruel
   *
   */
  
  class PowerLawSuperExpCutoff4 : public optimizers::Function {
    
  public:
    
    PowerLawSuperExpCutoff4(double Prefactor=1.,
			    double IndexS=-1.5,
			    double Scale=1000.,
			    double ExpfactorS=0.4,
			    double Index2=0.5);
    
    virtual optimizers::Function * clone() const {
      return new PowerLawSuperExpCutoff4(*this);
    }
    
  protected:

    double value(const optimizers::Arg &) const;
    
    double derivByParamImp(const optimizers::Arg & x,
                           const std::string &paramName) const;
    
  };
  
} // namespace Likelihood

#endif // Likelihood_PowerLawSuperExpCutoff4_h
