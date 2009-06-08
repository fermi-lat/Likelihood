/*
 * @file SmoothBrokenPowerLaw.h
 * @brief Declaration for the SmoothBrokenPowerLaw Function class
 * @author Benoit Lott
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SmoothBrokenPowerLaw.h,v 1.2 2007/03/20 23:46:22 jchiang Exp $
 */

#ifndef Likelihood_SmoothBrokenPowerLaw_h
#define Likelihood_SmoothBrokenPowerLaw_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class SmoothBrokenPowerLaw
   *
   * @brief Power Law with Super Exponential Cutoff function
   *
   * @author Benoit Lott
   *
   */
  
  class SmoothBrokenPowerLaw : public optimizers::Function {
    
  public:
    
    SmoothBrokenPowerLaw(){ init( 10., -2.1,100, -2.1, 1000., 0.2 );}
    SmoothBrokenPowerLaw(
			   double Prefactor,
			   double Index1, 
			   double Scale,
			   double Index2,
			   double BreakValue,
			   double Beta)
    {init(Prefactor, Index1, Scale, Index2, BreakValue, Beta);}
    
    double value(optimizers::Arg&) const;
    
    double derivByParam(optimizers::Arg &x, const std::string &paramName) const;
    
    virtual optimizers::Function *clone() const {
      return new SmoothBrokenPowerLaw(*this);
    }
    
  private:
    
    void init(double Prefactor,  
	      double Index1, 
	      double Scale,
	      double Index2,
	      double BreakValue,
	      double Beta);
  };
  
} // namespace Likelihood

#endif // Likelihood_SmoothBrokenPowerLaw_h
