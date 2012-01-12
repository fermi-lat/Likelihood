/*
 * @file SmoothDoubleBrokenPowerLaw.h
 * @brief Declaration for the SmoothDoubleBrokenPowerLaw Function class
 * @author Keith Bechtol
 *
 * $Header: /usr/local/CVS/SLAC/Likelihood/Likelihood/SmoothDoubleBrokenPowerLaw.h,v 1.1 2009/06/08 06:05:48 jchiang Exp $
 */

#ifndef Likelihood_SmoothDoubleBrokenPowerLaw_h
#define Likelihood_SmoothDoubleBrokenPowerLaw_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class SmoothDoubleBrokenPowerLaw
   *
   * @brief Power Law with two smooth breaks each corresponding to a steepening of the spectral index
   *
   * @author Keith Bechtol
   *
   */
  
  class SmoothDoubleBrokenPowerLaw : public optimizers::Function {
    
  public:
    
    SmoothDoubleBrokenPowerLaw(){ init( 10., -1.5, 100, -2.0, 1000., 0.1, -2.5, 10000., 0.1 );}
    SmoothDoubleBrokenPowerLaw(
			   double Prefactor,
			   double Index1, 
			   double Scale,
			   double Index2,
			   double BreakValue12,
			   double Beta12,
			   double Index3,
			   double BreakValue23,
			   double Beta23
			   )
    {init(Prefactor, Index1, Scale, Index2, BreakValue12, Beta12, Index3, BreakValue23, Beta23);}
    
    double value(optimizers::Arg&) const;
    
    double derivByParam(optimizers::Arg &x, const std::string &paramName) const;
    
    virtual optimizers::Function *clone() const {
      return new SmoothDoubleBrokenPowerLaw(*this);
    }
    
  private:
    
    void init(double Prefactor,  
	      double Index1, 
	      double Scale,
	      double Index2,
	      double BreakValue12,
	      double Beta12,
	      double Index3,
	      double BreakValue23,
	      double Beta23);
  };
  
} // namespace Likelihood

#endif // Likelihood_SmoothDoubleBrokenPowerLaw_h
