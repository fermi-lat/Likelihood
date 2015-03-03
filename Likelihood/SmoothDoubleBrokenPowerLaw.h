/*
 * @file SmoothDoubleBrokenPowerLaw.h
 * @brief Declaration for the SmoothDoubleBrokenPowerLaw Function class
 * @author Keith Bechtol
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SmoothDoubleBrokenPowerLaw.h,v 1.1 2012/01/12 16:46:56 jchiang Exp $
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

    SmoothDoubleBrokenPowerLaw(double Prefactor=10.,
                               double Index1=-1.5, 
                               double Scale=100.,
                               double Index2=-2.0,
                               double BreakValue12=1000.,
                               double Beta12=0.1,
                               double Index3=-2.5,
                               double BreakValue23=10000.,
                               double Beta23=0.1);

    virtual optimizers::Function * clone() const {
      return new SmoothDoubleBrokenPowerLaw(*this);
    }
    
  protected:

    double value(optimizers::Arg &) const;
    
    double derivByParamImp(optimizers::Arg & x,
                           const std::string & paramName) const;

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
