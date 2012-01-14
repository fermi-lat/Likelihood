/**
 * @file ExpCutoffSEDPeak.cxx
 * @brief Implementation for the ExpCutoff with SED peak energy and flux as variables
 * @author Rolf Buehler
 *
 * $Header$
 */

#ifndef Likelihood_ExpCutoffSEDPeak_h
#define Likelihood_ExpCutoffSEDPeak_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {

  /*
   * @class ExpCutoffSEDPeak
   *
   * @brief Power Law with Super Exponential Cutoff function
   *
   * @author Rolf Buehler
   *
   */
  
  class ExpCutoffSEDPeak : public optimizers::Function {
    
  public:
    
    ExpCutoffSEDPeak(){ init( 10., -2.1, 1000.);}
    ExpCutoffSEDPeak(
			   double Fpeak,
			   double Index, 
			   double Epeak)
      {init(Fpeak, Index, Epeak);}
    
    double value(optimizers::Arg&) const;
    
    double derivByParam(optimizers::Arg &x, const std::string &paramName) const;
    
    virtual optimizers::Function *clone() const {
      return new ExpCutoffSEDPeak(*this);
    }
    
  private:
    
    void init(double Fpeak, double Index, double Epeak);
    
  };
  
} // namespace Likelihood

#endif // Likelihood_ExpCutoffSEDPeak_h
