/**
 * @file ExpCutoffSEDPeak.cxx
 * @brief Implementation for the ExpCutoff with SED peak energy and flux as variables
 * @author Rolf Buehler
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExpCutoffSEDPeak.h,v 1.2 2015/03/03 18:05:36 jchiang Exp $
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
   */
  
  class ExpCutoffSEDPeak : public optimizers::Function {
    
  public:
    
    ExpCutoffSEDPeak(double Fpeak=10.,
                     double Index=-2.1, 
                     double Epeak=1000.);
    
    virtual optimizers::Function * clone() const {
      return new ExpCutoffSEDPeak(*this);
    }
    
  protected:

    double value(const optimizers::Arg &) const;
    
    double derivByParamImp(const optimizers::Arg & x,
                           const std::string & paramName) const;
    
  };
  
} // namespace Likelihood

#endif // Likelihood_ExpCutoffSEDPeak_h
