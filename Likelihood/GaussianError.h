/** 
 * @file GaussianError.h
 * @brief Declaration for the GaussianError Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/GaussianError.h,v 1.3 2015/03/03 18:05:36 jchiang Exp $
 */

#ifndef Likelihood_GaussianError_h
#define Likelihood_GaussianError_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class GaussianError
 *
 * @brief A power-law function that uses integrated flux and index
 * as free parameters and upper and lower bounds of integration as 
 * fixed parameters.
 *
 */
    
  class GaussianError : public optimizers::Function {
    
  public:
    
    typedef enum { Norm, Mean, Sigma, Offset } ParamTypes;

    GaussianError(double norm=1, double mean=0, double sigma=1, double offset=0);
    
    double derivative(const optimizers::Arg & x) const;
    
    virtual Function * clone() const {
      return new GaussianError(*this);
    }
    
  protected:
    
    double value(const optimizers::Arg & arg) const;
    
    double derivByParamImp(const optimizers::Arg & x, 
			   const std::string & paramName) const;
    
  };
  
} // namespace Likelihood

#endif // Likelihood_GaussianError_h
