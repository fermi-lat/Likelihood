/** 
 * @file PowerLaw2.h
 * @brief Declaration for the PowerLaw2 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw2.h,v 1.2 2004/12/22 06:00:40 jchiang Exp $
 */

#ifndef Likelihood_PowerLaw2_h
#define Likelihood_PowerLaw2_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class PowerLaw2
 *
 * @brief A power-law function that uses integrated flux and index
 * as free parameters and upper and lower bounds of integration as 
 * fixed parameters.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw2.h,v 1.2 2004/12/22 06:00:40 jchiang Exp $
 */
    
class PowerLaw2 : public optimizers::Function {

public:

   PowerLaw2(double Integral=1., double Index=-2., 
             double LowerLimit=100., double UpperLimit=2e5);

   double value(optimizers::Arg & x) const;

   double derivByParam(optimizers::Arg & x, 
                       const std::string & paramName) const;

   double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const;

   virtual Function * clone() const {
      return new PowerLaw2(*this);
   }

};

} // namespace Likelihood

#endif // Likelihood_PowerLaw2_h
