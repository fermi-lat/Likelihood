/** 
 * @file PowerLaw.h
 * @brief Declaration for the PowerLaw Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw.h,v 1.13 2003/07/19 04:38:03 jchiang Exp $
 */

#ifndef Likelihood_PowerLaw_h
#define Likelihood_PowerLaw_h

#include "optimizers/Function.h"
#include "optimizers/Arg.h"

namespace Likelihood {
/** 
 * @class PowerLaw
 *
 * @brief A power-law function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw.h,v 1.13 2003/07/19 04:38:03 jchiang Exp $
 */
    
class PowerLaw : public optimizers::Function {
public:

   PowerLaw(){init(0, -2, 1);}
   PowerLaw(double Prefactor, double Index, double Scale)
      {init(Prefactor, Index, Scale);}

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &x, const std::string &paramName) const
      throw(optimizers::ParameterNotFound);

   double integral(optimizers::Arg &xmin, optimizers::Arg &xmax) const;

   virtual optimizers::Function *clone() const {
      return new PowerLaw(*this);
   }

private:

   void init(double Prefactor, double Index, double Scale);

};

} // namespace Likelihood

#endif // Likelihood_PowerLaw_h
