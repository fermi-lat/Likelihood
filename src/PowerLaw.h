/** 
 * @file PowerLaw.h
 * @brief Declaration for the PowerLaw Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw.h,v 1.11 2003/04/25 18:32:19 jchiang Exp $
 */

#ifndef PowerLaw_h
#define PowerLaw_h

#include "Likelihood/Function.h"
#include "Likelihood/Arg.h"

namespace Likelihood {
/** 
 * @class PowerLaw
 *
 * @brief A power-law function
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw.h,v 1.11 2003/04/25 18:32:19 jchiang Exp $
 */
    
class PowerLaw : public Function {
public:

   PowerLaw(){init(0, -2, 1);}
   PowerLaw(double Prefactor, double Index, double Scale)
      {init(Prefactor, Index, Scale);}

   double value(Arg&) const;

   double derivByParam(Arg &x, const std::string &paramName) const
      throw(ParameterNotFound);

   double integral(Arg &xmin, Arg &xmax) const;

   virtual Function *clone() const {
      return new PowerLaw(*this);
   }

private:

   void init(double Prefactor, double Index, double Scale);

};

} // namespace Likelihood

#endif // PowerLaw_h
