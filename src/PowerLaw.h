/** @file PowerLaw.h
 * @brief Declaration for the PowerLaw Function class
 * @author J. Chiang
 *
 * $Header$
 */

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw.h,v 1.5 2003/03/16 21:53:26 jchiang Exp $
 */
    
class PowerLaw : public Function {
public:

   PowerLaw(){init(0, -2, 1);}
   PowerLaw(double Prefactor, double Index, double Scale)
      {init(Prefactor, Index, Scale);}

   double value(Arg&) const;

   double derivByParam(Arg &x, const std::string &paramName) const;

   double integral(Arg &xmin, Arg &xmax) const;

private:

   void init(double Prefactor, double Index, double Scale);

};

} // namespace Likelihood

