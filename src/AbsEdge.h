/** @file AbsEdge.h
 * @brief AbsEdge class declaration
 * @author J. Chiang
 *
 * $Header$
 */

#include "Likelihood/Function.h"
#include "Likelihood/Arg.h"

namespace Likelihood {
/** 
 * @class AbsEdge
 *
 * @brief This Function models an absorption edge as a multiplicative 
 * spectral component.
 *
 * @author J. Chiang
 *    
 * $Header$
 */
    
class AbsEdge : public Function {
public:

   AbsEdge(){init(1, 1, -3);}
   AbsEdge(double Tau0, double E0, double Index = -3)
      {init(Tau0, E0, Index);}

   double value(Arg &) const;

   double derivByParam(Arg &, const std::string &paramName) const;

private:

   // disable this
   double integral(Arg &, Arg &) const {return 0;}

   void init(double Tau0, double E0, double Index);

};

} // namespace Likelihood

