/** 
 * @file AbsEdge.h
 * @brief AbsEdge class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AbsEdge.h,v 1.6 2003/05/29 20:10:45 jchiang Exp $
 */

#ifndef Likelihood_AbsEdge_h
#define Likelihood_AbsEdge_h

#include "Likelihood/Function.h"

namespace Likelihood {

class Arg;

/** 
 * @class AbsEdge
 *
 * @brief This Function models an absorption edge as a multiplicative 
 * spectral component.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AbsEdge.h,v 1.6 2003/05/29 20:10:45 jchiang Exp $
 */
    
class AbsEdge : public Function {
public:

   AbsEdge() {init(1, 1, -3);}
   AbsEdge(double Tau0, double E0, double Index = -3)
      {init(Tau0, E0, Index);}

   double value(Arg &) const;

   double derivByParam(Arg &, const std::string &paramName) const
      throw(ParameterNotFound);

   virtual Function *clone() const {
      return new AbsEdge(*this);
   }

private:

   // disable this
   double integral(Arg &, Arg &) const {return 0;}

   void init(double Tau0, double E0, double Index);

};

} // namespace Likelihood

#endif // Likelihood_AbsEdge_h

