/** 
 * @file AbsEdge.h
 * @brief AbsEdge class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AbsEdge.h,v 1.7 2003/07/19 04:38:02 jchiang Exp $
 */

#ifndef Likelihood_AbsEdge_h
#define Likelihood_AbsEdge_h

#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class AbsEdge
 *
 * @brief This Function models an absorption edge as a multiplicative 
 * spectral component.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AbsEdge.h,v 1.7 2003/07/19 04:38:02 jchiang Exp $
 */
    
class AbsEdge : public optimizers::Function {
public:

   AbsEdge() {init(1, 1, -3);}
   AbsEdge(double Tau0, double E0, double Index = -3)
      {init(Tau0, E0, Index);}

   double value(optimizers::Arg &) const;

   double derivByParam(optimizers::Arg &, const std::string &paramName) const
      throw(optimizers::ParameterNotFound);

   virtual optimizers::Function *clone() const {
      return new AbsEdge(*this);
   }

private:

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

   void init(double Tau0, double E0, double Index);

};

} // namespace Likelihood

#endif // Likelihood_AbsEdge_h

