/** 
 * @file SumFunction.h
 * @brief Declaration of SumFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SumFunction.h,v 1.6 2003/03/22 01:22:50 jchiang Exp $
 */

#ifndef SumFunction_h
#define SumFunction_h

#include "Likelihood/CompositeFunction.h"

namespace Likelihood {
/** 
 * @class SumFunction
 *
 * @brief A Function that returns the linear sum of two Functions
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SumFunction.h,v 1.6 2003/03/22 01:22:50 jchiang Exp $
 *
 */
    
class SumFunction : public CompositeFunction {
public:

   SumFunction(Function &a, Function &b);

   double value(Arg &x) const
      {return m_a->value(x) + m_b->value(x);}

   double integral(Arg &xmin, Arg &xmax) const
      {return m_a->integral(xmin, xmax) + m_b->integral(xmin, xmax);}

   virtual Function* clone() const {
      return new SumFunction(*this);
   }

protected:

   void fetchDerivs(Arg &x, std::vector<double> &derivs, bool getFree) const;

};

} // namespace Likelihood

#endif // SumFunction_h
