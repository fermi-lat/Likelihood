/** 
 * @file ProductFunction.h
 * @brief Declaration of ProductFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ProductFunction.h,v 1.6 2003/03/22 01:22:50 jchiang Exp $
 */

#ifndef ProductFunction_h
#define ProductFunction_h

#include "Likelihood/CompositeFunction.h"

namespace Likelihood {
/** 
 * @class ProductFunction
 *
 * @brief A Function that returns the product of two Functions
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ProductFunction.h,v 1.6 2003/03/22 01:22:50 jchiang Exp $
 *
 */
    
class ProductFunction : public CompositeFunction {
public:

   ProductFunction(Function &a, Function &b);

   double value(Arg &x) const
      {return m_a->value(x)*m_b->value(x);}

   virtual Function* clone() const {
      return new ProductFunction(*this);
   }

protected:

   void fetchDerivs(Arg &x, std::vector<double> &derivs, bool getFree) const;

};

} // namespace Likelihood

#endif // ProductFunction_h
