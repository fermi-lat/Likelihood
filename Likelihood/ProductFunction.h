/** 
 * @file ProductFunction.h
 * @brief Declaration of ProductFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ProductFunction.h,v 1.7 2003/06/11 17:08:02 jchiang Exp $
 */

#ifndef Likelihood_ProductFunction_h
#define Likelihood_ProductFunction_h

#include "Likelihood/CompositeFunction.h"

namespace Likelihood {
/** 
 * @class ProductFunction
 *
 * @brief A Function that returns the product of two Functions
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ProductFunction.h,v 1.7 2003/06/11 17:08:02 jchiang Exp $
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

#endif // Likelihood_ProductFunction_h
