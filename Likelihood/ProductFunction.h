/** @file ProductFunction.h
 * @brief Declaration of ProductFunction class
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef ProductFunction_h
#define ProductFunction_h

#include "Likelihood/CompositeFunction.h"

namespace Likelihood {
/** 
 * @class ProductFunction
 *
 * @brief A Function that returns the linear sum of two Functions
 *
 * @author J. Chiang
 *    
 * $Header$
 *
 */
    
class ProductFunction : public CompositeFunction {
public:

   ProductFunction(Function *a, Function *b) {
      m_a = a; m_b = b;}

   double value(Arg &x) const
      {return m_a->value(x)*m_b->value(x);}

protected:

  void fetchDerivs(Arg &x, std::vector<double> &derivs, bool getFree) const;

};

} // namespace Likelihood

#endif // ProductFunction_h
