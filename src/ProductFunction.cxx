/** @file ProductFunction.cxx
 * @brief ProductFunction class implementation
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ProductFunction.cxx,v 1.2 2003/03/17 00:53:44 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "Likelihood/ProductFunction.h"

namespace Likelihood {

ProductFunction::ProductFunction(Function &a, Function &b) :
   CompositeFunction(a, b) {
   assert( (a.funcType() == Addend && b.funcType() == Factor) || 
           (a.funcType() == Factor && b.funcType() == Addend) || 
           (a.funcType() == Factor && b.funcType() == Factor) );
   if (a.funcType() == Addend || b.funcType() == Addend) {
      m_funcType = Addend;
   } else {
      m_funcType = Factor;
   }
   m_a = &a;
   m_b = &b;
   syncParams();
}

void ProductFunction::fetchDerivs(Arg &x, std::vector<double> &derivs, 
                                  bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   std::vector<double> my_derivs;
   if (getFree) {
      m_a->getFreeDerivs(x, my_derivs);
   } else {
      m_a->getDerivs(x, my_derivs);
   }
   for (unsigned int i = 0; i < my_derivs.size(); i++)
      derivs.push_back(my_derivs[i]*m_b->value(x));

   if (getFree) {
      m_b->getFreeDerivs(x, my_derivs);
   } else {
      m_b->getDerivs(x, my_derivs);
   }
   for (unsigned int i = 0; i < my_derivs.size(); i++)
      derivs.push_back(my_derivs[i]*m_a->value(x));
}

} // namespace Likelihood
