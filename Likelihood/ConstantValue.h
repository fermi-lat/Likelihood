/** @file ConstantValue.h
 * @brief Declaration for the ConstantValue Function class
 * @author J. Chiang
 *
 * $Header$
 *
 */

#ifndef ConstantValue_h
#define ConstantValue_h

#include "Likelihood/Function.h"

namespace Likelihood {
/** 
 * @class ConstantValue
 *
 * @brief This returns a constant value, the sole Parameter of this class,
 * regardless of the value or type of Arg.
 *
 * @author J. Chiang
 *    
 * $Header$
 *
 */
    
class ConstantValue : public Function {
public:

   ConstantValue(double value) {
      setMaxNumParams(1);
      addParam("Value", value, true);

// need to double-check these...
      m_funcType = Factor;
      m_argType = "";
   }

   virtual ~ConstantValue() {}

   double value(Arg&) const {return m_parameter[0].getTrueValue();}

   double derivByParam(Arg &, const std::string &) const
      {return 0;}

   virtual Function *clone() const {
      return new ConstantValue(*this);
   }

private:

   // disable this
   double integral(Arg &, Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // ConstantValue_h
