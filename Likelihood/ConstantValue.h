/** 
 * @file ConstantValue.h
 * @brief Declaration for the ConstantValue Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ConstantValue.h,v 1.3 2003/07/19 04:38:01 jchiang Exp $
 *
 */

#ifndef Likelihood_ConstantValue_h
#define Likelihood_ConstantValue_h

#include "optimizers/Function.h"

namespace Likelihood {

/** 
 * @class ConstantValue
 *
 * @brief This returns a constant value, the sole Parameter of this class,
 * regardless of the value or type of Arg.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ConstantValue.h,v 1.3 2003/07/19 04:38:01 jchiang Exp $
 *
 */
    
class ConstantValue : public optimizers::Function {
public:

   ConstantValue(double value) {
      setMaxNumParams(1);
      addParam("Value", value, true);

// need to double-check these...
      m_funcType = Factor;
      m_argType = "";
   }

   virtual ~ConstantValue() {}

   double value(optimizers::Arg&) const {return m_parameter[0].getTrueValue();}

   double derivByParam(optimizers::Arg &, const std::string &) const
      {return 0;}

   virtual optimizers::Function *clone() const {
      return new ConstantValue(*this);
   }

private:

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // Likelihood_ConstantValue_h
