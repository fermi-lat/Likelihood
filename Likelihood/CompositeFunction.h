/** @file CompositeFunction.h
 * @brief Declaration of CompositeFunction class
 * @author J. Chiang
 * $Header$
 */

#ifndef CompositeFunction_h
#define CompositeFunction_h

#include "Likelihood/Function.h"

namespace Likelihood {
/** 
 * @class CompositeFunction
 *
 * @brief Base class for Functions that are composites (sum or product)
 * of two other Functions.
 *
 * At some point, some sort of type-checking mechanism should
 * implemented to ensure that only Functions that operate on the same
 * Arg subclasses are combined.
 *
 * @author J. Chiang
 *    
 * $Header$
 * 
 */
    
class CompositeFunction : public Function {
public:

   CompositeFunction(){};
   virtual ~CompositeFunction(){};

   //! setParam method to include function name checking
   virtual void setParam(const Parameter &param, const std::string &funcName);

   //! group parameter access (note name mangling for inheritance 
   //! from Function)
   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator it) {
      it = m_a->setParamValues_(it);
      it = m_b->setParamValues_(it);
      return it;
   }

   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator it) {
      it = m_a->setFreeParamValues_(it);
      it = m_b->setFreeParamValues_(it);
      return it;
   }

   //! Parameter access including Function name specification
   virtual Parameter* getParam(const std::string &paramName, 
                               const std::string &funcName) const;
   
protected:

   // pointers to the Functions forming the composite
   Function *m_a;
   Function *m_b;

   //! method to sync the m_parameter vector with those of the two Functions
   void syncParams();

private:

   //! disable this since Parameters may no longer have unique names
   double derivByParam(Arg &, const std::string &) const {return 0;}

};

} // namespace Likelihood

#endif // CompositeFunction_h
