/** 
 * @file BrokenPowerLaw.cxx
 * @brief Implementation for the BrokenPowerLaw Function class
 * @author J. Chiang
 *
 * $Header$
 */

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "optimizers/dArg.h"
#include "BrokenPowerLaw.h"

namespace Likelihood {

// initialization function used by constructors
void BrokenPowerLaw::init(double Prefactor, double Index1, double Index2,
                          double BreakValue) {
// Implement BrokenPowerLaw class with three named parameters, 
// "Prefactor", "Index1", "Index2", "BreakValue"

   int nParams = 4;
   setMaxNumParams(nParams);

   addParam(std::string("Prefactor"), Prefactor, true);
   addParam(std::string("Index1"), Index1, true);
   addParam(std::string("Index2"), Index2, true);
   addParam(std::string("BreakValue"), BreakValue, true);

// set FuncType and ArgType for use with CompositeFunction hierarchy
   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "BrokenPowerLaw";
}

double BrokenPowerLaw::value(optimizers::Arg &xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

// assume a standard ordering for the parameters
   enum ParamTypes {Prefactor, Index1, Index2, BreakValue};

   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   if (x < my_params[BreakValue].getTrueValue()) {
      return my_params[Prefactor].getTrueValue()
         *pow((x/my_params[BreakValue].getTrueValue()), 
              my_params[Index1].getTrueValue());
   } else {
      return my_params[Prefactor].getTrueValue()
         *pow((x/my_params[BreakValue].getTrueValue()), 
              my_params[Index2].getTrueValue());
   }
}

double BrokenPowerLaw::derivByParam(optimizers::Arg &xarg, 
                                    const std::string &paramName) const 
   throw(optimizers::ParameterNotFound) {

   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

   enum ParamTypes {Prefactor, Index1, Index2, BreakValue};

   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, 
                                          getName(), 
                                          "BrokenPowerLaw::derivByParam");
   }
   
   switch (iparam) {
   case Prefactor:
      return value(xarg)/my_params[Prefactor].getTrueValue()
         *my_params[Prefactor].getScale();
      break;
   case Index1:
      if (x < my_params[BreakValue].getTrueValue()) {
         return value(xarg)*log(x/my_params[BreakValue].getTrueValue())
            *my_params[Index1].getScale();
      } else {
         return 0;
      }
      break;
   case Index2:
      if (x < my_params[BreakValue].getTrueValue()) {
         return 0;
      } else {
         return value(xarg)*log(x/my_params[BreakValue].getTrueValue())
            *my_params[Index2].getScale();
      }
      break;
   case BreakValue:
      if (x < my_params[BreakValue].getTrueValue()) {
         return -value(xarg)*(my_params[Index1].getTrueValue())
            /(my_params[BreakValue].getTrueValue())
            *my_params[BreakValue].getScale();
      } else {
         return -value(xarg)*(my_params[Index2].getTrueValue())
            /(my_params[BreakValue].getTrueValue())
            *my_params[BreakValue].getScale();
      }
      break;
   default:
      break;
   }
   return 0;
}

} // namespace Likelihood
