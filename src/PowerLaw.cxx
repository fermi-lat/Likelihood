/** 
 * @file PowerLaw.cxx
 * @brief Implementation for the PowerLaw Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw.cxx,v 1.14 2003/08/06 20:52:07 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "optimizers/dArg.h"
#include "PowerLaw.h"

namespace Likelihood {

// initialization function used by constructors
void PowerLaw::init(double Prefactor, double Index, double Scale) {
// Implement PowerLaw class with three named parameters, 
// "Prefactor", "Scale", "Index"

   int nParams = 3;
   setMaxNumParams(nParams);

   addParam(std::string("Prefactor"), Prefactor, true);
   addParam(std::string("Index"), Index, true);
// scale should always be fixed
   addParam(std::string("Scale"), Scale, false);

// set FuncType and ArgType for use with CompositeFunction hierarchy
   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "PowerLaw";
}

double PowerLaw::value(optimizers::Arg &xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

//! assume a standard ordering for the parameters
   enum ParamTypes {Prefactor, Index, Scale};

   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   return my_params[Prefactor].getTrueValue()
      *pow((x/my_params[Scale].getTrueValue()), 
           my_params[Index].getTrueValue());
}

double PowerLaw::derivByParam(optimizers::Arg &xarg, const std::string &paramName) const 
   throw(optimizers::ParameterNotFound) {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

   enum ParamTypes {Prefactor, Index, Scale};

   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, 
                                          getName(), 
                                          "PowerLaw::derivByParam");
   }
   
   switch (iparam) {
   case Prefactor:
      return value(xarg)/my_params[Prefactor].getTrueValue()
         *my_params[Prefactor].getScale();
      break;
   case Index:
      return value(xarg)*log(x/my_params[Scale].getTrueValue())
         *my_params[Index].getScale();
      break;
   case Scale:  // shouldn't ever need this, nonetheless....
      return -value(xarg)*(my_params[Index].getTrueValue())
         /(my_params[Scale].getTrueValue())
         *my_params[Scale].getScale();
      break;
   default:
      break;
   }
   return 0;
}

double PowerLaw::integral(optimizers::Arg &xargmin, 
                          optimizers::Arg &xargmax) const {
   double xmin = dynamic_cast<optimizers::dArg &>(xargmin).getValue();
   double xmax = dynamic_cast<optimizers::dArg &>(xargmax).getValue();

   enum ParamTypes {Prefactor, Index, Scale};
   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   double f0 = my_params[Prefactor].getTrueValue();
   double Gamma = my_params[Index].getTrueValue();
   double x0 = my_params[Scale].getTrueValue();

   return f0/(Gamma+1.)*(pow((xmax/x0), Gamma+1.) - pow((xmin/x0), Gamma+1.));
}

} // namespace Likelihood
