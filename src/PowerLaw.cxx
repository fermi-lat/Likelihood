/** @file PowerLaw.cxx
 * @brief Implementation for the PowerLaw Function class
 * @author J. Chiang
 *
 * $Header$
 */

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Likelihood/dArg.h"
#include "PowerLaw.h"

namespace Likelihood {

//! Implement PowerLaw class with three named parameters, 
//! "Prefactor", "Scale", "Index"

void PowerLaw::init(double Prefactor, double Index, double Scale) {
//! initialization function used by constructors

   int nParams = 3;
   setMaxNumParams(nParams);

   addParam(std::string("Prefactor"), Prefactor, true);
   addParam(std::string("Index"), Index, true);
// scale should always be fixed
   addParam(std::string("Scale"), Scale, false);   
}

double PowerLaw::value(Arg &xarg) const {
   double x = dynamic_cast<dArg &>(xarg).getValue();

//! assume a standard ordering for the parameters
   enum ParamTypes {Prefactor, Index, Scale};

   std::vector<Parameter> my_params;
   getParams(my_params);

   return my_params[Prefactor].getTrueValue()
      *pow((x/my_params[Scale].getTrueValue()), 
           my_params[Index].getTrueValue());
}

double PowerLaw::derivByParam(Arg &xarg, const std::string &paramName) const {
   double x = dynamic_cast<dArg &>(xarg).getValue();

   enum ParamTypes {Prefactor, Index, Scale};

   std::vector<Parameter> my_params;
   getParams(my_params);

   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
   }

   if (iparam == -1) {
// should throw an exception here
      std::cerr << "PowerLaw::derivByParam: "
                << "Parameter " << paramName << " is not found."
                << std::endl;
      return 0.;
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

double PowerLaw::integral(Arg &xargmin, Arg &xargmax) const {
   double xmin = dynamic_cast<dArg &>(xargmin).getValue();
   double xmax = dynamic_cast<dArg &>(xargmax).getValue();

   enum ParamTypes {Prefactor, Index, Scale};
   std::vector<Parameter> my_params;
   getParams(my_params);

   double f0 = my_params[Prefactor].getTrueValue();
   double Gamma = my_params[Index].getTrueValue();
   double x0 = my_params[Scale].getTrueValue();

   return f0/(Gamma+1.)*(pow((xmax/x0), Gamma+1.) - pow((xmin/x0), Gamma+1.));
}

} // namespace Likelihood
