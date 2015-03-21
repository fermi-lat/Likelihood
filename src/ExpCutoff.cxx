/**
 * @file ExpCutoff.cxx
 * @brief Implementation for the ExpCutoff Function class
 * @author Luis C. Reyes
 *
 */

#include <cmath>

#include <iostream>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/ExpCutoff.h"

namespace Likelihood {

ExpCutoff::
ExpCutoff(double Prefactor, double Index, double Scale,
          double Ebreak, double P1, double P2, double P3)
   : optimizers::Function("ExpCutoff", 7, "Prefactor") {
   addParam("Prefactor", Prefactor, true);
   addParam("Index", Index, true);
   addParam("Scale", Scale, false);
   addParam("Ebreak", Ebreak, true);
   addParam("P1", P1, true);
   addParam("P2", P2, true);
   addParam("P3", P3, true);
}

double ExpCutoff::value(const optimizers::Arg &xarg) const {
   double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();

// assume a standard ordering for the parameters
   enum ParamTypes {Prefactor, Index, Scale, Ebreak, P1, P2, P3};

   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   double prefactor = my_params[Prefactor].getTrueValue();
   double index = my_params[Index].getTrueValue();
   double scale = my_params[Scale].getTrueValue();
   double ebreak = my_params[Ebreak].getTrueValue();
   double p1 = my_params[P1].getTrueValue();
   double p2 = my_params[P2].getTrueValue();
   double p3 = my_params[P3].getTrueValue();

   if (x < ebreak) {
      return prefactor*pow(x/scale, index);
   } else {
      double ln = log(x/ebreak);
      return prefactor*pow(x/scale, index)*exp(-( (x-ebreak)/p1 
                                                  + p2*ln + p3*ln*ln));
   }
}

double ExpCutoff::derivByParamImp(const optimizers::Arg & xarg,
                                  const std::string & paramName) const {
   double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();

   enum ParamTypes {Prefactor, Index, Scale, Ebreak, P1, P2, P3};

   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   double prefactor = my_params[Prefactor].getTrueValue();
   double index = my_params[Index].getTrueValue();
   double scale = my_params[Scale].getTrueValue();
   double ebreak = my_params[Ebreak].getTrueValue();
   double p1 = my_params[P1].getTrueValue();
   double p2 = my_params[P2].getTrueValue();
   double p3 = my_params[P3].getTrueValue();

   double ln = log(x/ebreak);


   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "ExpCutoff::derivByParam");
   }

   switch (iparam) {
   case Prefactor:
        return value(xarg)/prefactor*my_params[Prefactor].getScale();
      break;
   case Index:
      return value(xarg)*log(x/scale)*my_params[Index].getScale();
      break;
   case Scale:
      return -value(xarg)*(index/scale)*my_params[Scale].getScale();
      break;
   case Ebreak:
      if(x < ebreak) {
         return 0.;
      }	else {
         return value(xarg)*(1./p1 + p2/ebreak + 
                             (2.*p3*ln)/ebreak)*my_params[Ebreak].getScale();
      }
      break;
   case P1:
      if (x < ebreak) {
         return 0.;
      } else {
         return value(xarg)*(x-ebreak)/p1/p1*my_params[P1].getScale();
      }
   case P2:
      if (x < ebreak) {
         return 0.;
      } else {
         return -value(xarg)*ln*my_params[P2].getScale();
      }
   case P3:
      if (x < ebreak) {
         return 0.;
      } else {
         return -value(xarg)*ln*ln*my_params[P3].getScale();
      }
   default:
      break;
   }
   return 0.;
}

} // namespace optimizers
