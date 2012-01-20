/*
 * @file SmoothBrokenDoublePowerLaw.cxx
 * @brief Implementation for the SmoothDoubleBrokenPowerLaw Function class
 * @author Keith Bechtol
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/SmoothDoubleBrokenPowerLaw.cxx,v 1.1 2012/01/12 16:46:56 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/SmoothDoubleBrokenPowerLaw.h"

namespace Likelihood {
  
void SmoothDoubleBrokenPowerLaw::init(double Prefactor,
                                      double Index1,
                                      double Scale,
                                      double Index2,
                                      double BreakValue12,
                                      double Beta12,
                                      double Index3,
                                      double BreakValue23,
                                      double Beta23) {
    
   // Implement PowerLaw class with nine named parameters,
   // "Prefactor", "Index1", "Scale", "Index2", "BreakValue12(Energy)",
   // "Beta12", "Index3", "BreakValue23(Energy)", "Beta23"
    
   int nParams = 9;
   setMaxNumParams(nParams);
    
   addParam(std::string("Prefactor"), Prefactor, true);
   addParam(std::string("Index1"), Index1, true); 
   addParam(std::string("Scale"), Scale, false);
   
   addParam(std::string("Index2"), Index2, true);
   addParam(std::string("BreakValue12"), BreakValue12, true);
   addParam(std::string("Beta12"), Beta12, true);

   addParam(std::string("Index3"), Index3, true);
   addParam(std::string("BreakValue23"), BreakValue23, true);
   addParam(std::string("Beta23"), Beta23, true);
    
   // Set FuncType and ArgType for use with CompositeFunction hierarchy.
   m_funcType = Addend;
   m_argType = "dArg";
   m_genericName = "SmoothDoubleBrokenPowerLaw";
   m_normParName = "Prefactor";
} 
  
double SmoothDoubleBrokenPowerLaw::value(optimizers::Arg &xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   
   // assume a standard ordering for the parameters
   enum ParamTypes {Prefactor, Index1, Scale, Index2, BreakValue12,
                    Beta12, Index3, BreakValue23, Beta23};
   
   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);
   
   double prefactor = my_params[Prefactor].getTrueValue();
   double index1 = my_params[Index1].getTrueValue();
   double scale = my_params[Scale].getTrueValue();
   double index2 = my_params[Index2].getTrueValue();
   double breakvalue12 = my_params[BreakValue12].getTrueValue();
   double beta12 = my_params[Beta12].getTrueValue();
   double index3 = my_params[Index3].getTrueValue();
   double breakvalue23 = my_params[BreakValue23].getTrueValue();
   double beta23 = my_params[Beta23].getTrueValue();
    
   return prefactor*pow(x/scale, index1) 
      *pow(1 + pow(x/breakvalue12, (index1 - index2)/beta12), -beta12) 
      *pow(1 + pow(x/breakvalue23, (index2 - index3)/beta23), -beta23);
}
  
double SmoothDoubleBrokenPowerLaw::
derivByParam(optimizers::Arg &xarg,
             const std::string &paramName) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   
   enum ParamTypes {Prefactor, Index1, Scale, Index2, BreakValue12,
                    Beta12, Index3, BreakValue23, Beta23};
    
   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);
    
   double prefactor = my_params[Prefactor].getTrueValue();
   double index1 = my_params[Index1].getTrueValue();
   double scale = my_params[Scale].getTrueValue();
   double index2 = my_params[Index2].getTrueValue();
   double breakvalue12 = my_params[BreakValue12].getTrueValue();
   double beta12 = my_params[Beta12].getTrueValue();
   double index3 = my_params[Index3].getTrueValue();
   double breakvalue23 = my_params[BreakValue23].getTrueValue();
   double beta23 = my_params[Beta23].getTrueValue();
   
   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) {
         iparam = i;
      }
   }
    
   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "SmoothDoubleBrokenPowerLaw::derivByParam");
   }
   double f1 = pow(x/breakvalue12, (index1 - index2)/beta12);
   double f2 = pow(x/breakvalue23, (index2 - index3)/beta23);

   switch (iparam) {
   case Prefactor:
      return value(xarg)/prefactor*my_params[Prefactor].getScale();
      break;
   case Index1:
      return value(xarg)*(log(x/scale) - log(x/breakvalue12)*f1/(1 + f1))
         *my_params[Index1].getScale();
      break;
   case Scale:
      return -value(xarg)*index1/scale*my_params[Scale].getScale();
      break;
   case Index2:
      return value(xarg)*(log(x/breakvalue12)*f1/(1 + f1)
                          - log(x/breakvalue23)*f2/(1 + f2));
      break;
   case BreakValue12:
      return value(xarg) *(index1-index2) *pow(x/breakvalue12,(index1-index2)/beta12)/(1+ pow(x/breakvalue12,(index1-index2)/beta12))/breakvalue12 *my_params[BreakValue12].getScale();
      break;
   case Beta12:
      return  value(xarg) * (-log(1+ pow(x/breakvalue12,(index1-index2)/beta12))+(index1-index2)*log(x/breakvalue12)/beta12/(1+ pow(x/breakvalue12,-(index1-index2)/beta12))) * my_params[Beta12].getScale();
      break;
   case Index3:
      return value(xarg) *log(x/breakvalue23)*pow(x/breakvalue23,(index2-index3)/beta23)/(1+ pow(x/breakvalue23,(index2-index3)/beta23))* my_params[Index3].getScale();
      break;
   case BreakValue23:
      return value(xarg) *(index2-index3) *pow(x/breakvalue23,(index2-index3)/beta23)/(1+ pow(x/breakvalue23,(index2-index3)/beta23))/breakvalue23 *my_params[BreakValue23].getScale();
      break;
   case Beta23:
      return  value(xarg) * (-log(1+ pow(x/breakvalue23,(index2-index3)/beta23))+(index2-index3)*log(x/breakvalue23)/beta23/(1+ pow(x/breakvalue23,-(index2-index3)/beta23))) * my_params[Beta23].getScale();
      break;
   default:
      break;
   }
    
   return 0;
}

} // namespace Likelihood
