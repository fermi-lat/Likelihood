/**
 * @file ExpCutoff2.cxx
 * @brief Implementation for the BrokenPowerLawExpCutoff Function class
 * @author Jennifer Carson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BrokenPowerLawExpCutoff.cxx,v 1.4 2007/07/13 15:35:11 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/BrokenPowerLawExpCutoff.h"

namespace Likelihood {

BrokenPowerLawExpCutoff::
BrokenPowerLawExpCutoff(double Prefactor, double Index1, double Index2, 
                        double BreakValue, double Eabs, double P1) 
   : optimizers::Function("BPLExpCutoff", 6, "Prefactor", "dArg", Addend) {
   addParam(std::string("Prefactor"), Prefactor, true);
   addParam(std::string("Index1"), Index1, true);
   addParam(std::string("Index2"), Index2, true);
   addParam(std::string("BreakValue"), BreakValue, false);
   addParam(std::string("Eabs"), Eabs, true);
   addParam(std::string("P1"), P1, true);
}

double BrokenPowerLawExpCutoff::value(optimizers::Arg & xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

// assume a standard ordering for the parameters
   enum ParamTypes {Prefactor, Index1, Index2, BreakValue, Eabs, P1};

   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   double prefactor = my_params[Prefactor].getTrueValue();
   double index1 = my_params[Index1].getTrueValue();
   double index2 = my_params[Index2].getTrueValue();
   double breakvalue = my_params[BreakValue].getTrueValue();
   double eabs = my_params[Eabs].getTrueValue();
   double p1 = my_params[P1].getTrueValue();

   if (x < breakvalue && x < eabs) {
      return prefactor*pow(x/breakvalue, index1);
   } else if (x > breakvalue && x < eabs) {
      return prefactor*pow(x/breakvalue, index2);
   } else if (x < breakvalue && x > eabs) {
      return prefactor*pow(x/breakvalue, index1)*exp(-( (x-eabs)/p1));
   } else if (x > breakvalue && x > eabs) {
      return prefactor*pow(x/breakvalue, index2)*exp(-( (x-eabs)/p1));
   }
   return 0;
}

double BrokenPowerLawExpCutoff::
derivByParamImp(optimizers::Arg & xarg,
                const std::string & paramName) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

   enum ParamTypes {Prefactor, Index1, Index2, BreakValue, Eabs, P1};

   std::vector<optimizers::Parameter> my_params;
   getParams(my_params);

   double pref = my_params[Prefactor].getTrueValue();
   double prefscale = my_params[Prefactor].getScale();
   double ind1 = my_params[Index1].getTrueValue();
   double ind1scale = my_params[Index1].getScale();
   double ind2 = my_params[Index2].getTrueValue();
   double ind2scale = my_params[Index2].getScale();
   double breakvalue = my_params[BreakValue].getTrueValue();
   double breakscale = my_params[BreakValue].getScale();
   double eabs = my_params[Eabs].getTrueValue();
   double eabsscale = my_params[Eabs].getScale();
   double p1 = my_params[P1].getTrueValue();
   double p1scale = my_params[P1].getScale();

   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "BrokenPowerLawExpCutoff::derivByParam");
   }

   switch (iparam) {
   case Prefactor:
      if (pref != 0) {
        return value(xarg)/pref*prefscale;
      } else {
	if (x < breakvalue && x < eabs) {
	  return std::pow(x/breakvalue, ind1)*prefscale;
	} else if (x > breakvalue && x < eabs) {
	  return std::pow(x/breakvalue, ind2)*prefscale;
	} else if (x < breakvalue && x > eabs) {
	  return std::pow(x/breakvalue, ind1)*prefscale 
	    * exp(-(x-eabs*eabsscale)/p1*p1scale);
	} else if (x > breakvalue && x > eabs) {
	  return std::pow(x/breakvalue, ind2)*prefscale
	    *exp(-(x-eabs*eabsscale)/p1*p1scale);
	}
      }
      break;
   case Index1:
     if (x < breakvalue) {
     return value(xarg)*std::log(x/breakvalue)*ind1scale;
     } else {
       return 0;
     }
     break;
   case Index2:
     if (x > breakvalue) {
     return value(xarg)*std::log(x/breakvalue)*ind2scale;
     } else {
       return 0;
     }
     break;
   case BreakValue:
     if (x < breakvalue) {
       return -value(xarg)*ind1/breakvalue*breakscale;
     } else {
       return -value(xarg)*ind2/breakvalue*breakscale;
     }
      break;
   case Eabs:
      if(x < eabs) {
         return 0.;
      }	else {
         return value(xarg)*(1./p1)*eabsscale;
      }
      break;
   case P1:
      if (x < eabs) {
         return 0.;
      } else {
         return value(xarg)*(x-eabs)/p1/p1*p1scale;
      }
   default:
      break;
   }
   return 0.;
}

} // namespace Likelihood
