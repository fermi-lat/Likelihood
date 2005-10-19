/** 
 * @file PowerLaw2.cxx
 * @brief Implementation for the PowerLaw2 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw2.cxx,v 1.3 2005/06/09 16:11:40 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/PowerLaw2.h"

namespace Likelihood {

PowerLaw2::PowerLaw2(double Integral, double Index, 
                     double LowerLimit, double UpperLimit) {
// Implement PowerLaw2 class with four named parameters, 
// "Integral", "Index", "LowerLimit", "UpperLimit"

   setMaxNumParams(4);

   addParam("Integral", Integral, true);
   addParam("Index", Index, true);
   addParam("LowerLimit", LowerLimit, false);
   addParam("UpperLimit", UpperLimit, false);

   setParamAlwaysFixed("LowerLimit");
   setParamAlwaysFixed("UpperLimit");

// Set FuncType and ArgType for use with CompositeFunction hierarchy.
   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "PowerLaw2";
}

double PowerLaw2::value(optimizers::Arg & xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

// Assume a standard ordering for the parameters.
   enum ParamTypes {Integral, Index, LowerLimit, UpperLimit};

   double NN = m_parameter[Integral].getTrueValue();
   double gamma = -m_parameter[Index].getTrueValue();
   double x1 = m_parameter[LowerLimit].getTrueValue();
   double x2 = m_parameter[UpperLimit].getTrueValue();
   
   if (gamma == 1.) {
      return NN/x/std::log(x2/x1);
   }
   double one_m_gam = 1. - gamma;

   return (NN*one_m_gam*std::pow(x, -gamma)
           /(std::pow(x2, one_m_gam) - std::pow(x1, one_m_gam)));
}

double PowerLaw2::
derivByParam(optimizers::Arg & xarg, const std::string & paramName) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

   int iparam(-1);
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         iparam = i;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "PowerLaw2::derivByParam");
   }
   
   enum ParamTypes {Integral, Index, LowerLimit, UpperLimit};

   double NN = m_parameter[Integral].getTrueValue();
   double gamma = -m_parameter[Index].getTrueValue();
   double x1 = m_parameter[LowerLimit].getTrueValue()/x;
   double x2 = m_parameter[UpperLimit].getTrueValue()/x;

   double one_m_gam(1. - gamma);
   double pow_x1(std::pow(x1, one_m_gam));
   double pow_x2(std::pow(x2, one_m_gam));

   switch (iparam) {
   case Integral:
      return one_m_gam/x/(pow_x2 - pow_x1)*m_parameter[Integral].getScale();
      break;
   case Index:
      return -( NN/x*(pow_x1*(1. - one_m_gam*std::log(x1)) -
                      pow_x2*(1. - one_m_gam*std::log(x2)))
                /(pow_x2 - pow_x1)/(pow_x2 - pow_x1) )
         *m_parameter[Index].getScale();
      break;
   case LowerLimit:
   case UpperLimit:
      throw std::runtime_error("PowerLaw2::derivByParam: attempt to "
                               "take derivative wrt a fixed parameter.");
      break;
   default:
      break;
   }
   return 0;
}

double PowerLaw2::integral(optimizers::Arg & x_min, 
                           optimizers::Arg & x_max) const {
   double xmin = dynamic_cast<optimizers::dArg &>(x_min).getValue();
   double xmax = dynamic_cast<optimizers::dArg &>(x_max).getValue();

   enum ParamTypes {Integral, Index, LowerLimit, UpperLimit};

   double NN = m_parameter[Integral].getTrueValue();
   double gamma = -m_parameter[Index].getTrueValue();
   double x1 = m_parameter[LowerLimit].getTrueValue();
   double x2 = m_parameter[UpperLimit].getTrueValue();

   double one_m_gam(1. - gamma);

   return NN*( (std::pow(xmax, one_m_gam) - std::pow(xmin, one_m_gam))/
               (std::pow(x2, one_m_gam) - std::pow(x1, one_m_gam)) );
}

} // namespace Likelihood
