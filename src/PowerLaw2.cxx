/** 
 * @file PowerLaw2.cxx
 * @brief Implementation for the PowerLaw2 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLaw2.cxx,v 1.8 2015/03/03 18:05:37 jchiang Exp $
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

double PowerLaw2::s_cX = 0;
double PowerLaw2::s_cLogX = 0;

PowerLaw2::PowerLaw2(double Integral, double Index, 
                     double LowerLimit, double UpperLimit) 
   : optimizers::Function("PowerLaw2", 4, "Integral") {
   addParam("Integral", Integral, true);
   addParam("Index", Index, true);
   addParam("LowerLimit", LowerLimit, false);
   addParam("UpperLimit", UpperLimit, false);

   setParamAlwaysFixed("LowerLimit");
   setParamAlwaysFixed("UpperLimit");

   m_cGamma = Index+1;
   m_cX = 0;
   updateCache(1.0, Index, LowerLimit, UpperLimit);
}

double PowerLaw2::value(const optimizers::Arg & xarg) const {
   double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();

// Assume a standard ordering for the parameters.
   enum ParamTypes {Integral, Index, LowerLimit, UpperLimit};

   double NN = m_parameter[Integral].getTrueValue();
   double gamma = m_parameter[Index].getTrueValue();
   double x1 = m_parameter[LowerLimit].getTrueValue();
   double x2 = m_parameter[UpperLimit].getTrueValue();

   updateCache(x, gamma, x1, x2);

   if (gamma == -1.0)
     return NN/x/(m_cLogXHi - m_cLogXLo);
   else
     return NN*m_cGXFact*m_cPowX;
}

double PowerLaw2::
derivByParamImp(const optimizers::Arg & xarg, const std::string & paramName) const {
   double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();

   int iparam(-1);
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         iparam = i;
	 break;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "PowerLaw2::derivByParam");
   }
   
   enum ParamTypes {Integral, Index, LowerLimit, UpperLimit};

   double NN = m_parameter[Integral].getTrueValue();
   double gamma = m_parameter[Index].getTrueValue();
   double x1 = m_parameter[LowerLimit].getTrueValue();
   double x2 = m_parameter[UpperLimit].getTrueValue();

   updateCache(x, gamma, x1, x2);

   double val = 0;
   switch (iparam) {
   case Integral:
      if (gamma == -1.) {
	 val = 1./x/(m_cLogXHi - m_cLogXLo);
      } else {
	 val = m_cGXFact*m_cPowX;
      }
      return val * m_parameter[Integral].getScale();
      break;
   case Index:
      if (gamma == -1.) {
	 val = -NN/(2.*x)*(m_cLogXHi+m_cLogXLo-2.*s_cLogX)/
	   (m_cLogXHi-m_cLogXLo);
      } else {
	 double one_p_gamma = 1.+gamma;
	 val = NN*(m_cPowXHi*(1.-one_p_gamma*(m_cLogXHi-s_cLogX)) -
		   m_cPowXLo*(1.-one_p_gamma*(m_cLogXLo-s_cLogX)))/
	   std::pow(m_cPowXHi - m_cPowXLo,2)*m_cPowX;
      }
      return val * m_parameter[Index].getScale();
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
   double xmin = dynamic_cast<const optimizers::dArg &>(x_min).getValue();
   double xmax = dynamic_cast<const optimizers::dArg &>(x_max).getValue();

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
