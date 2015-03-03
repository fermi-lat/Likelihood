/** 
 * @file BrokenPowerLaw3.cxx
 * @brief Implementation for the BrokenPowerLaw3 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BrokenPowerLaw3.cxx,v 1.1 2012/07/31 19:38:35 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/BrokenPowerLaw3.h"

namespace Likelihood {

BrokenPowerLaw3::BrokenPowerLaw3(double Integral1, double Index1,
                                 double Integral2, double Index2, 
                                 double LowerLimit1, double UpperLimit1,
                                 double LowerLimit2, double UpperLimit2)
   : optimizers::Function("BrokenPowerLaw3", 8, "Integral1") {
   addParam("Integral1", Integral1, true);
   addParam("Index1", Index1, true);
   addParam("Integral2", Integral2, true);
   addParam("Index2", Index2, true);
   addParam("LowerLimit1", LowerLimit1, false);
   addParam("UpperLimit1", UpperLimit1, false);
   addParam("LowerLimit2", LowerLimit2, false);
   addParam("UpperLimit2", UpperLimit2, false);

   setParamAlwaysFixed("LowerLimit1");
   setParamAlwaysFixed("UpperLimit1");
   setParamAlwaysFixed("LowerLimit2");
   setParamAlwaysFixed("UpperLimit2");
}

double BrokenPowerLaw3::value(optimizers::Arg & xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

   enum ParamTypes {Integral1, Index1, Integral2, Index2, 
                    LowerLimit1, UpperLimit1, LowerLimit2, UpperLimit2};

   double index1(m_parameter[Index1].getTrueValue());
   double index2(m_parameter[Index2].getTrueValue());
   double x0(x0_value());
   double N0(N0_value(x0));

   if (x < x0) {
      return N0*std::pow(x/x0, index1);
   } else {
      return N0*std::pow(x/x0, index2);
   }
}

double BrokenPowerLaw3::
derivByParamImp(optimizers::Arg & xarg, const std::string & paramName) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();

   int iparam(-1);
   for (size_t i(0); i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         iparam = i;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "BrokenPowerLaw3::derivByParam");
   }
   
   enum ParamTypes {Integral1, Index1, Integral2, Index2, 
                    LowerLimit1, UpperLimit1, LowerLimit2, UpperLimit2};

   double F1(m_parameter[Integral1].getTrueValue());
   double F2(m_parameter[Integral2].getTrueValue());
   double index1(m_parameter[Index1].getTrueValue());
   double index2(m_parameter[Index2].getTrueValue());
   double x1min(m_parameter[LowerLimit1].getTrueValue());
   double x1max(m_parameter[UpperLimit1].getTrueValue());
   double x2min(m_parameter[LowerLimit2].getTrueValue());
   double x2max(m_parameter[UpperLimit2].getTrueValue());

   double x0(x0_value());
   double N0(N0_value(x0));

   double one_p_ind1 = 1. + index1;
   double one_p_ind2 = 1. + index2;

   double my_deriv(0);

   double dN0dF1, dx0dF1, dx0dF2;
   double arg, dargdI1, dx0dI1, dN0dI1, dargdI2, dx0dI2, dN0dI2;

   switch (iparam) {
   case Integral1:
      if (x >= x0) {
         return 0;
      }
      my_deriv = std::pow(x/x0, index1)*N0/F1;
      return my_deriv*m_parameter[Integral1].getScale();
      break;
   case Index1:
      if (x >= x0) {
         return 0;
      }
      dN0dI1 = (N0*(1./one_p_ind1 + std::log(x0)
                    - ((std::pow(x1max, one_p_ind1)*std::log(x1max) - 
                        std::pow(x1min, one_p_ind1)*std::log(x1min))
                       /(std::pow(x1max, one_p_ind1) - 
                         std::pow(x1min, one_p_ind1)))));
      my_deriv = std::pow(x/x0, index1)*(dN0dI1 + N0*std::log(x/x0));
      return my_deriv*m_parameter[Index1].getScale();
      break;
   case Integral2:
      if (x < x0) {
         return 0;
      }
      my_deriv = N0/F2*std::pow(x/x0, index2);
      return my_deriv*m_parameter[Integral2].getScale();
   case Index2:
      if (x < x0) {
         return 0;
      }
      dN0dI2 = (N0*(1./one_p_ind2 + std::log(x0)
                    - ((std::pow(x2max, one_p_ind2)*std::log(x2max) - 
                        std::pow(x2min, one_p_ind2)*std::log(x2min))
                       /(std::pow(x2max, one_p_ind2) - 
                         std::pow(x2min, one_p_ind2)))));
      my_deriv = std::pow(x/x0, index2)*(dN0dI2 + N0*std::log(x/x0));
      return my_deriv*m_parameter[Index2].getScale();
      break;
   case LowerLimit1:
   case UpperLimit1:
   case LowerLimit2:
   case UpperLimit2:
      throw std::runtime_error("BrokenPowerLaw3::derivByParam: attempt to "
                               "take derivative wrt a fixed parameter.");
      break;
   default:
      break;
   }
   return 0;
}

double BrokenPowerLaw3::integral(optimizers::Arg & x_min, 
                                 optimizers::Arg & x_max) const {
   double xmin = dynamic_cast<optimizers::dArg &>(x_min).getValue();
   double xmax = dynamic_cast<optimizers::dArg &>(x_max).getValue();

   double prefactor(1.);
   if (xmin > xmax) {
      double tmp = xmin;
      xmin = xmax;
      xmax = tmp;
      prefactor = -1;
   }

   enum ParamTypes {Integral1, Index1, Integral2, Index2, 
                    LowerLimit1, UpperLimit1, LowerLimit2, UpperLimit2};

   double F1(m_parameter[Integral1].getTrueValue());
   double F2(m_parameter[Integral2].getTrueValue());
   double index1(m_parameter[Index1].getTrueValue());
   double index2(m_parameter[Index2].getTrueValue());
   double x1min(m_parameter[LowerLimit1].getTrueValue());
   double x1max(m_parameter[UpperLimit1].getTrueValue());
   double x2min(m_parameter[LowerLimit2].getTrueValue());
   double x2max(m_parameter[UpperLimit2].getTrueValue());

   double x0(x0_value());
   double N0(N0_value(x0));

   double one_p_ind1(1. + index1);
   double one_p_ind2(1. + index2);

   double my_integral(0);
   if (xmin < x0 && xmax >= x0) {
      my_integral = N0*((1. - std::pow(xmin/x0, one_p_ind1))/one_p_ind1 +
                        (std::pow(xmax/x0, one_p_ind2) - 1.)/one_p_ind2);
   } else if (x0 < xmin) {
      my_integral = N0*(std::pow(xmax/x0, one_p_ind2) -
                        std::pow(xmin/x0, one_p_ind2))/one_p_ind2;
   } else { // xmax <= x0
      my_integral = N0*(std::pow(xmax/x0, one_p_ind1) -
                        std::pow(xmin/x0, one_p_ind1))/one_p_ind1;
   }
   return prefactor*my_integral;
}

double BrokenPowerLaw3::N0_value(double x0) const {
   enum ParamTypes {Integral1, Index1, Integral2, Index2, 
                    LowerLimit1, UpperLimit1, LowerLimit2, UpperLimit2};
   double F1(m_parameter[Integral1].getTrueValue());
   double index1(m_parameter[Index1].getTrueValue());
   double x1min(m_parameter[LowerLimit1].getTrueValue());
   double x1max(m_parameter[UpperLimit1].getTrueValue());
   
   double one_p_ind1(1. + index1);

   double N0 = (F1*one_p_ind1*std::pow(x0, index1)
                /(std::pow(x1max, one_p_ind1) - std::pow(x1min, one_p_ind1)));
   return N0;
}

double BrokenPowerLaw3::x0_value() const {
   enum ParamTypes {Integral1, Index1, Integral2, Index2, 
                    LowerLimit1, UpperLimit1, LowerLimit2, UpperLimit2};
   double F1(m_parameter[Integral1].getTrueValue());
   double F2(m_parameter[Integral2].getTrueValue());
   double index1(m_parameter[Index1].getTrueValue());
   double index2(m_parameter[Index2].getTrueValue());
   double x1min(m_parameter[LowerLimit1].getTrueValue());
   double x1max(m_parameter[UpperLimit1].getTrueValue());
   double x2min(m_parameter[LowerLimit2].getTrueValue());
   double x2max(m_parameter[UpperLimit2].getTrueValue());
   
   double one_p_ind1(1. + index1);
   double one_p_ind2(1. + index2);
   double x0(std::pow((F1/F2*one_p_ind1/one_p_ind2
                       *(std::pow(x2max,one_p_ind2)-std::pow(x2min,one_p_ind2))
                       /(std::pow(x1max,one_p_ind1)-std::pow(x1min,one_p_ind1))),
                      1./(index2 - index1)));
   return x0;
}

} // namespace Likelihood
