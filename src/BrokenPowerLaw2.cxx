/** 
 * @file BrokenPowerLaw2.cxx
 * @brief Implementation for the BrokenPowerLaw2 Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BrokenPowerLaw2.cxx,v 1.3 2015/03/03 18:05:37 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/BrokenPowerLaw2.h"

namespace Likelihood {

BrokenPowerLaw2::BrokenPowerLaw2(double Integral, double Index1,
                                 double Index2, double BreakValue, 
                                 double LowerLimit, double UpperLimit) 
   : optimizers::Function("BrokenPowerLaw2", 6, "Integral") {
   addParam("Integral", Integral, true);
   addParam("Index1", Index1, true);
   addParam("Index2", Index2, true);
   addParam("BreakValue", BreakValue, true);
   addParam("LowerLimit", LowerLimit, false);
   addParam("UpperLimit", UpperLimit, false);

   setParamAlwaysFixed("LowerLimit");
   setParamAlwaysFixed("UpperLimit");
}

double BrokenPowerLaw2::value(const optimizers::Arg & xarg) const {
   double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();

// Assume a standard ordering for the parameters.
   enum ParamTypes {Integral, Index1, Index2, BreakValue, 
                    LowerLimit, UpperLimit};

   double gamma1 = -m_parameter[Index1].getTrueValue();
   double gamma2 = -m_parameter[Index2].getTrueValue();
   double x0 = m_parameter[BreakValue].getTrueValue();

   double N0 = N0_value();

   if (x < x0) {
      return N0*std::pow(x/x0, -gamma1);
   } else {
      return N0*std::pow(x/x0, -gamma2);
   }
}

double BrokenPowerLaw2::
derivByParamImp(const optimizers::Arg & xarg,
                const std::string & paramName) const {
   double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();

   int iparam(-1);
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         iparam = i;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "BrokenPowerLaw2::derivByParam");
   }
   
   enum ParamTypes {Integral, Index1, Index2, BreakValue,
                    LowerLimit, UpperLimit};

   double NN = m_parameter[Integral].getTrueValue();
   double gamma1 = -m_parameter[Index1].getTrueValue();
   double gamma2 = -m_parameter[Index2].getTrueValue();
   double x0 = m_parameter[BreakValue].getTrueValue();
   double x1 = m_parameter[LowerLimit].getTrueValue();
   double x2 = m_parameter[UpperLimit].getTrueValue();

   double N0 = N0_value();

   double one_m_gam1 = 1. - gamma1;
   double one_m_gam2 = 1. - gamma2;

   switch (iparam) {
   case Integral:
      double dN0dN;
      if (x1 > x0) {
         dN0dN = one_m_gam2/x0/(std::pow(x2/x0, one_m_gam2) - 
                                std::pow(x1/x0, one_m_gam2));
      } else if (x2 < x0) {
         dN0dN = one_m_gam1/x0/(std::pow(x2/x0, one_m_gam1) - 
                                std::pow(x1/x0, one_m_gam1));
      } else {
         dN0dN = 1./x0/((1. - std::pow(x1/x0, one_m_gam1))/one_m_gam1 + 
                        (std::pow(x2/x0, one_m_gam2) - 1.)/one_m_gam2);
      }
      if (x < x0) {
         return dN0dN*std::pow(x/x0, -gamma1)*m_parameter[Integral].getScale();
      }
      return dN0dN*std::pow(x/x0, -gamma2)*m_parameter[Integral].getScale();
      break;
   case Index1:
      double dN0dI1;
      if (x1 > x0) {
         return 0;
      } else if (x2 < x0) {
         double A1(std::pow(x1/x0, one_m_gam1));
         double A2(std::pow(x2/x0, one_m_gam1));
         double denom = A2 - A1;
         dN0dI1 = NN/x0/denom/denom*
            (denom - one_m_gam1*(A2*std::log(x2/x0) - A1*std::log(x1/x0)));
      } else {
         dN0dI1 = -N0*N0*x0/NN/one_m_gam1/one_m_gam1
            *(std::pow(x1/x0, one_m_gam1)*(1.-one_m_gam1*std::log(x1/x0))-1.);
      }
      if (x < x0) {
         return (dN0dI1 + N0*std::log(x/x0))*std::pow(x/x0, -gamma1)
            *m_parameter[Index1].getScale();
      }
      return dN0dI1*std::pow(x/x0, -gamma2)*m_parameter[Index1].getScale();
      break;
   case Index2:
      double dN0dI2;
      if (x2 < x0) {
         return 0;
      } else if (x1 > x0) {
         double A1(std::pow(x1/x0, one_m_gam2));
         double A2(std::pow(x2/x0, one_m_gam2));
         double denom = A2 - A1;
         dN0dI2 = NN/x0/denom/denom*
            (denom - one_m_gam2*(A2*std::log(x2/x0) - A1*std::log(x1/x0)));
      } else {
         dN0dI2 = N0*N0*x0/NN/one_m_gam2/one_m_gam2
            *(std::pow(x2/x0, one_m_gam2)*(1.-one_m_gam2*std::log(x2/x0))-1.);
      }
      if (x > x0) {
         return (dN0dI2 + N0*std::log(x/x0))*std::pow(x/x0, -gamma2)
            *m_parameter[Index2].getScale();
      }
      return dN0dI2*std::pow(x/x0, -gamma1)*m_parameter[Index2].getScale();
      break;
   case BreakValue:
      double dN0dx0;
      if (x2 < x0) {
         dN0dx0 = -N0*N0/NN*(1./one_m_gam1 - 1.)*(std::pow(x2/x0, one_m_gam1) -
                                                  std::pow(x1/x0, one_m_gam1));
      } else if (x1 > x0) {
         dN0dx0 = -N0*N0/NN*(1./one_m_gam2 - 1.)*(std::pow(x2/x0, one_m_gam2) -
                                                  std::pow(x1/x0, one_m_gam2));
      } else {
         dN0dx0 = -N0*N0/NN*( (1. - std::pow(x1/x0, one_m_gam1))/one_m_gam1 
                              + (std::pow(x2/x0, one_m_gam2) - 1.)/one_m_gam2
                              + std::pow(x1/x0, one_m_gam1) 
                              - std::pow(x2/x0, one_m_gam2) );
      }
      if (x < x0) {
         return std::pow(x/x0, -gamma1)*(dN0dx0 + N0*gamma1/x0)
            *m_parameter[BreakValue].getScale();
      }
      return std::pow(x/x0, -gamma2)*(dN0dx0 + N0*gamma2/x0)
         *m_parameter[BreakValue].getScale();
      break;
   case LowerLimit:
   case UpperLimit:
      throw std::runtime_error("BrokenPowerLaw2::derivByParam: attempt to "
                               "take derivative wrt a fixed parameter.");
      break;
   default:
      break;
   }
   return 0;
}

double BrokenPowerLaw2::integral(const optimizers::Arg & x_min, 
                                 const optimizers::Arg & x_max) const {
   double xmin = dynamic_cast<const optimizers::dArg &>(x_min).getValue();
   double xmax = dynamic_cast<const optimizers::dArg &>(x_max).getValue();
   double prefactor(1.);
   if (xmin > xmax) {
      double tmp = xmin;
      xmin = xmax;
      xmax = tmp;
      prefactor = -1;
   }

   enum ParamTypes {Integral, Index1, Index2, BreakValue,
                    LowerLimit, UpperLimit};

   double gamma1 = -m_parameter[Index1].getTrueValue();
   double gamma2 = -m_parameter[Index2].getTrueValue();
   double x0 = m_parameter[BreakValue].getTrueValue();

   double one_m_gam1(1. - gamma1);
   double one_m_gam2(1. - gamma2);

   double N0 = N0_value();

   double my_integral(0);
   if (xmin < x0 && xmax >= x0) {
      my_integral = N0*((1. - std::pow(xmin/x0, one_m_gam1))/one_m_gam1 +
                        (std::pow(xmax/x0, one_m_gam2) - 1.)/one_m_gam2);
   } else if (x0 < xmin) {
      my_integral = N0*(std::pow(xmax/x0, one_m_gam2) -
                        std::pow(xmin/x0, one_m_gam2))/one_m_gam2;
   } else { // xmax <= x0
      my_integral = N0*(std::pow(xmax/x0, one_m_gam1) -
                        std::pow(xmin/x0, one_m_gam1))/one_m_gam1;
   }
   return prefactor*my_integral;
}

double BrokenPowerLaw2::N0_value() const {
   enum ParamTypes {Integral, Index1, Index2, BreakValue, 
                    LowerLimit, UpperLimit};

   double NN = m_parameter[Integral].getTrueValue();
   double gamma1 = -m_parameter[Index1].getTrueValue();
   double gamma2 = -m_parameter[Index2].getTrueValue();
   double x0 = m_parameter[BreakValue].getTrueValue();
   double x1 = m_parameter[LowerLimit].getTrueValue();
   double x2 = m_parameter[UpperLimit].getTrueValue();

   double one_m_gam1 = 1. - gamma1;
   double one_m_gam2 = 1. - gamma2;

   if (x1 > x0) {
      return NN*one_m_gam2/x0/(std::pow(x2/x0, one_m_gam2) - 
                               std::pow(x1/x0, one_m_gam2));
   }
   if (x2 < x0) {
      return NN*one_m_gam1/x0/(std::pow(x2/x0, one_m_gam1) - 
                               std::pow(x1/x0, one_m_gam1));
   }
   return NN/x0/((1. - std::pow(x1/x0, one_m_gam1))/one_m_gam1 + 
                 (std::pow(x2/x0, one_m_gam2) - 1.)/one_m_gam2);
}

} // namespace Likelihood
