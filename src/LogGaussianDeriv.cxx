/** 
 * @file LogGaussianDeriv.cxx
 * @brief Implementation for the LogGaussianDeriv Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/LogGaussianDeriv.cxx,v 1.7 2009/01/19 15:18:18 sfegan Exp $
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/LogGaussianDeriv.h"

namespace Likelihood {

LogGaussianDeriv::LogGaussianDeriv(double norm, double mean, double sigma) {
   setMaxNumParams(3);

   addParam("Norm", norm, true);
   addParam("Mean", mean, true);
   addParam("Sigma", sigma, true);

// Set FuncType and ArgType for use with CompositeFunction hierarchy.
   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "LogGaussianDeriv";
   m_normParName = "Norm";
}

double LogGaussianDeriv::value(optimizers::Arg & xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   enum ParamTypes {Norm, Mean, Sigma};
   double norm = m_parameter[Norm].getTrueValue();
   double mean = m_parameter[Mean].getTrueValue();
   double sigma = m_parameter[Sigma].getTrueValue();
   return norm*(-(x - mean)/sigma/sigma);
}

double LogGaussianDeriv::
derivByParam(optimizers::Arg & xarg, const std::string & paramName) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   
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

   enum ParamTypes {Norm, Mean, Sigma};

   double norm = m_parameter[Norm].getTrueValue();
   double mean = m_parameter[Mean].getTrueValue();
   double sigma = m_parameter[Sigma].getTrueValue();
   double s3 = sigma*sigma*sigma;

   switch (iparam) {
   case Norm:
      return value(xarg)/norm;
      break;
   case Mean:
      return norm/sigma/sigma;
      break;
   case Sigma:
      return norm*(x - mean)/2./s3;
      break;
   default:
      break;
   }
   return 0;
}

double LogGaussianDeriv::integral(optimizers::Arg & x_min, 
                           optimizers::Arg & x_max) const {
   return 0;
}

} // namespace Likelihood
