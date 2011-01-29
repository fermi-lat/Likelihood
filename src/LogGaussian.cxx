/** 
 * @file LogGaussian.cxx
 * @brief Implementation for the LogGaussian Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/LogGaussian.cxx,v 1.7 2009/01/19 15:18:18 sfegan Exp $
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/LogGaussian.h"

namespace Likelihood {

LogGaussian::LogGaussian(double norm, double mean, double sigma) {
   setMaxNumParams(3);

   addParam("Norm", norm, true);
   addParam("Mean", mean, true);
   addParam("Sigma", sigma, true);

// Set FuncType and ArgType for use with CompositeFunction hierarchy.
   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "LogGaussian";
   m_normParName = "Norm";

   m_log_term = std::log(std::sqrt(2.*M_PI)*sigma);
}

double LogGaussian::value(optimizers::Arg & xarg) const {
   double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
   enum ParamTypes {Norm, Mean, Sigma};
   double norm = m_parameter[Norm].getTrueValue();
   double mean = m_parameter[Mean].getTrueValue();
   double sigma = m_parameter[Sigma].getTrueValue();
   return norm*(-(x - mean)*(x - mean)/2./sigma/sigma - m_log_term);
}

double LogGaussian::
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
      return -norm*mean/sigma/sigma;
      break;
   case Sigma:
      return norm*((x - mean)*(x - mean)/6./s3 - 1./sigma);
      break;
   default:
      break;
   }
   return 0;
}

double LogGaussian::integral(optimizers::Arg & x_min, 
                           optimizers::Arg & x_max) const {
   return 0;
}

} // namespace Likelihood
