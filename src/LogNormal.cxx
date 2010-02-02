/**
 * @file LogNormal.cxx
 * @brief Log-Normal distribution for spectral modeling.  log10 is applied
 * to the energy value.
 *
 * $Header$
 */

#include <cmath>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "LogNormal.h"

namespace Likelihood {

LogNormal::LogNormal(double prefactor, double log10_mean,
                     double log10_sigma) {

   setMaxNumParams(3);
   addParam("Prefactor", prefactor, true);
   addParam("Log10_Mean", log10_mean, true);
   addParam("Log10_Sigma", log10_sigma, true);

   m_funcType = Addend;
   m_argType = "dArg";
   m_genericName = "LogNormal";
   m_normParName = "Prefactor";
}

double LogNormal::value(optimizers::Arg & xarg) const {
   double x(dynamic_cast<optimizers::dArg &>(xarg).getValue());
   double log10x(std::log10(x));

   enum ParamTypes {Prefactor, Log10_Mean, Log10_Sigma};

   std::vector<optimizers::Parameter> pars;
   getParams(pars);
   
   double foo((log10x - pars[Log10_Mean].getTrueValue())
              /pars[Log10_Sigma].getTrueValue());

   return pars[Prefactor].getTrueValue()/std::sqrt(2.*M_PI)/x
      /pars[Log10_Mean].getTrueValue()*std::exp(-foo*foo/2.);
}

double LogNormal::derivByParam(optimizers::Arg & xarg,
                               const std::string & parName) const {
   double x(dynamic_cast<optimizers::dArg &>(xarg).getValue());
   double log10x(std::log10(x));

   enum ParamTypes {Prefactor, Log10_Mean, Log10_Sigma};

   std::vector<optimizers::Parameter> pars;
   getParams(pars);

   int iparam(-1);
   for (size_t i(0); i < pars.size(); i++) {
      if (parName == pars[i].getName()) {
         iparam = i;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(parName, getName(), 
                                          "LogNormal::derivByParam");
   }

   double prefactor(pars[Prefactor].getTrueValue());
   double log10mean(pars[Log10_Mean].getTrueValue());
   double log10sigma(pars[Log10_Sigma].getTrueValue());

   double foo((log10x - log10mean)/log10sigma);
   switch (iparam) {
   case Prefactor:
      return pars[Prefactor].getScale()*value(xarg)/prefactor;
      break;
   case Log10_Mean:
      return pars[Log10_Mean].getScale()*value(xarg)*foo/log10sigma;
      break;
   case Log10_Sigma:
      return pars[Log10_Sigma].getScale()*value(xarg)/log10sigma
         *(foo*foo/log10sigma - 1.);
      break;
   default:
      break;
   }
   return 0;
}

} // namespace Likelihood
