/** 
 * @file BandFunction.cxx
 * @brief Implementation for the BandFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BandFunction.cxx,v 1.1 2005/01/26 06:53:39 jchiang Exp $
 */

#include <cmath>

#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/BandFunction.h"

namespace {
   class Pars {
   public:
      Pars(const std::vector<optimizers::Parameter> & params) 
         : m_params(params) {}

      /// @return The true value of the parameter
      double operator[](unsigned int i) {
         return m_params[i].getTrueValue();
      }

      unsigned int size() const {
         return m_params.size();
      }

      const optimizers::Parameter & operator()(unsigned int i) const {
         return m_params[i];
      }
   private:
      const std::vector<optimizers::Parameter> & m_params;
   };
}

namespace Likelihood {

void BandFunction::init(double norm, double alpha, double beta, double Ep) {
   setMaxNumParams(4);

   addParam("norm", norm, true);
   addParam("alpha", alpha, true);
   addParam("beta", beta, true);
   addParam("Ep", Ep, true);

   m_funcType = Addend;
   m_argType = "dArg";

   m_genericName = "BandFunction";
}

double BandFunction::value(optimizers::Arg &xarg) const {
   ::Pars pars(m_parameter);

// energy is scaled by 100 keV:
   double energy = dynamic_cast<optimizers::dArg &>(xarg).getValue()/0.1;
   double epeak(pars[3]/0.1);
   double ebreak = epeak*(pars[1] - pars[2])/(pars[1] + 2.);

   double my_value;
   if (energy < ebreak) {
      my_value = pars[0]*std::pow(energy, pars[1])
         *std::exp(-energy*(2. + pars[1])/epeak);
   } else {
      my_value = pars[0]*std::pow(epeak*(pars[1] - pars[2])/(pars[1] + 2.),
                                  pars[1] - pars[2])
         *std::exp(pars[2] - pars[1])*std::pow(energy, pars[2]);
   }
   return my_value;
}

double BandFunction::derivByParam(optimizers::Arg & xarg,
                          const std::string & paramName) const {
   ::Pars pars(m_parameter);

   double energy = dynamic_cast<optimizers::dArg &>(xarg).getValue()/0.1;
   double epeak(pars[3]/0.1);
   double ebreak = epeak*(pars[1] - pars[2])/(pars[1] + 2.);

   int iparam = -1;
   for (unsigned int i = 0; i < pars.size(); i++) {
      if (paramName == pars(i).getName()) {
         iparam = i;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(), 
                                          "BandFunction::derivByParam");
   }
   
   enum ParamTypes {norm, alpha, beta, Ep};
   switch (iparam) {
   case norm:
      return value(xarg)/pars[0]*m_parameter[norm].getScale();
   case alpha:
      if (energy < ebreak) {
         return value(xarg)*(std::log(energy) - energy/epeak)
            *m_parameter[alpha].getScale();
      } else {
         return value(xarg)*(std::log(ebreak)*(epeak*(pars[2] + 2.)
                                               /(pars[1]+2.)/(pars[1]+2.)-1.))
            *m_parameter[alpha].getScale();
      }
   case beta:
      if (energy < ebreak) {
         return 0;
      } else {
         return value(xarg)*(1. + std::log(energy) + epeak/(pars[1] + 2.)
                             *std::log(ebreak))
            *m_parameter[beta].getScale();
      }
   case Ep:
      if (energy < ebreak) {
         return value(xarg)*(2. + pars[1])*energy/epeak
            *m_parameter[Ep].getScale();
      } else {
         return value(xarg)*(pars[1] - pars[2])/epeak
            *m_parameter[Ep].getScale();
      }
   default:
      break;
   }
   return 0;
}

} // namespace Likelihood
