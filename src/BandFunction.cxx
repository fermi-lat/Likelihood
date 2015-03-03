/** 
 * @file BandFunction.cxx
 * @brief Implementation for the BandFunction class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BandFunction.cxx,v 1.5 2009/05/13 23:47:19 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
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

BandFunction::
BandFunction(double norm, double alpha, double beta, double Ep, double scale)
   : optimizers::Function("BandFunction", 5, "norm", "dArg", Addend) {
   addParam("norm", norm, true);
   addParam("alpha", alpha, true);
   addParam("beta", beta, true);
   addParam("Ep", Ep, true);
   addParam("Scale", scale, false);
}

double BandFunction::value(optimizers::Arg &xarg) const {
   ::Pars pars(m_parameter);

   double N0(pars[0]);
   double alpha(pars[1]);
   double beta(pars[2]);
   double scale(pars[4]);
   double energy = dynamic_cast<optimizers::dArg &>(xarg).getValue()/scale;
   double epeak(pars[3]/scale);
   double ebreak = epeak*(alpha - beta);

   if (alpha <= beta) {
      return N0*std::pow(energy, beta);
   }

   double my_value;
   if (energy < ebreak) {
      my_value = N0*std::pow(energy, alpha)*std::exp(-energy/epeak);
   } else {
      my_value = N0*std::pow(energy, beta)*std::pow(ebreak, alpha - beta)
         *std::exp(beta - alpha);
   }
   return my_value;
}

double BandFunction::derivByParamImp(optimizers::Arg & xarg,
                                     const std::string & paramName) const {
   ::Pars pars(m_parameter);

   double N0(pars[0]);
   double alpha(pars[1]);
   double beta(pars[2]);
   double scale(pars[4]);
   double energy = dynamic_cast<optimizers::dArg &>(xarg).getValue()/scale;
   double epeak(pars[3]/scale);
   double ebreak = epeak*(alpha - beta);

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
   
   enum ParamTypes {Norm, Alpha, Beta, Ep, Scale};
   switch (iparam) {
   case Norm:
      return value(xarg)/N0*m_parameter[Norm].getScale();
   case Alpha:
      if (alpha <= beta) {
         return N0*std::pow(energy, alpha)*std::log(energy)
            *m_parameter[Alpha].getScale();
      }
      if (energy < ebreak) {
         return value(xarg)*std::log(energy)
            *m_parameter[Alpha].getScale();
      } else {
         return value(xarg)*std::log(epeak*(alpha - beta))
            *m_parameter[Alpha].getScale();
      }
   case Beta:
      if (alpha <= beta) {
         return N0*std::pow(energy, beta)*std::log(energy)
            *m_parameter[Beta].getScale();
      }
      if (energy < ebreak) {
         return 0;
      } else {
         return value(xarg)*std::log(energy/epeak/(alpha - beta))*m_parameter[Beta].getScale();
      }
   case Ep:
      if (alpha <= beta) {
         return 0;
      }
      if (energy < ebreak) {
         return value(xarg)*energy/epeak/(epeak*scale)*m_parameter[Ep].getScale();
      } else {
         return value(xarg)*(alpha - beta)/(epeak*scale)
            *m_parameter[Ep].getScale();
      }
   case Scale:
      throw std::runtime_error("BandFunction::derivByParam: attempt to "
                               "take derivative wrt a fixed parameter.");
      break;
   default:
      break;
   }
   return 0;
}

} // namespace Likelihood
