/** 
 * @file LogParabola.cxx
 * @brief Implementation for the LogParabola class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LogParabola.cxx,v 1.5 2015/03/03 18:05:37 jchiang Exp $
 */

#include <cmath>

#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/LogParabola.h"

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

LogParabola::
LogParabola(double norm, double alpha, double beta, double Eb)
   : optimizers::Function("LogParabola", 4, "norm") {
   addParam("norm", norm, true);
   addParam("alpha", alpha, true);
   addParam("beta", beta, true);
   addParam("Eb", Eb, true);
}

double LogParabola::value(const optimizers::Arg & xarg) const {
   ::Pars pars(m_parameter);

   double energy = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
   double x = energy/pars[3];
   double my_value = pars[0]*std::pow(x, -(pars[1] + pars[2]*std::log(x)));
   return my_value;
}

double LogParabola::derivByParamImp(const optimizers::Arg & xarg,
                                    const std::string & paramName) const {
   ::Pars pars(m_parameter);

   double energy = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
   double x = energy/pars[3];
   double logx = std::log(x);
   double dfdnorm = std::pow(x, -(pars[1] + pars[2]*logx));

   int iparam = -1;
   for (unsigned int i = 0; i < pars.size(); i++) {
      if (paramName == pars(i).getName()) {
         iparam = i;
      }
   }

   if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(), 
                                          "LogParabola::derivByParam");
   }
   
   enum ParamTypes {norm, alpha, beta, Eb};
   switch (iparam) {
   case norm:
      return dfdnorm*m_parameter[norm].getScale();
   case alpha:
      return -pars[0]*logx*dfdnorm*m_parameter[alpha].getScale();
   case beta:
      return -pars[0]*logx*logx*dfdnorm*m_parameter[beta].getScale();
   case Eb:
      return value(xarg)/pars[3]*(pars[1] + 2.*pars[2]*logx)
         *m_parameter[Eb].getScale();
   default:
      break;
   }
   return 0;
}

} // namespace Likelihood
