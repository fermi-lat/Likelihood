/** 
 * @file Rosen.cxx
 * @brief Implementation for the 2D Rosenbrock objective function
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Rosen.cxx,v 1.7 2003/05/29 20:10:46 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include "Rosen.h"
#include "Likelihood/dArg.h"

namespace Likelihood {

void Rosen::init() {
   int nParams = 2;
   setMaxNumParams(nParams);

   addParam(std::string("x"), 1, true);
   addParam(std::string("y"), 1, true);
//   std::cout << "Rosen::init: " << m_parameter.size() << std::endl;
}

double Rosen::value(const std::vector<double> &params) {
   double x = params[0];
   double y = params[1];
   m_parameter[0].setValue(x);
   m_parameter[1].setValue(y);
//   std::cout << "Rosen::value: " << m_parameter.size() << std::endl;

   return -(m_prefactor*pow((y - x*x), 2) + pow((1 - x), 2));
}

void Rosen::getFreeDerivs(std::vector<double> &freeDerivs) {
   if (!freeDerivs.empty()) freeDerivs.clear();

//   std::cout << "Rosen::getFreeDervis: " << m_parameter.size() << std::endl;

   dArg xarg(0);
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      freeDerivs.push_back(-derivByParam(xarg, m_parameter[i].getName()));
   }
}

double Rosen::derivByParam(Arg &, 
                           const std::string &paramName) 
   throw(ParameterNotFound) {
   std::vector<double> params;
   getParamValues(params);

   double x = params[0];
   double y = params[1];

   if (paramName == "x") {
      return -4.*m_prefactor*(y - x*x)*x - 2.*(1. - x);
   } else if (paramName == "y") {
      return 2.*m_prefactor*(y - x*x);
   }
   throw ParameterNotFound(paramName, getName(), "Rosen::derivByParam");
}

} // namespace Likelihood
