#include <vector>
#include <string>
#include <cmath>
#include "MyFun.h"

namespace Likelihood {

/* implement MyFun as a polynomial of degree m_maxNumParams-1 */

double MyFun::value(const double x) const {

   double my_val = 0.;
   std::vector<double> params = getParamValues();

   for (int i = 0; i < params.size(); i++) {
      my_val += params[i]*pow(x, i);
//      std::cerr << params[i] << "  " << i;
   }
//   std::cerr << std::endl;
   
   return my_val;
}

double MyFun::derivByParam(const double x, 
			   const std::string paramName) const {

   std::vector<string> my_paramName = getParamNames();
   std::vector<double> my_param = getParamValues();

   for (int i = 0; i < my_paramName.size(); i++) {
      if (paramName == my_paramName[i]) 
         return pow(x, i);
   }
   std::cerr << "Parameter " << paramName << " is not found."
             << std::endl;
   return 0.;
}

std::vector<double> MyFun::getDerivs(const double x) const {

   std::vector<string> my_paramName = getParamNames();
   std::vector<double> my_derivs;

   for (int i = 0; i < my_paramName.size(); i++) {
      my_derivs.push_back(derivByParam(x, my_paramName[i]));
   }
   return my_derivs;
}

} // namespace Likelihood
