#include <vector>
#include <string>
#include <cmath>
#include "MyFun.h"

namespace Likelihood {

/* implement MyFun as a polynomial of degree m_maxNumParams-1 */

MyFun::MyFun() {
   setMaxNumParams(3);
   addParam(std::string("Ruthie"), 0.);
   addParam(std::string("Mary"), 0.);
   addParam(std::string("Jane"), 0.);
   addParam(std::string("Plain"), 3.14159);
}

double MyFun::value(double x) const {

   double my_val = 0.;
   std::vector<double> params;
   getParamValues(params);

   for (unsigned int i = 0; i < params.size(); i++) {
      my_val += params[i]*pow(x, i);
//      std::cerr << params[i] << "  " << i;
   }
//   std::cerr << std::endl;
   
   return my_val;
}

double MyFun::derivByParam(double x, const std::string &paramName) const {

   std::vector<string> my_paramName;
   getParamNames(my_paramName);
   std::vector<double> my_param;
   getParamValues(my_param);

   for (unsigned int i = 0; i < my_paramName.size(); i++) {
      if (paramName == my_paramName[i]) 
         return pow(x, i);
   }
   std::cerr << "Parameter " << paramName << " is not found."
             << std::endl;
   return 0.;
}

} // namespace Likelihood
