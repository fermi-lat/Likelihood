/** @file MyFun.cxx
 * @brief Implementation a simple test function 
 * @author J. Chiang
 *
 * $Header$
 */

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include "Likelihood/dArg.h"
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

double MyFun::value(Arg &xarg) const {
   double x = dynamic_cast<dArg &>(xarg).getValue();

   double my_val = 0.;
   std::vector<Parameter> params;
   getParams(params);

   for (unsigned int i = 0; i < params.size(); i++) {
      my_val += params[i].getTrueValue()*pow(x, i);
   }
   
   return my_val;
}

double MyFun::derivByParam(Arg &xarg, const std::string &paramName) const {
   double x = dynamic_cast<dArg &>(xarg).getValue();

   std::vector<Parameter> params;
   getParams(params);

   for (unsigned int i = 0; i < params.size(); i++) {
      if (paramName == params[i].getName()) 
         return params[i].getScale()*pow(x, i);
   }
   std::cerr << "Parameter " << paramName << " is not found."
             << std::endl;
   return 0.;
}

} // namespace Likelihood
