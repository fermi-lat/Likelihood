/** 
 * @file MyFun.cxx
 * @brief Implementation of a simple test function 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/MyFun.cxx,v 1.7 2003/05/21 23:12:13 jchiang Exp $
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

// built-in unit test (justifies the existence of MyFun)
   try {
      addParam(std::string("Plain"), 3.14159);
   } catch (LikelihoodException &eObj) {
      std::cout << eObj.what() << std::endl;
   }
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

double MyFun::derivByParam(Arg &xarg, const std::string &paramName) const 
   throw(ParameterNotFound) {
   double x = dynamic_cast<dArg &>(xarg).getValue();

   std::vector<Parameter> params;
   getParams(params);

   for (unsigned int i = 0; i < params.size(); i++) {
      if (paramName == params[i].getName()) 
         return params[i].getScale()*pow(x, i);
   }
   throw ParameterNotFound(paramName, getName(), "MyFun::deriveByParam");
}

} // namespace Likelihood
