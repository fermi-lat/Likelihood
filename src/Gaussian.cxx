#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Likelihood/dArg.h"
#include "Gaussian.h"

namespace Likelihood {

//! Implement Gaussian class with three named parameters, 
//! "Prefactor", "Mean", "Sigma"

void Gaussian::m_init(double Prefactor, double Mean, double Sigma) {
//! initialization function used by constructors

   int nParams = 3;
   setMaxNumParams(nParams);

   addParam(std::string("Prefactor"), Prefactor, true);
   addParam(std::string("Mean"), Mean, true);
   addParam(std::string("Sigma"), Sigma, true);
}

double Gaussian::integral(Arg &xargmin, Arg &xargmax) const {
   double xmin = dynamic_cast<dArg &>(xargmin).getValue();
   double xmax = dynamic_cast<dArg &>(xargmax).getValue();

   std::vector<Parameter> my_params;
   getParams(my_params);
   enum paramTypes {Prefactor, Mean, Sigma};

   double f0 = my_params[Prefactor].getTrueValue();
   double x0 = my_params[Mean].getTrueValue();
   double sigma = my_params[Sigma].getTrueValue();

   double zmin = (xmin - x0)/sqrt(2.)/sigma;
   double zmax = (xmax - x0)/sqrt(2.)/sigma;

   return f0*(m_erfcc(zmin) - m_erfcc(zmax))/2.;
}

double Gaussian::m_erfcc(double x) const {
/* (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y.. */
   double t, z, ans;

   z=fabs(x);
   t=1.0/(1.0+0.5*z);
   ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
         t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
         t*(-0.82215223+t*0.17087277)))))))));
   return x >= 0.0 ? ans : 2.0-ans;
}

double Gaussian::value(Arg &xarg) const {
   double x = dynamic_cast<dArg &>(xarg).getValue();

//! assume a standard ordering for the parameters
   enum paramTypes {Prefactor, Mean, Sigma};

   std::vector<Parameter> my_params;
   getParams(my_params);

   return my_params[Prefactor].getTrueValue()/sqrt(2.*M_PI)
      /my_params[Sigma].getTrueValue()
      *exp(-pow( (x - my_params[Mean].getTrueValue())
                 /my_params[Sigma].getTrueValue(), 2 )/2.);
}

double Gaussian::derivByParam(Arg &xarg, 
                              const std::string &paramName) const {
   double x = dynamic_cast<dArg &>(xarg).getValue();

   enum paramTypes {Prefactor, Mean, Sigma};

   std::vector<Parameter> my_params;
   getParams(my_params);

   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
   }

   if (iparam == -1) {
// should throw an exception here
      std::cerr << "Gaussian::derivByParam: "
                << "Parameter " << paramName << " is not found."
                << std::endl;
      return 0.;
   }
   
   switch (iparam) {
   case Prefactor:
      return value(xarg)/my_params[Prefactor].getTrueValue()
         *my_params[Prefactor].getScale();
      break;
   case Mean:
      return value(xarg)*(x - my_params[Mean].getTrueValue())
         /my_params[Sigma].getTrueValue()
         *my_params[Mean].getScale();
      break;
   case Sigma:
      return value(xarg)/my_params[Sigma].getTrueValue()
         *( pow((x - my_params[Mean].getTrueValue())
                /my_params[Sigma].getTrueValue(), 2) - 1. )
         *my_params[Sigma].getScale();
      break;
   default:
      break;
   }
   return 0;
}

} // namespace Likelihood
