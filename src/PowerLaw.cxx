#include <vector>
#include <string>
#include <cmath>
#include "PowerLaw.h"

namespace Likelihood {

//! Implement PowerLaw class with three named parameters, 
//! "Prefactor", "Scale", "Index"

void PowerLaw::m_init(double Prefactor, double Index, double Scale) {
//! initialization function used by constructors

   int nParams = 3;
   setMaxNumParams(nParams);

   addParam(string("Prefactor"), Prefactor, true);
   addParam(string("Index"), Index, true);
   addParam(string("Scale"), Scale, false);   // scale should always be fixed
}

double PowerLaw::value(double x) const {
//! assume a standard ordering for the parameters

   enum paramTypes {Prefactor, Index, Scale};

   std::vector<Parameter> my_params;
   getParams(my_params);

   return my_params[Prefactor].getValue()*pow((x/my_params[Scale].getValue()), 
					      my_params[Index].getValue());
}

double PowerLaw::derivByParam(double x, const std::string &paramName) const {

   enum paramTypes {Prefactor, Index, Scale};

   std::vector<Parameter> my_params;
   getParams(my_params);

   int iparam = -1;
   for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
   }

   if (iparam == -1) {
// should throw an exception here
      std::cerr << "PowerLaw::derivByParam: "
		<< "Parameter " << paramName << " is not found."
		<< std::endl;
      return 0.;
   }
   
   switch (iparam) {
   case Prefactor:
      return value(x)/my_params[Prefactor].getValue();
      break;
   case Index:
      return value(x)*log(x/my_params[Scale].getValue());
      break;
   case Scale:  // shouldn't ever need this, nonetheless....
      return -value(x)*(my_params[Index].getValue())
	 /(my_params[Scale].getValue());
      break;
   default:
      break;
   }
   return 0;
}

double PowerLaw::integral(double xmin, double xmax) {

   enum paramTypes {Prefactor, Index, Scale};
   std::vector<Parameter> my_params;
   getParams(my_params);

   double f0 = my_params[Prefactor].getValue();
   double Gamma = my_params[Index].getValue();
   double x0 = my_params[Scale].getValue();

   return f0/(Gamma+1.)*(pow((xmax/x0), Gamma+1.) - pow((xmin/x0), Gamma+1.));
}

} // namespace Likelihood
