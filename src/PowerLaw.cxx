#include <vector>
#include <string>
#include <cmath>
#include "PowerLaw.h"

namespace Likelihood {

//! Implement PowerLaw class with three named parameters, 
//! "Prefactor", "Scale", "Index"

void PowerLaw::m_init(const double Prefactor, const double Index,
		      const double Scale) {
//! initialization function used by constructors

   int nParams = 3;
   setMaxNumParams(nParams);

   setParam(string("Prefactor"), Prefactor, true);
   setParam(string("Index"), Index, true);
   setParam(string("Scale"), Scale, false);   // scale should always be fixed
}

//! copy constructor
PowerLaw::PowerLaw(const PowerLaw &func) {
   setMaxNumParams(func.getMaxNumParams());

   std::vector<Parameter> params = func.getParams();

   for (int i = 0; i < params.size(); i++)
      setParam(params[i]);
}

double PowerLaw::value(const double x) const {
//! assume a standard ordering for the parameters

   enum paramTypes {Prefactor, Index, Scale};

   std::vector<Parameter> my_params = getParams();

   return my_params[Prefactor].getValue()*pow((x/my_params[Scale].getValue()), 
					      my_params[Index].getValue());
}

double PowerLaw::derivByParam(const double x, 
			      const std::string paramName) const {

   enum paramTypes {Prefactor, Index, Scale};

   std::vector<Parameter> my_params = getParams();

   int iparam = -1;
   for (int i = 0; i < my_params.size(); i++) {
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
}

std::vector<double> PowerLaw::getDerivs(const double x) const {

   std::vector<string> my_paramName = getParamNames();
   std::vector<double> my_derivs;

   for (int i = 0; i < my_paramName.size(); i++) {
      my_derivs.push_back(derivByParam(x, my_paramName[i]));
   }
   return my_derivs;
}

std::vector<double> PowerLaw::getFreeDerivs(const double x) const {

   std::vector<Parameter> my_params = getParams();
   std::vector<double> my_derivs;

   for (int i = 0; i < my_params.size(); i++) {
      if (my_params[i].isFree())
	 my_derivs.push_back(derivByParam(x, my_params[i].getName()));
   }
   return my_derivs;
}

} // namespace Likelihood
