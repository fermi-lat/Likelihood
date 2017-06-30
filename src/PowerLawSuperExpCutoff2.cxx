/*
 * @file PowerLawSuperExpCutoff2.cxx
 * @brief Implementation for the PowerLawSuperExpCutoff2 Function class
 * @author Matthew Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PowerLawSuperExpCutoff2.cxx,v 1.2 2017/06/20 01:50:43 mdwood Exp $
 */

#include <cmath>

#include <iostream>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/PowerLawSuperExpCutoff2.h"

namespace Likelihood {
  
  PowerLawSuperExpCutoff2::
  PowerLawSuperExpCutoff2(double Prefactor,  
                         double Index1,
                         double Scale,
                         double Expfactor,
                         double Index2)
     : optimizers::Function("PLSuperExpCutoff2", 5, "Prefactor") {
    addParam(std::string("Prefactor"), Prefactor, true);
    addParam(std::string("Index1"), Index1, true);
    addParam(std::string("Scale"), Scale, false);
    addParam(std::string("Expfactor"), Expfactor, true);
    addParam(std::string("Index2"), Index2, true);
  } 
  
  double PowerLawSuperExpCutoff2::value(const optimizers::Arg &xarg) const {
    double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
    
    // assume a standard ordering for the parameters
    enum ParamTypes {Prefactor, Index1, Scale, Expfactor, Index2};
    
    std::vector<optimizers::Parameter> my_params;
    getParams(my_params);
    
    double prefactor = my_params[Prefactor].getTrueValue();
    double index1 = my_params[Index1].getTrueValue();
    double scale = my_params[Scale].getTrueValue();
    double expfactor = my_params[Expfactor].getTrueValue();
    double index2 = my_params[Index2].getTrueValue();

    // max value for exp argument to avoid overflow
    const double max_expfactor = std::log(1E100);
    return prefactor * pow(x/scale,index1) * 
      exp( std::min(- expfactor * pow(x,index2), max_expfactor) );
  }
  
  
  double PowerLawSuperExpCutoff2::derivByParamImp(const optimizers::Arg &xarg,
                                                 const std::string &paramName) const {
    double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
    
    enum ParamTypes {Prefactor, Index1, Scale, Expfactor, Index2};
    
    std::vector<optimizers::Parameter> my_params;
    getParams(my_params);
    
    double prefactor = my_params[Prefactor].getTrueValue();
    double index1 = my_params[Index1].getTrueValue();
    double scale = my_params[Scale].getTrueValue();
    double expfactor = my_params[Expfactor].getTrueValue();
    double index2 = my_params[Index2].getTrueValue();
    
    int iparam = -1;
    for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
    }
    
    if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "PowerLawSuperExpCutoff2::derivByParam");
    }
    
    switch (iparam) 
      {
      case Prefactor:
	return value(xarg) / prefactor * my_params[Prefactor].getScale();
	break;
      case Index1:
	return value(xarg) * log(x/scale) * my_params[Index1].getScale();
	break;
      case Scale:
	return -value(xarg) * index1 / scale * my_params[Scale].getScale();
	break;
      case Expfactor:
	return value(xarg) * ( - pow(x,index2) ) * my_params[Expfactor].getScale();
	break;
      case Index2:
	return -value(xarg) * expfactor * pow(x,index2) * log(x) * my_params[Index2].getScale();
	break;
      default:
	break;
      }
    
    return 0;
    
  }
  
} // namespace Likelihood
