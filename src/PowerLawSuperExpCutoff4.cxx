/*
 * @file PowerLawSuperExpCutoff4.cxx
 * @brief Implementation for the PowerLawSuperExpCutoff4 Function class
 * @author Jean Ballet, Philippe Bruel
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/PowerLawSuperExpCutoff4.cxx,v 1.3 2017/06/30 23:50:46 mdwood Exp $
 */

#include <cmath>

#include <iostream>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/PowerLawSuperExpCutoff4.h"

namespace Likelihood {
  
  PowerLawSuperExpCutoff4::
  PowerLawSuperExpCutoff4(double Prefactor,  
                         double IndexS,
                         double Scale,
                         double ExpfactorS,
                         double Index2)
     : optimizers::Function("PLSuperExpCutoff4", 5, "Prefactor") {
    addParam(std::string("Prefactor"), Prefactor, true);
    addParam(std::string("IndexS"), IndexS, true);
    addParam(std::string("Scale"), Scale, false);
    addParam(std::string("ExpfactorS"), ExpfactorS, true);
    addParam(std::string("Index2"), Index2, true);
  } 
  
  double PowerLawSuperExpCutoff4::value(const optimizers::Arg &xarg) const {
    double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
    
    // assume a standard ordering for the parameters
    enum ParamTypes {Prefactor, IndexS, Scale, ExpfactorS, Index2};
    
    std::vector<optimizers::Parameter> my_params;
    getParams(my_params);
    
    double prefactor = my_params[Prefactor].getTrueValue();
    double indexS = my_params[IndexS].getTrueValue();
    double scale = my_params[Scale].getTrueValue();
    double expfactorS = my_params[ExpfactorS].getTrueValue();
    double index2 = my_params[Index2].getTrueValue();

    // max value for exp argument to avoid overflow
    const double max_expfactor = std::log(1E100);
    return prefactor * exp( std::min((indexS+expfactorS/index2) * log(x/scale) + expfactorS/index2/index2 * (1-pow(x/scale,index2) ), max_expfactor) );
  }
  
  
  double PowerLawSuperExpCutoff4::derivByParamImp(const optimizers::Arg &xarg,
                                                 const std::string &paramName) const {
    double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();

    enum ParamTypes {Prefactor, IndexS, Scale, ExpfactorS, Index2};
    
    std::vector<optimizers::Parameter> my_params;
    getParams(my_params);
    
    double prefactor = my_params[Prefactor].getTrueValue();
    double indexS = my_params[IndexS].getTrueValue();
    double scale = my_params[Scale].getTrueValue();
    double expfactorS = my_params[ExpfactorS].getTrueValue();
    double index2 = my_params[Index2].getTrueValue();
    
    int iparam = -1;
    for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
    }
    
    if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "PowerLawSuperExpCutoff4::derivByParam");
    }

    switch (iparam) 
      {
      case Prefactor:
	return value(xarg) / prefactor * my_params[Prefactor].getScale();
	break;
      case IndexS:
	return value(xarg) * log(x/scale) * my_params[IndexS].getScale();
	break;
      case Scale:
	return value(xarg) * ( -(indexS+expfactorS/index2) / scale + expfactorS/index2*pow(x/scale,index2) / scale) * my_params[Scale].getScale();
	break;
      case ExpfactorS:
	return value(xarg) * ( log(x/scale)/index2 + (1.0-pow(x/scale,index2))/index2/index2 ) * my_params[ExpfactorS].getScale();
	break;
      case Index2:
	return - value(xarg) * expfactorS/index2/index2 * ( 2.0/index2 + log(x/scale) + (log(x/scale)-2.0/index2)*pow(x/scale,index2) ) * my_params[Index2].getScale();
	break;
      default:
	break;
      }
    
    return 0;
    
  }
  
} // namespace Likelihood
