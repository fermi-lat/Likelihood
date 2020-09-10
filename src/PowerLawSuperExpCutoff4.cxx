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

    double y = log(x/scale);
    double y2 = index2*y;
    if(fabs(y2)<1e-2)
      return prefactor * exp( y*(indexS-expfactorS*y/2*(1+y2/3*(1+y2/4))) );
    return prefactor * exp( std::min((indexS+expfactorS/index2) * y + expfactorS/index2/index2 * (1-exp(y2) ), max_expfactor) );
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

    double y = log(x/scale);
    double y2 = index2*y;

    switch (iparam) 
      {
      case Prefactor:
	return value(xarg) / prefactor * my_params[Prefactor].getScale();
	break;
      case IndexS:
	return value(xarg) * y * my_params[IndexS].getScale();
	break;
      case Scale:
        if(fabs(y2)<1e-2)
          return value(xarg) / scale * ( -indexS+expfactorS*y*(1 + y2/2 * (1 + y2/3))) * my_params[Scale].getScale();
	return value(xarg) / scale * ( -(indexS+expfactorS/index2) + expfactorS/index2*exp(y2) ) * my_params[Scale].getScale();
	break;
      case ExpfactorS:
        if(fabs(y2)<1e-2)
          return value(xarg) * (-y*y/2 * (1 + y2/3 * (1 + y2/4))) * my_params[ExpfactorS].getScale();
        return value(xarg) / index2 *( y + (1-exp(y2))/index2 ) * my_params[ExpfactorS].getScale();
	break;
      case Index2:
        if(fabs(y2)<1e-2)
          return - value(xarg) * expfactorS*y*y*y/6 * (1 + y2/2) * my_params[Index2].getScale();
        return - value(xarg) * expfactorS/index2/index2 * ( 2.0/index2 + y + (y-2.0/index2)*exp(y2) ) * my_params[Index2].getScale();
	break;
      default:
	break;
      }
    
    return 0;
    
  }
  
} // namespace Likelihood
