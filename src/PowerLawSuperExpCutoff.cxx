/*
 * @file PowerLawSuperExpCutoff.cxx
 * @brief Implementation for the PowerLawSuperExpCutoff Function class
 * @author Damien Parent
 *
 * $Header$
 */

#include <cmath>

#include <iostream>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/PowerLawSuperExpCutoff.h"

namespace Likelihood {
  
  // initialization function used by constructors
  void PowerLawSuperExpCutoff::init( 
				    double Prefactor,  
				    double Index1,
				    double Scale,
				    double Cutoff,
				    double Index2
				    ) {
    
    // Implement PowerLaw class with five named parameters,
    // "Prefactor", "Index1", "Scale(Energy)", "Cutoff(Energy)" , "Index2"
    
    int nParams = 5;
    setMaxNumParams(nParams);
    
    addParam(std::string("Prefactor"), Prefactor, true);
    addParam(std::string("Index1"), Index1, true);
    addParam(std::string("Scale"), Scale, false);
    addParam(std::string("Cutoff"), Cutoff, true);
    addParam(std::string("Index2"), Index2, true);
    
    // Set FuncType and ArgType for use with CompositeFunction hierarchy.
    m_funcType = Addend;
    m_argType = "dArg";
    m_genericName = "PLSuperExpCutoff";
  } 
  // end of initialization
  
  
  double PowerLawSuperExpCutoff::value(optimizers::Arg &xarg) const {
    double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
    
    // assume a standard ordering for the parameters
    enum ParamTypes {Prefactor, Index1, Scale, Cutoff, Index2};
    
    std::vector<optimizers::Parameter> my_params;
    getParams(my_params);
    
    double prefactor = my_params[Prefactor].getTrueValue();
    double index1 = my_params[Index1].getTrueValue();
    double scale = my_params[Scale].getTrueValue();
    double cutoff = my_params[Cutoff].getValue();
    double index2 = my_params[Index2].getTrueValue();
    
    // De Jager
    return prefactor * pow(x/scale,index1) * exp(-pow(x/cutoff,index2));
  }
  
  
  double PowerLawSuperExpCutoff::derivByParam(optimizers::Arg &xarg,
					      const std::string &paramName) const {
    double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
    
    enum ParamTypes {Prefactor, Index1, Scale, Cutoff, Index2};
    
    std::vector<optimizers::Parameter> my_params;
    getParams(my_params);
    
    double prefactor = my_params[Prefactor].getTrueValue();
    double index1 = my_params[Index1].getTrueValue();
    double scale = my_params[Scale].getTrueValue();
    double cutoff = my_params[Cutoff].getValue();
    double index2 = my_params[Index2].getTrueValue();
    
    int iparam = -1;
    for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
    }
    
    if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "PowerLawSuperExpCutoff::derivByParam");
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
      case Cutoff:
	return value(xarg) * ( pow(x/cutoff,index2) * index2 / cutoff ) * my_params[Cutoff].getScale();
	break;
      case Index2:
	return -value(xarg) * pow(x/cutoff,index2) * log(x/cutoff) * my_params[Index2].getScale();
	break;
      default:
	break;
      }
    
    return 0;
    
  }
  
} // namespace Likelihood
