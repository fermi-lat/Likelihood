/*
 * @file SmoothBrokenPowerLaw.cxx
 * @brief Implementation for the SmoothBrokenPowerLaw Function class
 * @author Benoit Lott
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SmoothBrokenPowerLaw.cxx,v 1.1 2009/06/08 06:05:49 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/SmoothBrokenPowerLaw.h"

namespace Likelihood {
  
  // initialization function used by constructors
  void SmoothBrokenPowerLaw::init( 
				    double Prefactor,  
				    double Index1,
				    double Scale,
				    double Index2,
				    double BreakValue,
				    double Beta
				    ) {
    
    // Implement PowerLaw class with five named parameters,
    // "Prefactor", "Index1", "Scale", "Index2", "BreakValue(Energy)" , "Beta"
    
    int nParams = 6;
    setMaxNumParams(nParams);
    
    addParam(std::string("Prefactor"), Prefactor, true);
    addParam(std::string("Index1"), Index1, true); 
    addParam(std::string("Scale"), Scale, false);
    addParam(std::string("Index2"), Index2, true);
    addParam(std::string("BreakValue"), BreakValue, true);
    addParam(std::string("Beta"), Beta, true);
    
    // Set FuncType and ArgType for use with CompositeFunction hierarchy.
    m_funcType = Addend;
    m_argType = "dArg";
    m_genericName = "SmoothBrokenPowerLaw";
    m_normParName = "Prefactor";
  } 
  // end of initialization
  
  
  double SmoothBrokenPowerLaw::value(optimizers::Arg &xarg) const {
    double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
    
    // assume a standard ordering for the parameters
    enum ParamTypes {Prefactor, Index1, Scale, Index2, BreakValue, Beta};
    
    std::vector<optimizers::Parameter> my_params;
    getParams(my_params);
    
    double prefactor = my_params[Prefactor].getTrueValue();
    double index1 = my_params[Index1].getTrueValue();
    double scale = my_params[Scale].getTrueValue();
    double index2 = my_params[Index2].getTrueValue();
    double breakvalue = my_params[BreakValue].getTrueValue();
    double beta = my_params[Beta].getTrueValue();
    
    // 
    return prefactor * pow(x/scale,index1) * pow(1+ pow(x/breakvalue,(index1-index2)/beta),-beta);
  }
  
  
  double SmoothBrokenPowerLaw::derivByParam(optimizers::Arg &xarg,
					      const std::string &paramName) const {
    double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
    
    enum ParamTypes {Prefactor, Index1, Scale, Index2, BreakValue, Beta};
    
    std::vector<optimizers::Parameter> my_params;
    getParams(my_params);
    
    double prefactor = my_params[Prefactor].getTrueValue();
    double index1 = my_params[Index1].getTrueValue();
    double scale = my_params[Scale].getTrueValue();
    double index2 = my_params[Index2].getTrueValue();
    double breakvalue = my_params[BreakValue].getTrueValue();
    double beta = my_params[Beta].getTrueValue();
    
    int iparam = -1;
    for (unsigned int i = 0; i < my_params.size(); i++) {
      if (paramName == my_params[i].getName()) iparam = i;
    }
    
    if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "SmoothBrokenPowerLaw::derivByParam");
    }
    switch (iparam) 
      {
      case Prefactor:
	return value(xarg) / prefactor * my_params[Prefactor].getScale();
	break;
      case Index1:
	return value(xarg) * (log(x/scale)-log(x/breakvalue)*pow(x/breakvalue,(index1-index2)/beta)/(1+ pow(x/breakvalue,(index1-index2)/beta))) * my_params[Index1].getScale();
	break;
      case Scale:
	return -value(xarg) * index1/ scale * my_params[Scale].getScale();
	break;
      case Index2:
	return value(xarg) *log(x/breakvalue)*pow(x/breakvalue,(index1-index2)/beta)/(1+ pow(x/breakvalue,(index1-index2)/beta))* my_params[Index2].getScale();
	break;
      case BreakValue:
	return value(xarg) *(index1-index2) *pow(x/breakvalue,(index1-index2)/beta)/(1+ pow(x/breakvalue,(index1-index2)/beta))/breakvalue *my_params[BreakValue].getScale();
	break;
      case Beta:
	return  value(xarg) * (-log(1+ pow(x/breakvalue,(index1-index2)/beta))+(index1-index2)*log(x/breakvalue)/beta/(1+ pow(x/breakvalue,-(index1-index2)/beta))) * my_params[Beta].getScale();
	break;
      default:
	break;
      }
    
    return 0;
    
  }
  
} // namespace Likelihood
