/** 
 * @file GaussianError.cxx
 * @brief Implementation for the GaussianError Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/GaussianError.cxx,v 1.3 2015/03/03 18:05:37 jchiang Exp $
 */

#include <cmath>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "optimizers/dArg.h"
#include "optimizers/ParameterNotFound.h"

#include "Likelihood/GaussianError.h"

#define SQRT_2_PI 2.5066282746310002
#define LOG_SQRT_2_PI 0.9189385332046727

namespace Likelihood {

  GaussianError::GaussianError(double norm, double mean, double sigma, double offset) 
   : optimizers::Function("GaussianError", 4, "Norm") {
    addParam("Norm", norm, true);
    addParam("Mean", mean, true);
    addParam("Sigma", sigma, true);
    addParam("Offset", offset, true);    
  }

  double GaussianError::value(const optimizers::Arg & xarg) const {
    double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
    double norm = m_parameter[Norm].getTrueValue();
    double mean = m_parameter[Mean].getTrueValue();
    double sigma = m_parameter[Sigma].getTrueValue();
    double offset = m_parameter[Sigma].getTrueValue();
    
    double retVal = -(x - mean)*(x - mean)/(2.*sigma*sigma);
    retVal += offset;
    retVal *= norm;      
    return retVal;
  }
  
  double GaussianError::derivByParamImp(const optimizers::Arg & xarg,
					const std::string & paramName) const {
    double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
    
    int iparam(-1);
    for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
	iparam = i;
	break;
      }
    }
    
    if (iparam == -1) {
      throw optimizers::ParameterNotFound(paramName, getName(),
                                          "GaussianError::derivByParam");
    }
    
    double norm = m_parameter[Norm].getTrueValue();
    double mean = m_parameter[Mean].getTrueValue();
    double sigma = m_parameter[Sigma].getTrueValue();
    double offset = m_parameter[Offset].getTrueValue();
    double s3 = sigma*sigma*sigma;
    
    switch (iparam) {
    case Norm:
      return value(xarg)/norm;
      break;
    case Mean:
      return norm*(x - mean)/(sigma*sigma);
      break;
    case Sigma:
      return norm*(x - mean)*(x - mean)/(6.*s3);
      break;
    case Offset:
      return norm;
      break;
    default:
      break;
    }
    return 0;
  }
  
  double GaussianError::derivative(const optimizers::Arg & xarg) const {
    double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
    double norm = m_parameter[Norm].getTrueValue();
    double mean = m_parameter[Mean].getTrueValue();
    double sigma = m_parameter[Sigma].getTrueValue();
    return norm*(-(x - mean)/(sigma*sigma));
  }
  
} // namespace Likelihood
