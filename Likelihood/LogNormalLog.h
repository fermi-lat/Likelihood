/** 
 * @file LogNormalLog.h
 * @brief class for returning the log of a LogNormal function
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LogNormalLog.h,v 1.3 2015/03/19 17:35:06 jchiang Exp $
 */

#ifndef Likelihood_LogNormalLog_h
#define Likelihood_LogNormalLog_h

#include "LogNormal.h"
#include "optimizers/Arg.h"
#include "optimizers/dArg.h"

namespace Likelihood {

/** 
 * @class LogNormalLog
 * Returns the Log of a LogNormal function
 * This is needed in case of logNormal prior,
 * as the prior defition has to include the Log from the likelihood 
 * function definition
 */
    
class LogNormalLog : public LogNormal {

public:

   LogNormalLog(double prefactor=1, double log10_mean=3,
		double log10_sigma=2) {
     LogNormal(prefactor, log10_mean,log10_sigma);
     setGenericName("LogNormalLog");
     //Note: I am not redefining the prefactor here....so it is now additive
     // instead of multiplicative
   }

   virtual Function * clone() const {
      return new LogNormalLog(*this);
   }

   double derivative(const optimizers::Arg & xarg) const {
     double x = dynamic_cast<const optimizers::dArg &>(xarg).getValue();
     enum ParamTypes {Norm, Mean, Sigma};
     double mean = m_parameter[Mean].getTrueValue();
     double sigma = m_parameter[Sigma].getTrueValue();
     return -(1 + (std::log10(x)-mean)/std::log(10)/sigma/sigma)/x;
   }

protected:

   double value(const optimizers::Arg & arg) const {
     return std::log(LogNormal::value(arg));
   }

   double derivByParamImp(const optimizers::Arg & x, 
                          const std::string & paramName) const {
     return LogNormal::derivByParam(x,paramName)/LogNormal::value(x);
   }

};

} // namespace Likelihood

#endif // Likelihood_LogNormalLog_h
