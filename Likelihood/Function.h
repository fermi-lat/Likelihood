/** @file Function.h
 * @brief Declaration of Function class
 * $Header:
 */

#ifndef Function_h
#define Function_h

#include "Parameter.h"
#include <vector>
#include <string>

namespace Likelihood {

/** 
 * @class Function
 *
 * @brief Base class for Science Tools Functions, i.e., things that are
 * not data containers per se, such as source models, fit statistics, etc..
 * This class uses the Parameter class.
 *
 * @authors J. Chiang, P. Nolan, T. Burnett
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/src/Function.h,v 1.1.1.1 2003/01/30 23:23:03 burnett Exp $
 */

class Function {

public:
    
   Function(){};
   Function(const Function&);

   virtual ~Function(){};

   //! parameter access
   int getMaxNumParams() const {return m_maxNumParams;};
   void setParam(const std::string paramName, const double paramValue,
		 const bool isFree);
   void setParam(const std::string paramName, const double paramValue);
   void setParam(const Parameter param) 
      {setParam(param.getName(), param.getValue(), param.isFree());};
   double getParamValue(const std::string paramName) const;
   Parameter* getParam(const std::string paramName);
   
   //! function value at current parameters, 
   //! pure virtual so each subclass must implement
   virtual double value(const double) const = 0;
   
   //! allow any subclass to behave like a function object 
   virtual double operator()(const double) const = 0;

   //! partial derivative wrt param (virtual, with default)
   virtual double derivByParam(const double, 
			       const std::string paramName) {return 0.;};

   //! derivatives as a group
   virtual std::vector<double> getDerivs(const double) const 
      {return std::vector<double> ();}

   //! derivatives wrt free parameters
   virtual std::vector<double> getFreeDerivs(const double) const 
      {return std::vector<double> ();}

   //! provide a string identifier
   void setMyName(std::string functionName) {m_functionName = functionName;};
   std::string getMyName() const {return m_functionName;};

   //! parameter access in groups
   int getNumParams() const {return m_parameter.size();};
   std::vector<std::string> getParamNames() const;
   std::vector<double> getParamValues() const;
   void setParamValues(const std::vector<double> paramVec);
   std::vector<Parameter> getParams() const {return m_parameter;};
   
   //! free parameter access
   int getNumFreeParams() const;
   std::vector<std::string> getFreeParamNames() const;
   std::vector<double> getFreeParamValues() const;
   void setFreeParamValues(const std::vector<double> paramVec);
   std::vector<Parameter> getFreeParams() const;

protected:

   void setMaxNumParams(const int nParams) {m_maxNumParams = nParams;};

   int m_maxNumParams;

   std::string m_functionName;

   std::vector<Parameter> m_parameter;

};

} // namespace Likelihood

#endif // Function_h
