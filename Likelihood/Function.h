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
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/Likelihood/Function.h,v 1.1 2003/02/19 01:34:33 jchiang Exp $
 */

class Function {

public:
    
   Function(){};
   Function(const Function&);

   virtual ~Function(){};

   //! parameter access
   unsigned int getMaxNumParams() const {return m_maxNumParams;};
   void setParam(const std::string &paramName, double paramValue, bool isFree);
   void setParam(const std::string &paramName, double paramValue);
   void setParam(const Parameter &param) 
      {setParam(param.getName(), param.getValue(), param.isFree());};
   double getParamValue(const std::string &paramName) const;
   Parameter* getParam(const std::string &paramName);
   
   //! function value at current parameters, 
   //! pure virtual so each subclass must implement
   virtual double value(double) const = 0;
   
   //! allow any subclass to behave like a function object 
   virtual double operator()(double) const = 0;

   //! partial derivative wrt param (virtual, with default)
   virtual double derivByParam(double, const std::string &) const {return 0.;};

   //! derivatives as a group
   virtual void getDerivs(double x, std::vector<double> &derivs) const
      {fetchDerivs(x, derivs, false);}

   //! derivatives wrt free parameters
   virtual void getFreeDerivs(double x, std::vector<double> &derivs) const
      {fetchDerivs(x, derivs, true);}

   //! integral of function wrt data variable
   virtual double integral(double, double) {return 0.;};

   //! provide a string identifier
   void setMyName(std::string functionName) {m_functionName = functionName;};
   std::string getMyName() const {return m_functionName;};

   //! parameter access in groups
   unsigned int getNumParams() const {return m_parameter.size();};
   void getParamNames(std::vector<std::string> &names) const
      {fetchParamNames(names, false);}
   void getParamValues(std::vector<double> &values) const
      {fetchParamValues(values, false);}
   void setParamValues(const std::vector<double> &paramVec);
   void getParams(std::vector<Parameter> &params) const
      {params = m_parameter;}
   
   //! free parameter access
   unsigned int getNumFreeParams() const;
   void getFreeParamNames(std::vector<std::string> &names) const
      {fetchParamNames(names, true);}
   void getFreeParamValues(std::vector<double> &values) const
      {fetchParamValues(values, true);}
   void setFreeParamValues(const std::vector<double> &paramVec);
   void getFreeParams(std::vector<Parameter> &) const;

protected:

   void setMaxNumParams(int nParams) {m_maxNumParams = nParams;};

   void addParam(const std::string &paramName, double paramValue, bool isFree);
   void addParam(const std::string &paramName, double paramValue)
      {addParam(paramName, paramValue, true);};
   void addParam(const Parameter &param) 
      {addParam(param.getName(), param.getValue(), param.isFree());};

   void fetchParamValues(std::vector<double> &values, bool getFree) const;
   void fetchParamNames(std::vector<std::string> &names, bool getFree) const;
   void fetchDerivs(double x, std::vector<double> &derivs, bool getFree) const;

   unsigned int m_maxNumParams;

   std::string m_functionName;

   std::vector<Parameter> m_parameter;

};

} // namespace Likelihood

#endif // Function_h
