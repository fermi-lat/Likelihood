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
 * @brief Base class for Science Tools Functions, such as spectral models, 
 * fit statistics, etc..
 *
 * The implementation is based on Hippodraw's FunctionBase class.
 *
 * This class uses the Parameter class.
 *
 * @authors J. Chiang, P. Nolan, T. Burnett 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools/Likelihood/Likelihood/Function.h,v 1.2 2003/02/23 22:28:58 jchiang Exp $
 */

class Function {

public:
    
   Function(){};
   Function(const Function&);

   virtual ~Function(){};

   //! parameter access
   unsigned int getMaxNumParams() const {return m_maxNumParams;}

   //! set Parameter value and free state
   virtual void setParam(const std::string &paramName, 
			 double paramValue, bool isFree)
      {setParameter(paramName, paramValue, isFree);}

   //! set Parameter value, preserving current free state
   virtual void setParam(const std::string &paramName, double paramValue) 
      {setParameter(paramName, paramValue);}

   //! set Parameter using a Parameter object
   virtual void setParam(const Parameter &param) 
      {setParameter(param.getName(), param.getValue(), param.isFree());}

   virtual double getParamValue(const std::string &paramName) const;
   Parameter* getParam(const std::string &paramName);
   
   //! partial derivative wrt param (with default)
   virtual double derivByParam(double, const std::string &) const 
      {return 0.;}

   //! derivatives as a group
   virtual void getDerivs(double x, std::vector<double> &derivs) const
      {fetchDerivs(x, derivs, false);}

   //! derivatives wrt free parameters
   virtual void getFreeDerivs(double x, std::vector<double> &derivs) const
      {fetchDerivs(x, derivs, true);}

   //! parameter access in groups
   unsigned int getNumParams() const {return m_parameter.size();}

   virtual void getParamNames(std::vector<std::string> &names) const
      {fetchParamNames(names, false);}
   virtual void getParamValues(std::vector<double> &values) const
      {fetchParamValues(values, false);}
   virtual void getParams(std::vector<Parameter> &params) const
      {params = m_parameter;}

   void setParamValues(const std::vector<double> &paramVec);
   //! do a bit of name mangling to allow for inheritance of setParamValues
   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   
   //! free parameter access
   unsigned int getNumFreeParams() const;

   virtual void getFreeParamNames(std::vector<std::string> &names) const
      {fetchParamNames(names, true);}
   virtual void getFreeParamValues(std::vector<double> &values) const
      {fetchParamValues(values, true);}
   virtual void getFreeParams(std::vector<Parameter> &) const;

   void setFreeParamValues(const std::vector<double> &paramVec);
   //! note name mangling here too
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   //! integral of function wrt data variable
   virtual double integral(double, double) {return 0.;}

   //! function value at current parameters, pure virtual 
   virtual double value(double) const = 0;
   
   //! allow any subclass to behave like a function object 
   virtual double operator()(double) const = 0;

   //! provide a string identifier
   void setMyName(std::string functionName) {m_functionName = functionName;};
   std::string getMyName() const {return m_functionName;};

protected:

   void setMaxNumParams(int nParams) {m_maxNumParams = nParams;}

   void setParameter(const std::string &paramName, double paramValue, 
		     int isFree = -1);

   //! for subclass constructor use
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
