/** @file Function.h
 * @brief Declaration of Function class
 * $Header:
 */

#ifndef Function_h
#define Function_h

#include <vector>
#include <string>

#include "../Likelihood/Parameter.h"
#include "../Likelihood/Arg.h"

namespace Likelihood {

/** 
 * @class Function
 *
 * @brief Base class for Science Tools Functions, such as spectral models, 
 * fit statistics, etc..
 *
 * The implementation is based on Hippodraw's FunctionBase class.
 *
 * This class uses the Parameter classe.
 *
 * @authors J. Chiang, P. Nolan, T. Burnett 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Function.h,v 1.4 2003/02/27 18:47:43 jchiang Exp $
 */

class Function {

public:
    
   Function(){};
   Function(const Function&);

   virtual ~Function(){};

   //! provide a string identifier
   void setMyName(std::string functionName) {m_functionName = functionName;};
   std::string getMyName() const {return m_functionName;};

   ///////////////////////
   //! parameter access //
   ///////////////////////
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

   double getParamValue(const std::string &paramName) const;
   Parameter* getParam(const std::string &paramName);
   
   /////////////////////////////////
   //! parameter access in groups //
   /////////////////////////////////
   unsigned int getNumParams() const {return m_parameter.size();}

   void setParamValues(const std::vector<double> &paramVec);

   //! do a bit of name mangling to allow for inheritance of setParamValues
   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   
   void getParamNames(std::vector<std::string> &names) const
      {fetchParamNames(names, false);}
   void getParamValues(std::vector<double> &values) const
      {fetchParamValues(values, false);}
   void getParams(std::vector<Parameter> &params) const
      {params = m_parameter;}

   //! free parameter access
   unsigned int getNumFreeParams() const;

   void setFreeParamValues(const std::vector<double> &paramVec);

   //! note name mangling here too
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   void getFreeParamNames(std::vector<std::string> &names) const
      {fetchParamNames(names, true);}
   void getFreeParamValues(std::vector<double> &values) const
      {fetchParamValues(values, true);}
   void getFreeParams(std::vector<Parameter> &) const;
   
   /////////////////////////////////////////
   //! Arg-dependent member functions //
   /////////////////////////////////////////
   virtual double value(Arg &) const = 0;
   double operator()(Arg &x) const {return value(x);}
   
   virtual double derivByParam(Arg &x, 
			       const std::string &paramName) const = 0;
   
   virtual void getDerivs(Arg &x, std::vector<double> &derivs) const
      {fetchDerivs(x, derivs, false);}
   
   virtual void getFreeDerivs(Arg &x, std::vector<double> &derivs) const
      {fetchDerivs(x, derivs, true);}

   //! integral of function wrt data variable
   virtual double integral(Arg &, Arg &) const {return 0;}

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

   void fetchDerivs(Arg &x , std::vector<double> &derivs, 
		    bool getFree) const;

   unsigned int m_maxNumParams;

   std::string m_functionName;

   std::vector<Parameter> m_parameter;

};

} // namespace Likelihood

#endif // Function_h
