/** @file Function.h
 * @brief Declaration of Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Function.h,v 1.10 2003/03/17 00:53:43 jchiang Exp $
 */

#ifndef Function_h
#define Function_h

#include <iostream>
#include <vector>
#include <string>

#include "Likelihood/Parameter.h"
#include "Likelihood/Arg.h"

namespace Likelihood {

class SumFunction;
class ProductFunction;

/** 
 * @class Function
 *
 * @brief Base class for Science Tools Functions, such as spectral models, 
 * fit statistics, etc..
 *
 * The implementation is based on Hippodraw's FunctionBase class.
 *
 * This class uses the Parameter and Arg classes.
 *
 * @authors J. Chiang, P. Nolan, T. Burnett 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Function.h,v 1.10 2003/03/17 00:53:43 jchiang Exp $
 */

class Function {

public:
    
   Function() {}

// need only member-wise copying
//   Function(const Function&);

   virtual ~Function(){}

   //! provide a string identifier
   void setMyName(std::string functionName) {m_functionName = functionName;}
   std::string getMyName() const {return m_functionName;}

   ///////////////////////
   //! parameter access 
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
   //! parameter access in groups 
   unsigned int getNumParams() const {return m_parameter.size();}

   void setParamValues(const std::vector<double> &paramVec);

   //! do a bit of name mangling to allow for inheritance of setParamValues
   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   
   void setParams(std::vector<Parameter> &params) {
      if (params.size() == m_parameter.size()) 
         m_parameter = params;
   }

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
   
   /////////////////////////////////////
   //! Arg-dependent member functions 
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

   //! clone function, with default
   virtual Function *clone() const {return 0;}

   ////////////////////////////////////////////////////
   //! +, * operator overloading (and memory leaking)
   SumFunction &operator+(Function &);
   ProductFunction &operator*(Function &);

   enum FuncType {None, Addend, Factor};

   FuncType funcType() {return m_funcType;}
   std::string &argType() {return m_argType;}

protected:

   FuncType m_funcType;

   std::string m_argType;

   void setMaxNumParams(int nParams) {m_maxNumParams = nParams;}

   void setParameter(const std::string &paramName, double paramValue, 
                     int isFree = -1);

   //! for subclass constructor use
   void addParam(const std::string &paramName, double paramValue, bool isFree);
   void addParam(const std::string &paramName, double paramValue)
      {addParam(paramName, paramValue, true);}
   void addParam(const Parameter &param) 
      {addParam(param.getName(), param.getValue(), param.isFree());}

   void fetchParamValues(std::vector<double> &values, bool getFree) const;
   void fetchParamNames(std::vector<std::string> &names, bool getFree) const;

   virtual void fetchDerivs(Arg &x ,std::vector<double> &derivs, 
                            bool getFree) const;

   unsigned int m_maxNumParams;

   std::string m_functionName;

   mutable std::vector<Parameter> m_parameter;

};

} // namespace Likelihood

#endif // Function_h
