/** 
 * @file Function.h
 * @brief Declaration of Function and ParameterNotFound classes
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Function.h,v 1.23 2003/05/29 20:10:30 jchiang Exp $
 */

#ifdef _MSC_VER
#pragma warning(disable:4290)
#endif

#ifndef Function_h
#define Function_h

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "Likelihood/Parameter.h"
#include "Likelihood/LikelihoodException.h"

namespace Likelihood {

class Arg;
class ParameterNotFound;

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Function.h,v 1.23 2003/05/29 20:10:30 jchiang Exp $
 */

class Function {

public:
    
   Function() {}

// need only member-wise copying
//   Function(const Function&);

   virtual ~Function() {}

   //! provide a string identifier
   void setName(std::string functionName) {m_functionName = functionName;}
   std::string getName() const {return m_functionName;}

   ///////////////////////
   //! parameter access 
   //! set Parameter value and free state
   virtual void setParam(const std::string &paramName, 
                         double paramValue, bool isFree)
      throw(ParameterNotFound) {setParameter(paramName, paramValue, isFree);}

   //! set Parameter value, preserving current free state
   virtual void setParam(const std::string &paramName, double paramValue) 
      throw(ParameterNotFound) {setParameter(paramName, paramValue);}

   //! set Parameter using a Parameter object
   virtual void setParam(const Parameter &param) throw(ParameterNotFound);

   virtual double getParamValue(const std::string &paramName) const
      throw(ParameterNotFound);

   virtual Parameter getParam(const std::string &paramName) const
      throw(ParameterNotFound);
   
   virtual void setParamBounds(const std::string &paramName, double lower,
                               double upper) throw(ParameterNotFound);
   
   virtual void setParamScale(const std::string &paramName, double scale)
      throw(ParameterNotFound);

   virtual void setParamTrueValue(const std::string &paramName, 
                                  double paramValue) throw(ParameterNotFound);

   /////////////////////////////////
   //! parameter access in groups 
   unsigned int getNumParams() const {return m_parameter.size();}

   void setParamValues(const std::vector<double> &paramVec)
      throw(LikelihoodException);

   //! do a bit of name mangling to allow for inheritance of setParamValues
   virtual std::vector<double>::const_iterator setParamValues_(
      std::vector<double>::const_iterator);
   
   virtual void setParams(std::vector<Parameter> &params) 
      throw(LikelihoodException) {
      if (params.size() == m_parameter.size()) {
         m_parameter = params;
      } else {
         throw LikelihoodException
            ("Function::setParams: incompatible number of parameters.");
      }
   }

   void getParamNames(std::vector<std::string> &names) const
      {fetchParamNames(names, false);}
   void getParamValues(std::vector<double> &values) const
      {fetchParamValues(values, false);}
   void getParams(std::vector<Parameter> &params) const
      {params = m_parameter;}

   //! free parameter access
   unsigned int getNumFreeParams() const;

   void setFreeParamValues(const std::vector<double> &paramVec)
      throw(LikelihoodException);

   //! note name mangling here too
   virtual std::vector<double>::const_iterator setFreeParamValues_(
      std::vector<double>::const_iterator);

   void getFreeParamNames(std::vector<std::string> &names) const
      {fetchParamNames(names, true);}
   void getFreeParamValues(std::vector<double> &values) const
      {fetchParamValues(values, true);}
   void getFreeParams(std::vector<Parameter> &) const;

   virtual void setFreeParams(std::vector<Parameter> &) 
      throw(LikelihoodException);
   
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

   enum FuncType {None, Addend, Factor};

   FuncType funcType() {return m_funcType;}
   std::string &argType() {return m_argType;}

protected:

   FuncType m_funcType;

   std::string m_argType;

   unsigned int m_maxNumParams;

   std::string m_functionName;

   mutable std::vector<Parameter> m_parameter;

   void setMaxNumParams(int nParams) {m_maxNumParams = nParams;}

   void setParameter(const std::string &paramName, double paramValue, 
                     int isFree = -1) throw(ParameterNotFound);

   //! for subclass constructor use
   void addParam(const std::string &paramName, 
                 double paramValue, bool isFree) throw(LikelihoodException);

   void addParam(const std::string &paramName, double paramValue)
      throw(LikelihoodException) {addParam(paramName, paramValue, true);}

   void addParam(const Parameter &param) throw(LikelihoodException);

   void fetchParamValues(std::vector<double> &values, bool getFree) const;
   void fetchParamNames(std::vector<std::string> &names, bool getFree) const;

   virtual void fetchDerivs(Arg &x ,std::vector<double> &derivs, 
                            bool getFree) const;
};

/**
 * @class ParameterNotFound
 * @brief Returns a standard error message for Parameters 
 * looked for but not found in the desired Function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Function.h,v 1.23 2003/05/29 20:10:30 jchiang Exp $
 */
class ParameterNotFound : public LikelihoodException {

public:
   ParameterNotFound(const std::string &paramName, 
                     const std::string &funcName,
                     const std::string &routineName) {
      std::ostringstream errorMessage;
      errorMessage << "Function::" << routineName << ": \n"
                   << "A Parameter named " << paramName
                   << " is not a Parameter of Function "
                   << funcName << "\n";
      m_what = errorMessage.str();
      m_code = 0;
   }
};   

} // namespace Likelihood

#endif // Function_h
