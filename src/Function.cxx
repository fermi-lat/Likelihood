/** @file Function.cxx
 * @brief Function class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Function.cxx,v 1.15 2003/04/25 21:33:40 jchiang Exp $
 */

#include <iostream>

#include "Likelihood/Function.h"
#include "Likelihood/SumFunction.h"
#include "Likelihood/ProductFunction.h"

namespace Likelihood {

#if 0 // disable these until assignment operator is provided
// use of these operators will no doubt leak memory....
SumFunction &Function::operator+(Function &a) {
   SumFunction *sumFunc = new SumFunction(*this, a);
   return *sumFunc;
}

ProductFunction &Function::operator*(Function &a) {
   ProductFunction *productFunc = new ProductFunction(*this, a);
   return *productFunc;
}
#endif
   
void Function::setParameter(const std::string &paramName, 
			    double paramValue,
                            int isFree) {
// check if parameter is present...
   for (unsigned int i=0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         m_parameter[i].setValue(paramValue);
// and update the free state if asked (yes, a bit kludgy...)
         if (isFree > -1) m_parameter[i].setFree(isFree);
         return;
      }
   }
   std::cerr << "Trying to set parameter " << paramName 
             << " equal to " << paramValue << ", \n"
             << "   but it isn't present."
             << std::endl;
}

void Function::addParam(const std::string &paramName,
			double paramValue, 
			bool isFree) {

// check if paramName is already present; if so, complain....
   for (unsigned int i=0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         std::cerr << "This parameter name already exists: "
                   << paramName << ";\n"
                   << "   you can't add another one." << std::endl;
         return;
      }
   }
// if there's room, add this guy onto the vector
   if (m_parameter.size() < m_maxNumParams) {
      Parameter my_param(paramName, paramValue, isFree);
      m_parameter.push_back(my_param);
   } else {
      std::cerr << "Can't add parameter " << paramName << ";\n"
                << "   The parameter list is full at " 
                << m_maxNumParams << "." << std::endl;
   }
}

//! search for parameter name, return value or zero if not found
double Function::getParamValue(const std::string &paramName) const {

   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName())
         return m_parameter[i].getValue();
   }
   std::cerr << "Function::getParamValue: "
             << "Parameter " << paramName << " is not found."
             << std::endl;
   return 0.;
}

//! search for parameter name, return Parameter or NULL pointer 
//! if not found
Parameter* Function::getParam(const std::string &paramName) {
   
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName())
         return &(m_parameter[i]);
   }
   std::cerr << "Function::getParam: "
             << "Parameter " << paramName << " is not found."
             << std::endl;
   return 0;
}

void Function::setParamBounds(const std::string &paramName, double lower,
                              double upper) {
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) { 
         m_parameter[i].setBounds(lower, upper);
//          std::cerr << "setting bounds for "
//                    << m_parameter[i].getName() << std::endl;
      }
   }
}

void Function::setParamScale(const std::string &paramName, double scale) {
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         m_parameter[i].setScale(scale);
//          std::cerr << "setting the scale for "
//                    << m_parameter[i].getName() << std::endl;
      }
   }
}

void Function::setParamTrueValue(const std::string &paramName, 
                                 double paramValue) {
   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
         m_parameter[i].setTrueValue(paramValue);
      }
   }
}

void Function::setParamValues(const std::vector<double> &paramVec) {
   if (paramVec.size() != m_parameter.size()) {
      std::cerr
         << "Function::setParamValues: \n"
         << "The input vector size does not match the number of parameters."
         << std::endl;
   } else {
      std::vector<double>::const_iterator it = paramVec.begin();
      setParamValues_(it);
   }
}
   
std::vector<double>::const_iterator Function::setParamValues_(
   std::vector<double>::const_iterator it) {
   for (unsigned int i = 0; i < m_parameter.size(); i++) 
      m_parameter[i].setValue(*it++);
   return it;
}

void Function::setFreeParamValues(
   const std::vector<double> &paramVec) {
   if (paramVec.size() != getNumFreeParams()) {
      std::cerr
         << "The input vector size does not match"
         << " the number of free parameters."
         << std::endl;
   } else {
      std::vector<double>::const_iterator it = paramVec.begin();
      setFreeParamValues_(it);
   }
}

std::vector<double>::const_iterator Function::setFreeParamValues_(
   std::vector<double>::const_iterator it) {
   for (unsigned int i = 0; i < m_parameter.size(); i++) 
      if (m_parameter[i].isFree()) m_parameter[i].setValue(*it++);
   return it;
}

unsigned int Function::getNumFreeParams() const {
   int j = 0;
   for (unsigned int i = 0; i < m_parameter.size(); i++)
      j += m_parameter[i].isFree();
   return j;
}

void Function::getFreeParams(std::vector<Parameter> &params) const {
   if (!params.empty()) params.clear();
   
   for (unsigned int i = 0; i < m_parameter.size(); i++)
      if (m_parameter[i].isFree()) params.push_back(m_parameter[i]);
}

void Function::setFreeParams(std::vector<Parameter> &params) 
   throw(LikelihoodException) {
   if (params.size() == getNumFreeParams()) {
      int j = 0;
      for (unsigned int i = 0; i < m_parameter.size(); i++) {
         if (m_parameter[i].isFree()) {
            m_parameter[i] = params[j];
            std::cerr << m_parameter[i].getName() << ": " 
                      << m_parameter[i].getValue() << std::endl;
            j++;
         }
      }
   } else {
      throw LikelihoodException
         ("Function::setFreeParams: incompatible number of parameters.");
   }
}      

void Function::fetchParamValues(std::vector<double> &values,
                                bool getFree) const {
   if (!values.empty()) values.clear();

   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (!getFree || m_parameter[i].isFree())
         values.push_back(m_parameter[i].getValue());
   }
}

void Function::fetchParamNames(std::vector<std::string> &names,
                               bool getFree) const {
   if (!names.empty()) names.clear();

   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (!getFree || m_parameter[i].isFree())
         names.push_back(m_parameter[i].getName());
   }
}

void Function::fetchDerivs(Arg &x, std::vector<double> &derivs, 
                           bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   for (unsigned int i = 0; i < m_parameter.size(); i++) {
      if (!getFree || m_parameter[i].isFree())
         derivs.push_back(derivByParam(x, m_parameter[i].getName()));
   }
}

} // namespace Likelihood
