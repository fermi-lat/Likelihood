/** @file Function.cxx
 * @brief Function class implementation
 *
 * $Header:
 */

#include "../Likelihood/Function.h"

#include <vector>
#include <string>

namespace Likelihood {

//! copy constructor 
Function::Function(const Function &func) {
  m_functionName = func.m_functionName;
  m_maxNumParams = func.m_maxNumParams;
  m_parameter = func.m_parameter;
}

//! resets Parameter values if named Parameter is present,
//! otherwise adds the Parameter if there is room.
void Function::setParam(const std::string paramName, 
                        const double paramValue, 
                        const bool isFree) {

/* check if parameter is already present; 
   if so, set it to the new value and free state */

   for (int i=0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
	 m_parameter[i].setValue(paramValue);
	 m_parameter[i].setFree(isFree);
	 return;
      }
   }

/* otherwise, add it to the list if there's space */
   if (m_parameter.size() < m_maxNumParams) {
      Parameter my_param(paramName, paramValue, isFree);
      m_parameter.push_back(my_param);
   } else {
      std::cerr << "Trying to add a parameter " << paramName 
		<< " equal to " << paramValue
		<< ", " << std::endl 
		<< "   but the parameter list is already full." 
		<< std::endl;
   }
}

void Function::setParam(const std::string paramName, 
			   const double paramValue) {

/* check if parameter is already present; 
   if so, set it to the new value and retain the previous free state */

   for (int i=0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName()) {
	 m_parameter[i].setValue(paramValue);
	 return;
      }
   }

/* otherwise, add it to the list if there's space and assume it's free */
   if (m_parameter.size() < m_maxNumParams) {
      Parameter my_param(paramName, paramValue, true);
      m_parameter.push_back(my_param);
   } else {
      std::cerr << "Trying to add a parameter " << paramName 
		<< " equal to " << paramValue
		<< ", " << std::endl 
		<< "   but the parameter list is already full." 
		<< std::endl;
   }
}

//! search for parameter name, return value or zero if not found
double Function::getParamValue(const std::string paramName) const {

   for (int i = 0; i < m_parameter.size(); i++) {
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
Parameter* Function::getParam(const std::string paramName) {

   for (int i = 0; i < m_parameter.size(); i++) {
      if (paramName == m_parameter[i].getName())
         return &(m_parameter[i]);
   }
   std::cerr << "Function::getParam: "
	     << "Parameter " << paramName << " is not found."
             << std::endl;
   return 0;
}

std::vector<std::string> Function::getParamNames() const {
   std::vector<std::string> myNames;
   for (int i = 0; i < m_parameter.size(); i++)
      myNames.push_back(m_parameter[i].getName());
   return myNames;
}

std::vector<double> Function::getParamValues() const {
   std::vector<double> myValues;
   for (int i = 0; i < m_parameter.size(); i++) 
      myValues.push_back(m_parameter[i].getValue());
   return myValues;
}

void Function::setParamValues(const std::vector<double> paramVec) {
   if (paramVec.size() != m_parameter.size()) {
/* should do some exception handling here */
      std::cerr
         << "The input vector size does not match the number of parameters."
         << std::endl;
   } else {
      for (int i = 0; i < m_parameter.size(); i++) {
	 m_parameter[i].setValue(paramVec[i]);
      }
   }
}

int Function::getNumFreeParams() const {
   int j = 0;
   for (int i = 0; i < m_parameter.size(); i++)
      j += m_parameter[i].isFree();
   return j;
}

std::vector<std::string> Function::getFreeParamNames() const {
   std::vector<std::string> myNames;
   for (int i = 0; i < m_parameter.size(); i++)
      if (m_parameter[i].isFree()) myNames.push_back(m_parameter[i].getName());
   return myNames;
}

std::vector<double> Function::getFreeParamValues() const {
   std::vector<double> myValues;
   for (int i = 0; i < m_parameter.size(); i++) 
      if (m_parameter[i].isFree()) 
	 myValues.push_back(m_parameter[i].getValue());
   return myValues;
}

void Function::setFreeParamValues(const std::vector<double> paramVec) {
   if (paramVec.size() != getNumFreeParams()) {
/* should do some exception handling here */
      std::cerr
         << "The input vector size does not match"
	 << " the number of free parameters."
         << std::endl;
   } else {
      int j = 0;
      for (int i = 0; i < m_parameter.size(); i++) 
	 if (m_parameter[i].isFree()) 
	    m_parameter[i].setValue(paramVec[j++]);
   }
}

std::vector<Parameter> Function::getFreeParams() const {
   std::vector<Parameter> myParams;
   for (int i = 0; i < m_parameter.size(); i++)
      if (m_parameter[i].isFree()) myParams.push_back(m_parameter[i]);
   return myParams;
}

} // namespace Likelihood
