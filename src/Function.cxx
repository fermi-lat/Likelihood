/** @file Function.cxx
 * @brief Function class implementation
 *
 * $Header:
 */

#include "../Likelihood/Function.h"

namespace Likelihood {

//! copy constructor
Function::Function(const Function &func) {
  m_functionName = func.m_functionName;
  m_maxNumParams = func.m_maxNumParams;
  m_parameter = func.m_parameter;
}

void Function::setParameter(const std::string &paramName, 
				     double paramValue,
				     int isFree = -1) {
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

//! add a Parameter
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

void Function::setParamValues(const std::vector<double> &paramVec) {
   if (paramVec.size() != m_parameter.size()) {
      std::cerr
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
