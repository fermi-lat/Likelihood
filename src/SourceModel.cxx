/** @file SourceModel.cxx
 * @brief SourceModel class implementation
 *
 * $Header:
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "../Likelihood/SourceModel.h"

namespace Likelihood {

SourceModel::SourceModel(const SourceModel &rhs) : Function(rhs) {
   m_functions = rhs.m_functions;
   m_sources = rhs.m_sources;
   m_paramIterators = rhs.m_paramIterators;
}

SourceModel::~SourceModel() {
// is this all good enough?
   setMaxNumParams(0);
   m_functions.erase(m_functions.begin(), m_functions.end());
   m_paramIterators.erase(m_paramIterators.begin(), m_paramIterators.end());
};

void SourceModel::setParam(const Parameter &param, const std::string &fName) {
/* loop over the Parameters associated with the function fName
   and update the matching one */
   
   for (unsigned int i = 0; i < m_functions.size(); i++) {
      if (fName == (*m_functions[i]).getMyName()) {
	 for (std::vector<Parameter>::iterator iter = m_paramIterators[i];
	      iter != m_paramIterators[i] 
		 + (*m_functions[i]).getNumFreeParams();
	      iter++) {
	    if ((*iter).getName() == param.getName()) {
	       (*iter).setValue(param.getValue());
	       (*iter).setFree(param.isFree());
/* update the corresponding Parameter in the named Function */
	       (*m_functions[i]).setParam(param.getName(), param.getValue(),
					  param.isFree());
	       return;
	    }
	 }
	 std::cerr << "SourceModel::setParam:  Parameter " << param.getName()
		   << " was not found." << std::endl;
	 return;
      }
   }
   std::cerr << "SourceModel::setParam:  Function " << fName 
	     << " was not found." << std::endl;
}
 
void SourceModel::setParamValues(const std::vector<double> &paramVec) {
   if (paramVec.size() != m_parameter.size()) {
/* should do some exception handling here */
      std::cerr
         << "The input vector size does not match the number of parameters."
         << std::endl;
   } else {
      for (unsigned int i = 0; i < m_parameter.size(); i++) {
	 m_parameter[i].setValue(paramVec[i]);
      }
   }
   m_syncParams();
}

Parameter* SourceModel::getParam(const std::string &paramName,
				 const std::string &fName) const {
   for (unsigned int i = 0; i < m_functions.size(); i++) {
      if (fName == (*m_functions[i]).getMyName()) {
	 std::vector<Parameter> params;
         (*m_functions[i]).getFreeParams(params);
	 for (unsigned int j = 0; j < params.size(); j++) 
	    if (paramName == params[j].getName())
	       return &(params[j]);
	 return 0; //if param not found
      }
   }
   return 0; //if function not found
}
void SourceModel::addSource(Function *func, const std::string &fName) {
   if (getMaxNumParams() < getNumParams() + (*func).getNumFreeParams())
      setMaxNumParams(getNumParams() + (*func).getNumParams());

//! loop over functions to ensure unique names

   for (unsigned int i = 0; i < m_functions.size(); i++) 
      assert(fName != (*m_functions[i]).getMyName());

//! set the function's name (overwriting the existing name in *func)
//! and add to the vector of Functions

   (*func).setMyName(fName);
   m_functions.push_back(func);

//! add the free ones to the Parameter vector
   std::vector<Parameter> params;
   (*func).getFreeParams(params);
   for (unsigned int i = 0; i < params.size() ; i++) {
      m_parameter.push_back(params[i]);

//! store the pointer to first Parameter to be added to the list
      if (i == 0) m_paramIterators.push_back(m_parameter.end()-1);
   }
}

void SourceModel::deleteSource(const std::string &fName) {
   
   for (unsigned int i = 0; i < m_functions.size(); i++) {
      if (fName == (*m_functions[i]).getMyName()) {
	 m_functions.erase(m_functions.begin() + i, 
			   m_functions.begin() + i + 1);
	 break;
      }
   }

// remake parameter vector from scratch
   m_parameter.erase(m_parameter.begin(), m_parameter.end());
   m_paramIterators.erase(m_paramIterators.begin(), m_paramIterators.end());

   for (unsigned int i = 0; i < m_functions.size(); i++) {
      std::vector<Parameter> params;
      (*m_functions[i]).getFreeParams(params);
      for (unsigned int j = 0; j < params.size() ; j++) {
	 m_parameter.push_back(params[j]);
	 if (j == 0) m_paramIterators.push_back(m_parameter.end()-1);
      }
   }
}

Function * SourceModel::getFunc(const std::string &fName) {
   for (unsigned int i = 0; i < m_functions.size(); i++) 
      if (fName == (*m_functions[i]).getMyName()) 
	 return m_functions[i];
   std::cerr << "SourceModel::getFunc: Function " 
	     << fName << " was not found.";
   return 0;
}

void SourceModel::getSrcNames(std::vector<std::string> &names) const {
// ensure names is empty
   if (!names.empty()) names.erase(names.begin(), names.end());

   for (unsigned int i = 0; i < m_functions.size(); i++) {
      names.push_back(m_functions[i]->getMyName());
   }
}

void SourceModel::getNumSrcParams(std::vector<int> &numParams) const {
   if (!numParams.empty()) numParams.erase(numParams.begin(), numParams.end());

   for (unsigned int i = 0; i < m_functions.size(); i++) {
      numParams.push_back(m_functions[i]->getNumFreeParams());
   }
}

double SourceModel::evaluate_at(double x) const {
   double my_val = 0.;
   for (unsigned int i = 0; i < m_functions.size(); i++) 
      my_val += (*m_functions[i])(x);
   return my_val;
}

void SourceModel::getDerivs(double x, std::vector<double> &derivs) const {
// ensure that derivs is empty
   if (!derivs.empty()) derivs.erase(derivs.begin(), derivs.end());

   for (unsigned int i = 0; i < m_functions.size(); i++) {
      std::vector<double> my_freeDerivs;
      (*m_functions[i]).getFreeDerivs(x, my_freeDerivs);
      for (unsigned int j = 0; j < my_freeDerivs.size(); j++)
	 derivs.push_back(my_freeDerivs[j]);
   }
}

void SourceModel::m_syncParams() {
   int k = 0;
   for (unsigned int i = 0; i < m_functions.size(); i++) {
      int n_params = (*m_functions[i]).getNumFreeParams();
      for (int j = 0; j < n_params; j++) {
	 (*m_functions[i]).setParam(m_parameter[k].getName(), 
				    m_parameter[k].getValue(),
				    m_parameter[k].isFree());
	 k++;
      }
   }
}

} // namespace Likelihood
