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
#include "../Likelihood/Arg.h"

namespace Likelihood {

SourceModel::SourceModel(const SourceModel &rhs) : 
   Function(rhs) {
   m_sources = rhs.m_sources;
}

void SourceModel::setParam(const Parameter &param, 
			   const std::string &funcName,
			   const std::string &srcName) {
   for (unsigned int i = 0; i < m_sources.size(); i++) {
      if (srcName == (*m_sources[i]).getName()) {
         Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
	 if (srcFuncs.count(funcName)) {
	    srcFuncs[funcName]->setParam(param.getName(), 
                                         param.getValue(),
                                         param.isFree());
// this seems inefficient, but necessary because of srcFuncs map(?)
            m_syncParams();
	    return;
	 }
      }
   }
   std::cerr << "SourceModel::setParam:  Function " 
	     << funcName << " for source "
	     << srcName << " was not found." << std::endl;
}
 
std::vector<double>::const_iterator SourceModel::setParamValues_(
   std::vector<double>::const_iterator it) { 
   for (unsigned int i = 0; i < m_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
	 it = (*func_it).second->setParamValues_(it);
   }
   m_syncParams();
   return it;
}

std::vector<double>::const_iterator SourceModel::setFreeParamValues_(
   std::vector<double>::const_iterator it) { 
   for (unsigned int i = 0; i < m_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
	 it = (*func_it).second->setFreeParamValues_(it);
   }
   m_syncParams();
   return it;
}

Parameter* SourceModel::getParam(const std::string &paramName,
				 const std::string &funcName,
                                 const std::string &srcName) const {
   for (unsigned int i = 0; i < m_sources.size(); i++) {
      if (srcName == (*m_sources[i]).getName()) {
	 std::vector<Parameter> params;
         Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
         if (srcFuncs.count(funcName)) {    //check for funcName
            srcFuncs[funcName]->getParams(params);
            for (unsigned int j = 0; j < params.size(); j++) 
               if (paramName == params[j].getName())
                  return &(params[j]);
            return 0; //if paramName not found
         }
         return 0; //if funcName not found
      }
   }
   return 0; //if srcName not found
}

void SourceModel::addSource(Source *src) {
// loop over sources to ensure unique names
   for (unsigned int i = 0; i < m_sources.size(); i++) 
      assert((*src).getName() != (*m_sources[i]).getName());

// add this one to the vector
   m_sources.push_back(src);

// add the Parameters to the m_parameter vector 
// (would it be better just to use m_syncParams() here?)
   Source::FuncMap srcFuncs = (*src).getSrcFuncs();
   Source::FuncMap::iterator func_it = srcFuncs.begin();
   for (; func_it != srcFuncs.end(); func_it++) {
      std::vector<Parameter> params;
      (*func_it).second->getParams(params);
      for (unsigned int ip = 0; ip < params.size(); ip++) 
	 m_parameter.push_back(params[ip]);
   }      
}
 
void SourceModel::deleteSource(const std::string &srcName) {
   for (unsigned int i = 0; i < m_sources.size(); i++) {
      if (srcName == (*m_sources[i]).getName()) {
	 m_sources.erase(m_sources.begin() + i, 
                         m_sources.begin() + i + 1);
         m_syncParams();
         return;
      }
   }
   std::cerr << "SourceModel::deleteSource: " 
             << srcName << " was not found." << std::endl;
}

// remake parameter vector from scratch 
void SourceModel::m_syncParams() {
// can we just do m_parameter.clear() here?
   m_parameter.erase(m_parameter.begin(), m_parameter.end());

   for (unsigned int i = 0; i < m_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<Parameter> params;
         (*func_it).second->getParams(params);
	 for (unsigned int ip = 0; ip < params.size(); ip++)
            m_parameter.push_back(params[ip]);
      }
   }
}

void SourceModel::getSrcNames(std::vector<std::string> &names) const {
// ensure names is empty
   if (!names.empty()) names.erase(names.begin(), names.end());

   for (unsigned int i = 0; i < m_sources.size(); i++) {
      names.push_back(m_sources[i]->getName());
   }
}

double SourceModel::evaluate_at(Arg &x) const {
   double my_val = 0.;
   for (unsigned int i = 0; i < m_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
	 my_val += (*func_it).second->value(x);
      }
   }
   return my_val;
}

void SourceModel::fetchDerivs(Arg &x, std::vector<double> &derivs, 
                              bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   for (unsigned int i = 0; i < m_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*m_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<double> my_derivs;
         if (getFree) {
            (*func_it).second->getFreeDerivs(x, my_derivs);
         } else {
            (*func_it).second->getDerivs(x, my_derivs);
         }
         for (unsigned int j = 0; j < my_derivs.size(); j++) 
            derivs.push_back(my_derivs[j]);
      }
   }
}

} // namespace Likelihood
