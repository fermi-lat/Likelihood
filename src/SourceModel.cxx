/** 
 * @file SourceModel.cxx
 * @brief SourceModel class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceModel.cxx,v 1.18 2003/06/10 23:58:52 jchiang Exp $
 */

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <cassert>
#include "Likelihood/SourceModel.h"
#include "Likelihood/Arg.h"

namespace Likelihood {

int SourceModel::s_refCount = 0;
std::vector<Source *> SourceModel::s_sources;

SourceModel::~SourceModel() {
   s_refCount--;
   if (s_refCount == 0) {
      std::vector<Source *>::iterator it = s_sources.begin();
      for (; it != s_sources.end(); it++)
         delete (*it);
   }
   s_sources.clear();
}

void SourceModel::setParam(const Parameter &param, 
                           const std::string &funcName,
                           const std::string &srcName) 
   throw(LikelihoodException) {
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (srcName == (*s_sources[i]).getName()) {
         Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
         if (srcFuncs.count(funcName)) {
            srcFuncs[funcName]->setParam(param.getName(), 
                                         param.getValue(),
                                         param.isFree());
// this seems inefficient, but necessary because of srcFuncs map(?)
            syncParams();
            return;
         }
      }
   }
   std::ostringstream errorMessage;
   errorMessage << "SourceModel::setParam:  Function " 
                << funcName << " for Source "
                << srcName << " was not found.\n";
   throw LikelihoodException(errorMessage.str());
}
 
std::vector<double>::const_iterator SourceModel::setParamValues_(
   std::vector<double>::const_iterator it) { 
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
         it = (*func_it).second->setParamValues_(it);
   }
   syncParams();
   return it;
}

std::vector<double>::const_iterator SourceModel::setFreeParamValues_(
   std::vector<double>::const_iterator it) { 
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) 
         it = (*func_it).second->setFreeParamValues_(it);
   }
   syncParams();
   return it;
}

Parameter SourceModel::getParam(const std::string &paramName,
                                const std::string &funcName,
                                const std::string &srcName) const 
   throw(LikelihoodException, ParameterNotFound) {
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (srcName == (*s_sources[i]).getName()) {
         std::vector<Parameter> params;
         Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
         if (srcFuncs.count(funcName)) {    //check for funcName
            srcFuncs[funcName]->getParams(params);
            for (unsigned int j = 0; j < params.size(); j++) {
               if (paramName == params[j].getName()) {
                  return params[j];
               }
            }
            throw ParameterNotFound(paramName, funcName, 
                                    "SourceModel::getParam");
         }
         std::ostringstream errorMessage;
         errorMessage << "SourceModel::getParam:\n"
                      << "Function " << funcName 
                      << " was not found in Source " 
                      << srcName << "\n";
         throw LikelihoodException(errorMessage.str());
      }
   }
   std::ostringstream errorMessage;
   errorMessage << "SourceModel::getParam: "
                << "Source " << srcName 
                << " was not found.\n";
   throw LikelihoodException(errorMessage.str());
}

void SourceModel::setParamBounds(const std::string &paramName,
                                 const std::string &funcName,
                                 const std::string &srcName,
                                 double lower, double upper)
   throw(ParameterNotFound, Parameter::OutOfBounds) {
   Parameter my_param = getParam(paramName, funcName, srcName);
   my_param.setBounds(lower, upper);
   setParam(my_param, funcName, srcName);
   syncParams();
}

void SourceModel::setParamScale(const std::string &paramName,
                                const std::string &funcName,
                                const std::string &srcName,
                                double scale) throw(ParameterNotFound) {
   Parameter my_param = getParam(paramName, funcName, srcName);
   my_param.setScale(scale);
   setParam(my_param, funcName, srcName);
   syncParams();
}

void SourceModel::setParamTrueValue(const std::string &paramName,
                                    const std::string &funcName,
                                    const std::string &srcName,
                                    double paramValue)
   throw(ParameterNotFound, Parameter::OutOfBounds) {
   Parameter my_param = getParam(paramName, funcName, srcName);
   my_param.setTrueValue(paramValue);
   setParam(my_param, funcName, srcName);
   syncParams();
}

void SourceModel::setParams_(std::vector<Parameter> &params, bool setFree) 
   throw(LikelihoodException, ParameterNotFound) {
// ensure the number of Parameters matches
   unsigned int numParams;
   if (setFree) {
      numParams = getNumFreeParams();
   } else {
      numParams = getNumParams();
   }
   if (params.size() != numParams) {
      std::string errorMessage = std::string("SourceModel::setParams:\n") 
         + std::string("Inconsistent number of Parameters.");
      throw LikelihoodException(errorMessage);
   }
// assume ordering of Parameters in params matches that given by the
// ordering of the Sources and their Functions
   int k = 0;  // params' index
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         unsigned int numParams;
         if (setFree) {
            numParams = func_it->second->getNumFreeParams();
         } else { 
            numParams = func_it->second->getNumParams();
         }
         for (unsigned int j = 0; j < numParams; j++) {
            func_it->second->setParam(params[k]);
            k++;
         }
      }
   }
   syncParams();
}

void SourceModel::addSource(Source *src) {
// loop over sources to ensure unique names
   for (unsigned int i = 0; i < s_sources.size(); i++) 
      assert((*src).getName() != (*s_sources[i]).getName());

// add a clone of this Source to the vector
   s_sources.push_back(src->clone());

// add the Parameters to the m_parameter vector 
// (would it be better just to use syncParams() here?)
   Source::FuncMap srcFuncs = (*src).getSrcFuncs();
   Source::FuncMap::iterator func_it = srcFuncs.begin();
   for (; func_it != srcFuncs.end(); func_it++) {
      std::vector<Parameter> params;
      (*func_it).second->getParams(params);
      for (unsigned int ip = 0; ip < params.size(); ip++) 
         m_parameter.push_back(params[ip]);
   }      
}
 
void SourceModel::deleteSource(const std::string &srcName) 
   throw(LikelihoodException) {
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (srcName == (*s_sources[i]).getName()) {
         s_sources.erase(s_sources.begin() + i, 
                         s_sources.begin() + i + 1);
         syncParams();
         return;
      }
   }
   std::ostringstream errorMessage;
   errorMessage << "SourceModel::deleteSource: " 
                << srcName << " was not found.\n";
   throw LikelihoodException(errorMessage.str());
}

void SourceModel::deleteAllSources() {
   std::vector<std::string> srcNames;
   getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++)
      deleteSource(srcNames[i]);
   m_parameter.clear();
}

Source * SourceModel::getSource(const std::string &srcName) {
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (srcName == s_sources[i]->getName())
         return s_sources[i];
   }
   return 0;
}

void SourceModel::getSrcNames(std::vector<std::string> &names) const {
// ensure names is empty
   if (!names.empty()) names.erase(names.begin(), names.end());

   for (unsigned int i = 0; i < s_sources.size(); i++) {
      names.push_back(s_sources[i]->getName());
   }
}

double SourceModel::evaluate_at(Arg &x) const {
   double my_val = 0.;
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         my_val += (*func_it).second->value(x);
      }
   }
   return my_val;
}

// remake parameter vector from scratch 
void SourceModel::syncParams() {
   m_parameter.clear();

   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
      Source::FuncMap::iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<Parameter> params;
         (*func_it).second->getParams(params);
         for (unsigned int ip = 0; ip < params.size(); ip++)
            m_parameter.push_back(params[ip]);
      }
   }
}

void SourceModel::fetchDerivs(Arg &x, std::vector<double> &derivs, 
                              bool getFree) const {
   if (!derivs.empty()) derivs.clear();

   for (unsigned int i = 0; i < s_sources.size(); i++) {
      Source::FuncMap srcFuncs = (*s_sources[i]).getSrcFuncs();
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
