/**
 * @file TiedParameter.cxx
 * @brief Implementation for class that handles parameters that are 
 * "tied together" in the Xspec sense.
 *
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/TiedParameter.cxx,v 1.1 2010/07/08 01:09:39 jchiang Exp $
 */

#include <algorithm>

#include "Likelihood/TiedParameter.h"

namespace Likelihood {

TiedParameter::TiedParameter(const std::string & name, double value,
                             double minValue, double maxValue,
                             bool isFree, double error) 
   : optimizers::Parameter(name, value, minValue, maxValue, isFree, error) {}

TiedParameter::TiedParameter(const TiedParameter & other) 
   : optimizers::Parameter(other), m_pars(other.m_pars) {}

TiedParameter::~TiedParameter() throw() {
   // nothing to delete.
}

TiedParameter & TiedParameter::operator=(const TiedParameter & rhs) {
   if (this != &rhs) {
      optimizers::Parameter::operator=(rhs);
      m_pars = rhs.m_pars;
   }
   return *this;
}

void TiedParameter::addParam(optimizers::Parameter & par) {
   if (!has_member(par)) {
      m_pars.push_back(&par);
   }
}

void TiedParameter::removeParam(optimizers::Parameter & par) {
   ParVector_t::iterator it = std::find(m_pars.begin(), m_pars.end(), &par);
   if (it != m_pars.end()) {
      m_pars.erase(it);
   }
}

bool TiedParameter::has_member(const optimizers::Parameter & par) const {
   if (std::count(m_pars.begin(), m_pars.end(), &par) == 0) {
      return false;
   }
   return true;
}

void TiedParameter::setName(const std::string name) {
   optimizers::Parameter::setName(name);
   for (size_t i(0); i < m_pars.size(); i++) {
      m_pars.at(i)->setName(name);
   }
}

void TiedParameter::setValue(double value) {
   optimizers::Parameter::setValue(value);
   for (size_t i(0); i < m_pars.size(); i++) {
      m_pars.at(i)->setValue(value);
   }
}

void TiedParameter::setScale(double scale) {
   optimizers::Parameter::setScale(scale);
   for (size_t i(0); i < m_pars.size(); i++) {
      m_pars.at(i)->setScale(scale);
   }
}

void TiedParameter::setTrueValue(double trueValue) {
   optimizers::Parameter::setTrueValue(trueValue);
   for (size_t i(0); i < m_pars.size(); i++) {
      m_pars.at(i)->setTrueValue(trueValue);
   }
}

void TiedParameter::setBounds(double minValue, double maxValue) {
   optimizers::Parameter::setBounds(minValue, maxValue);
   for (size_t i(0); i < m_pars.size(); i++) {
      m_pars.at(i)->setBounds(minValue, maxValue);
   }
}

void TiedParameter::setFree(bool free) {
   optimizers::Parameter::setFree(free);
   for (size_t i(0); i < m_pars.size(); i++) {
      m_pars.at(i)->setFree(free);
   }
}

void TiedParameter::setAlwaysFixed(bool flag) {
   optimizers::Parameter::setAlwaysFixed(flag);
   for (size_t i(0); i < m_pars.size(); i++) {
      m_pars.at(i)->setAlwaysFixed(flag);
   }
}

void TiedParameter::setError(double error) {
   optimizers::Parameter::setError(error);
   for (size_t i(0); i < m_pars.size(); i++) {
      m_pars.at(i)->setError(error);
   }
}

} //namespace Likelihood
