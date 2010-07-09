/**
 * @file TiedParameter.cxx
 * @brief Implementation for class that handles parameters that are 
 * "tied together" in the Xspec sense.
 *
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/TiedParameter.cxx,v 1.4 2010/07/09 03:57:50 jchiang Exp $
 */

#include <algorithm>

#include "Likelihood/LogLike.h"
#include "Likelihood/TiedParameter.h"

namespace Likelihood {

TiedParameter::TiedParameter() : optimizers::Parameter() {}

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

void TiedParameter::addParam(LogLike & like, size_t i) {
   if (m_pars.empty()) {
      optimizers::Parameter::operator=(like.parameters().at(i));
   }
   if (!has_member(like, i)) {
      m_pars.push_back(std::make_pair(&like, i));
      like.parameters().at(i) = *this;
   }
}

void TiedParameter::removeParam(LogLike & like, size_t i) {
   ParVector_t::iterator it = std::find(m_pars.begin(), m_pars.end(),
                                        std::make_pair(&like, i));
   if (it != m_pars.end()) {
      m_pars.erase(it);
   }
}

bool TiedParameter::has_member(const LogLike & like, size_t i) const {
   for (ParVectorConstIterator_t it = m_pars.begin();
        it != m_pars.end(); ++it) {
      if (it->first == &like && it->second == i) {
         return true;
      }
   }
   return false;
}

void TiedParameter::setName(const std::string name) {
   optimizers::Parameter::setName(name);
   for (ParVectorIterator_t it(m_pars.begin()); it != m_pars.end(); ++it) {
      it->first->parameters().at(it->second).setName(name);
   }
}

void TiedParameter::setValue(double value) {
   optimizers::Parameter::setValue(value);
   for (ParVectorIterator_t it(m_pars.begin()); it != m_pars.end(); ++it) {
      size_t indx = it->second;
      it->first->parameters().at(indx).setValue(value);
   }
}

void TiedParameter::setScale(double scale) {
   optimizers::Parameter::setScale(scale);
   for (ParVectorIterator_t it(m_pars.begin()); it != m_pars.end(); ++it) {
      it->first->parameters().at(it->second).setScale(scale);
   }
}

void TiedParameter::setTrueValue(double trueValue) {
   optimizers::Parameter::setTrueValue(trueValue);
   for (ParVectorIterator_t it(m_pars.begin()); it != m_pars.end(); ++it) {
      it->first->parameters().at(it->second).setTrueValue(trueValue);
   }
}

void TiedParameter::setBounds(double minValue, double maxValue) {
   optimizers::Parameter::setBounds(minValue, maxValue);
   for (ParVectorIterator_t it(m_pars.begin()); it != m_pars.end(); ++it) {
      it->first->parameters().at(it->second).setBounds(minValue, maxValue);
   }
}

void TiedParameter::setFree(bool free) {
   optimizers::Parameter::setFree(free);
   for (ParVectorIterator_t it(m_pars.begin()); it != m_pars.end(); ++it) {
      it->first->parameters().at(it->second).setFree(free);
   }
}

void TiedParameter::setAlwaysFixed(bool flag) {
   optimizers::Parameter::setAlwaysFixed(flag);
   for (ParVectorIterator_t it(m_pars.begin()); it != m_pars.end(); ++it) {
      it->first->parameters().at(it->second).setAlwaysFixed(flag);
   }
}

void TiedParameter::setError(double error) {
   optimizers::Parameter::setError(error);
   for (ParVectorIterator_t it(m_pars.begin()); it != m_pars.end(); ++it) {
      it->first->parameters().at(it->second).setError(error);
   }
}

} //namespace Likelihood
