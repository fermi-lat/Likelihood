/**
 * @file Composite2.cxx
 * @brief Composite log-likelihood that allows for tying together of
 * arbitrary combinations of parameters.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Composite2.cxx,v 1.7 2011/03/02 23:49:17 jchiang Exp $
 */

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Likelihood/Composite2.h"

namespace Likelihood {

Composite2::~Composite2() throw() {
   try {
      for (std::vector<TiedParameter *>::iterator it(m_tiedPars.begin());
           it != m_tiedPars.end(); ++it) {
         delete *it;
      }
   } catch (...) {
   }
}

void Composite2::addComponent(LogLike & like) {
   std::vector<size_t> tiedPars;
   m_components[&like] = tiedPars;
}

void Composite2::tieParameters(const TiedParameter::ParVector_t & pars) {
   TiedParameter * tiedPar = new TiedParameter();
   for (TiedParameter::ParVectorConstIterator_t it(pars.begin());
        it != pars.end(); ++it) {
      tiedPar->addParam(*it->first, it->second);
      m_components[it->first].push_back(it->second);
   }
   m_tiedPars.push_back(tiedPar);
}

double Composite2::value() const {
   double my_value(0);
   for (ComponentConstIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      my_value += it->first->value();
   }
   return my_value;
}

void Composite2::
getFreeParams(std::vector<optimizers::Parameter> & params) const {
   params.clear();
   if (m_components.empty()) {
      throw std::runtime_error("getFreeParams: empty composite list");
   }
// Loop over LogLike components and gather up free, untied parameters
   for (ComponentConstIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      const std::vector<size_t> & tiedPars(it->second);
      const std::vector<optimizers::Parameter> & pars(it->first->parameters());
      for (size_t i(0); i < pars.size(); i++) {
         if (pars.at(i).isFree() 
             && std::count(tiedPars.begin(), tiedPars.end(), i) == 0) {
            params.push_back(pars.at(i));
         }
      }
   }
// Add the free TiedParameters.  The copy operation will slice the
// TiedParameter attributes, but that's ok in this context.
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars.at(i)->isFree()) {
         params.push_back(*m_tiedPars.at(i));
      }
   }
}

  void Composite2::fetchParamValues(std::vector<double> &values, bool getFree) const {
    if (!values.empty()) values.clear();
    std::vector<optimizers::Parameter> params;
    if(getFree)
      getFreeParams(params);
    else
      getParams(params);
    for(int i=0;i<params.size();i++){
      values.push_back(params.at(i).getValue());
    }
  }

void Composite2::
setFreeParamValues(const std::vector<double> & values) {
   if (m_components.empty()) {
      throw std::runtime_error("setFreeParamValues: empty composite list");
   }
   size_t j(0);
// Loop over LogLike components and set free, untied parameters
   for (ComponentConstIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      const std::vector<size_t> & tiedPars(it->second);
      std::vector<optimizers::Parameter> & pars(it->first->parameters());
      for (size_t i(0); i < pars.size(); i++) {
         if (pars.at(i).isFree() 
             && std::count(tiedPars.begin(), tiedPars.end(), i) == 0) {
            pars.at(i).setValue(values.at(j++));
         }
      }
   }
// Set the free TiedParameters.
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars.at(i)->isFree()) {
         m_tiedPars.at(i)->setValue(values.at(j++));
      }
   }
   syncParams();
}

void Composite2::setErrors(const std::vector<double> & errors) {
   if (m_components.empty()) {
      throw std::runtime_error("setErrors: empty composite list");
   }
   size_t j(0);
// Loop over LogLike components and set free, untied parameters
   for (ComponentConstIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      const std::vector<size_t> & tiedPars(it->second);
      std::vector<optimizers::Parameter> & pars(it->first->parameters());
      for (size_t i(0); i < pars.size(); i++) {
         if (pars.at(i).isFree() 
             && std::count(tiedPars.begin(), tiedPars.end(), i) == 0) {
            pars.at(i).setError(errors.at(j++));
         }
      }
   }
// Set the free TiedParameters.
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars.at(i)->isFree()) {
         m_tiedPars.at(i)->setError(errors.at(j++));
      }
   }
   syncParams();
}

void Composite2::syncParams() {
   m_parameter.clear();
   for (ComponentIterator_t it(m_components.begin()); 
        it != m_components.end(); ++it) {
      const std::vector<size_t> & tiedPars(it->second);
      std::vector<optimizers::Parameter> freePars;
      const std::vector<optimizers::Parameter> & pars(it->first->parameters());
      for (size_t i(0); i < pars.size(); i++) {
         if (pars.at(i).isFree()) {
            freePars.push_back(pars.at(i));
         }
	 //Fill the optimizers::Function m_parameter vector component 
	 //by component, leaving the tied parameters for the end 
         if (std::count(tiedPars.begin(), tiedPars.end(), i) == 0) {
            m_parameter.push_back(pars.at(i));
         }
      }
      it->first->setFreeParams(freePars);
   }
   // Now fill m_parameter with the tied parameters.
   for (size_t i(0); i < m_tiedPars.size(); i++) {
     m_parameter.push_back(*m_tiedPars.at(i));
   }
}

unsigned int Composite2::getNumFreeParams() const {
   unsigned int npars(0);
// Loop over LogLike components and count free, untied parameters
   for (ComponentConstIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      const std::vector<size_t> & tiedPars(it->second);
      const std::vector<optimizers::Parameter> & pars(it->first->parameters());
      for (size_t i(0); i < pars.size(); i++) {
         if (pars.at(i).isFree() 
             && std::count(tiedPars.begin(), tiedPars.end(), i) == 0) {
            npars++;
         }
      }
   }
// Add the free TiedParameters.  The copy operation will slice the
// TiedParameter attributes, but that's ok in this context.
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars.at(i)->isFree()) {
         npars++;
      }
   }
   return npars;
}

void Composite2::getFreeDerivs(std::vector<double> & derivs) const {
   derivs.clear();
   std::vector<double> tp_derivs(m_tiedPars.size(), 0);
   for (ComponentConstIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      const std::vector<size_t> & tiedPars(it->second);
      const std::vector<optimizers::Parameter> & pars(it->first->parameters());
      std::vector<double> freeDerivs;
      it->first->getFreeDerivs(freeDerivs);
      size_t j(0);
      for (size_t i(0); i < pars.size(); i++) {
         if (pars.at(i).isFree()) {
            if (std::count(tiedPars.begin(), tiedPars.end(), i) == 0) {
               // normal free parameter
               derivs.push_back(freeDerivs.at(j));
            } else {
               // tied parameter: sum up contribution from individual params
               std::pair<LogLike *, size_t> item = std::make_pair(it->first, i);
               std::vector<TiedParameter *>::const_iterator tp 
                  = m_tiedPars.begin();
               for (size_t k(0) ; tp != m_tiedPars.end(); ++tp, k++) {
                  if (std::count((*tp)->pars().begin(), 
                                 (*tp)->pars().end(), item) != 0) {
                     tp_derivs.at(k) += freeDerivs.at(j);
                  }
               }
            }
            j++;
         }
      } // pars.at(i)
   } // m_components
   /// Append the derivative sums for the tied parameters.
   for (size_t k(0); k < tp_derivs.size(); k++) {
      if (m_tiedPars.at(k)->isFree()) {
         derivs.push_back(tp_derivs.at(k));
      }
   }
}

int Composite2::findIndex(const LogLike & like, size_t par_index) const {
   int j(0);
/// Loop over components and free, untied parameters, advancing index
/// until the desired par_index is found.
   for (ComponentConstIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      const std::vector<size_t> & tiedPars(it->second);
      const std::vector<optimizers::Parameter> & 
         pars(it->first->parameters());
      if (&like == it->first) {
         /// LogLike object matches, so loop over parameters and if
         /// index matches, return the incremented index value, j.
         for (size_t i(0); i < pars.size(); i++) {
            if (std::count(tiedPars.begin(), tiedPars.end(), i) != 0
                || !pars.at(i).isFree()) {
               continue;
            }
            if (par_index == i) {
               return j;
            }
            j++;
         }
      } else {
         /// LogLike object does not match, so just loop over free,
         /// untied parameters and increment the j index.
         for (size_t i(0); i < pars.size(); i++) {
            if (std::count(tiedPars.begin(), tiedPars.end(), i) != 0
                || !pars.at(i).isFree()) {
               continue;
            }
            j++;
         }
      }
   }
/// Loop over tied parameters.
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars.at(i)->containsIndex(like, par_index)) {
         return j;
      }
      j++;
   }
   return -1;
}

TiedParameter & Composite2::getTiedParam(const LogLike & like, size_t i) {
   std::vector<TiedParameter *>::const_iterator tp(m_tiedPars.begin());
   for (; tp != m_tiedPars.end(); ++tp) {
      if ((*tp)->has_member(like, i)) {
         return *(*tp);
      }
   }
}
  
} // namespace Likleihood
