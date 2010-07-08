/**
 * @file Composite2.cxx
 * @brief Statistic object that comprises LogLike objects that have a
 * source (e.g. DM source) with common fit parameters intended to be
 * tied together.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Composite2.cxx,v 1.9 2010/05/17 21:17:49 jchiang Exp $
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Likelihood/Composite2.h"

namespace Likelihood {

void Composite2::addComponent(LogLike & like) {
   std::vector<size_t> tiedPars;
   m_components[&like] = tiedPars;
}

void Composite2::tieParameters(const TiedParameter::ParVector_t & pars) {
}

double Composite2::value() const {
   double my_value(0);
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
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
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
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

void Composite2::
setFreeParamValues(const std::vector<double> & values) {
   if (m_components.empty()) {
      throw std::runtime_error("setFreeParamValues: empty composite list");
   }

   size_t j(0);

// Loop over LogLike components and set free, untied parameters
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
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

void Composite2::syncParams() {
   for (ComponentIterator_t it(m_components.begin()); 
        it != m_components.end(); ++it) {
      it->first->syncParams();
   }
}

unsigned int Composite2::getNumFreeParams() const {
   unsigned int npars(0);
// Loop over LogLike components and count free, untied parameters
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
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
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
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
      
   for (size_t k(0); k < tp_derivs.size(); k++) {
      derivs.push_back(tp_derivs.at(k));
   }
}

} // namespace Likleihood
