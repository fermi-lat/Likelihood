/**
 * @file SummedLikelihood.cxx
 * @brief Statistic object that comprises LogLike objects that have
 * identical source models.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/SummedLikelihood.cxx,v 1.4 2014/02/07 06:37:11 jchiang Exp $
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Likelihood/SummedLikelihood.h"

namespace Likelihood {

SummedLikelihood::~SummedLikelihood() throw() {
   try {
      for (std::vector<TiedParameter *>::iterator it(m_tiedPars.begin());
           it != m_tiedPars.end(); ++it) {
         delete *it;
      }
   } catch (...) {
   }
}

void SummedLikelihood::addComponent(LogLike & component) {
   m_components.push_back(&component);
   if (m_masterComponent == 0) {
      m_masterComponent = &component;
   }
}

double SummedLikelihood::value() const {
   double my_value(0);
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
      my_value += (*it)->value();
   }
   return my_value;
}

void SummedLikelihood::
getFreeParams(std::vector<optimizers::Parameter> & params) const {
   if (m_components.empty()) {
      throw std::runtime_error("SummedLikelihood::getFreeParams: "
                               "empty component list");
   }
   // All components have the same parameters, so we just need to get
   // the ones from the master component.
   const std::vector<optimizers::Parameter> 
      & pars(m_masterComponent->parameters());
   /// Gather free, untied parameters first.
   for (size_t par_index(0); par_index < pars.size(); par_index++) {
      if (pars[par_index].isFree() &&
          m_tiedIndices.find(par_index) == m_tiedIndices.end()) {
         params.push_back(pars[par_index]);
      }
   }
   // Now get the free TiedParameters.
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars[i]->isFree()) {
         params.push_back(*m_tiedPars[i]);
      }
   }
}

void SummedLikelihood::
setFreeParamValues(const std::vector<double> & values) {
   if (m_components.empty()) {
      throw std::runtime_error("SummedLikelihood::setFreeParamValues: "
                               "empty component list");
   }
   size_t j(0);  // This is the index over values.
   std::vector<optimizers::Parameter> & pars(m_masterComponent->parameters());
   for (size_t par_index(0); par_index < pars.size(); par_index++) {
      if (pars[par_index].isFree() &&
          m_tiedIndices.find(par_index) == m_tiedIndices.end()) {
         pars[par_index].setValue(values.at(j++));
      }
   }
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars[i]->isFree()) {
         m_tiedPars[i]->setValue(values.at(j++));
      }
   }
   // Set the parameters for all of the components to be the same.
   for (ComponentIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      // if (*it == m_masterComponent) { 
      //    // The parameters for the master should already have been set.
      //    continue;
      // }
      (*it)->setParams(pars);
   }
   syncParams();
}

void SummedLikelihood::syncParams() {
   for (ComponentIterator_t it(m_components.begin()); 
        it != m_components.end(); ++it) {
      (*it)->syncParams();
   }
}

unsigned int SummedLikelihood::getNumFreeParams() const {
   unsigned int npars(0);
   const std::vector<optimizers::Parameter> & 
      pars(m_masterComponent->parameters());
   for (size_t par_index(0); par_index < pars.size(); par_index++) {
      if (pars[par_index].isFree() &&
          m_tiedIndices.find(par_index) == m_tiedIndices.end()) {
         npars++;
      }
   }
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars[i]->isFree()) {
         npars++;
      }
   }
   return npars;
}

void SummedLikelihood::fetchParamValues(std::vector<double> & values,
                                        bool getFree) const {
   if (getFree) {
      values.clear();
      std::vector<optimizers::Parameter> pars;
      getFreeParams(pars);
      for (size_t i(0); i < pars.size(); i++) {
         values.push_back(pars[i].getValue());
      }
   } else {
      m_masterComponent->getParamValues(values);
   }
}

void SummedLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
   // The free derivatives are the same for all components, so the
   // composite values are just the free derivatives from any single
   // component multiplied by the total number of components.
   //
   // However, we need to extract these in an order that accounts for
   // the handling of the tied parameters.
   size_t numComponents(m_components.size());

   // Get the all of the free derivatives as viewed from the master
   // component.
   m_masterComponent->getFreeDerivs(derivs);

   // Just multiply by number of components if no tied parameter
   // handling is required.
   if (m_tiedIndices.size() == 0) {
      for (size_t i(0); i < derivs.size(); i++) {
         derivs[i] *= numComponents;
      }
      return;
   }
   // Loop over all parameters and use findIndex(par_index) to
   // determine the Minos indices and their order, accounting for the
   // tied parameters.
   std::vector<double> master_derivs(derivs);
   const std::vector<optimizers::Parameter> & 
      pars(m_masterComponent->parameters());
   size_t master_index(0);
   for (size_t par_index(0); par_index < pars.size(); par_index++) {
      int minos_index(findIndex(par_index));
      if (minos_index > -1) {
         derivs.at(minos_index) = numComponents*master_derivs[master_index++];
      }
   }
}

void SummedLikelihood::
tieParameters(const std::vector<size_t> & par_indices) {
   TiedParameter * tiedPar = new TiedParameter();
   for (std::vector<size_t>::const_iterator it(par_indices.begin());
        it != par_indices.end(); ++it) {
      tiedPar->addParam(*m_masterComponent, *it);
      if (m_tiedIndices.find(*it) != m_tiedIndices.end()) {
         throw std::runtime_error("A parameter can belong to one "
                                  "group of tied parameters at most.");
      }
      m_tiedIndices.insert(*it);
   }
   m_tiedPars.push_back(tiedPar);
}

void SummedLikelihood::
setErrors(const std::vector<double> & errors) {
   if (m_components.empty()) {
      throw std::runtime_error("SummedLikelihood::setErrors: "
                               "no log-likelihood components");
   }
   size_t j(0);
   std::vector<optimizers::Parameter> & pars(m_masterComponent->parameters());
   for (size_t par_index(0); par_index < pars.size(); par_index++) {
      if (pars[par_index].isFree() &&
          m_tiedIndices.find(par_index) == m_tiedIndices.end()) {
         pars[par_index].setError(errors.at(j++));
      }
   }
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars[i]->isFree()) {
         m_tiedPars[i]->setError(errors.at(j++));
      }
   }
   for (ComponentIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      (*it)->setParams(pars);
   }
   syncParams();
}

int SummedLikelihood::
findIndex(size_t par_index) const {
   int j(0);
   const std::vector<optimizers::Parameter> & 
      pars(m_masterComponent->parameters());
   for (size_t i(0); i < pars.size(); i++) {
      if (!pars[i].isFree() || m_tiedIndices.find(i) != m_tiedIndices.end()) {
         continue;
      }
      if (par_index == i) {
         return j;
      }
      j++;
   }
   for (size_t i(0); i < m_tiedPars.size(); i++) {
      if (m_tiedPars[i]->containsIndex(*m_masterComponent, par_index)) {
         return j;
      }
      j++;
   }
   return -1;
}

TiedParameter & SummedLikelihood::
getTiedParam(size_t par_index) {
   throw std::runtime_error("not implemented");
   return *m_tiedPars.front();
}

void SummedLikelihood::
setTiedParamValue(size_t i, double value) {
   throw std::runtime_error("not implemented");
}

} // namespace Likleihood
