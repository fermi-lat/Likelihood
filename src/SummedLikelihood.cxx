/**
 * @file SummedLikelihood.cxx
 * @brief Statistic object that comprises LogLike objects that have
 * identical source models.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SummedLikelihood.cxx,v 1.8 2008/11/26 23:35:02 jchiang Exp $
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Likelihood/SummedLikelihood.h"

namespace Likelihood {

void SummedLikelihood::addComponent(LogLike & component) {
   m_components.push_back(&component);
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
      throw std::runtime_error("getFreeParams: empty log-likelihood "
                               "component list");
   }
// All components have the same free parameters, so just need to
// get the ones from the first component.
   m_components.front()->getFreeParams(params);
}

void SummedLikelihood::
setFreeParamValues(const std::vector<double> & values) {
   if (m_components.empty()) {
      throw std::runtime_error("getFreeParams: empty composite list");
   }
   for (ComponentIterator_t it(m_components.begin());
        it != m_components.end(); ++it) {
      (*it)->setFreeParamValues(values);
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
   return m_components.front()->getNumFreeParams();
}

void SummedLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
// Build vector of component parameter names.  This list should have
// the same ordering as the derivatives wrt the free parameters.
   m_components.front()->getFreeDerivs(derivs);
   for (ComponentIterator_t it(m_components.begin() + 1);
        it != m_components.end(); ++it) {
      std::vector<double> freeDerivs;
      (*it)->getFreeDerivs(freeDerivs);
      for (size_t i(0); i < derivs.size(); i++) {
         derivs.at(i) += freeDerivs.at(i);
      }
   }
}

} // namespace Likleihood
