/**
 * @file CompositeLikelihood.cxx
 * @brief Statistic object that comprises LogLike objects that have a
 * source (e.g. DM source) with common fit parameters intended to be
 * tied together.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#include <sstream>
#include <stdexcept>

#include "Likelihood/CompositeLikelihood.h"

namespace Likelihood {

CompositeLikelihood::CompositeLikelihood() : optimizers::Statistic() {}

CompositeLikelihood::~CompositeLikelihood() throw() {}

void CompositeLikelihood::addComponent(const std::string & srcName, 
                                       LogLike & component) {
   if (!component.hasSrcNamed(srcName)) {
      std::ostringstream message;
      message << "Log-likelihood component does not have source named "
              << srcName;
      throw std::runtime_error(message.str());
   }
   m_components[srcName] = &component;
}

double CompositeLikelihood::value() const {
   double my_value(0);
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
      my_value += it->second->value();
   }
   return my_value;
}

void CompositeLikelihood::
getFreeParams(std::vector<optimizers::Parameter> & params) const {
   if (m_components.empty()) {
      throw std::runtime_error("getFreeParams: empty composite list");
   }

// Loop over LogLike components and gather up the free Parameters for
// all of the Sources, except for the common source types, which we
// reserve to the end.
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
      const std::string & commonSrcName(it->first);
      std::map<std::string, Source *>::const_iterator 
         src(it->second->sources().begin());
      for ( ; src != it->second->sources().end(); ++src) {
         if (src->first != commonSrcName) { 
            std::vector<optimizers::Parameter> my_params;
            src->second->spectrum().getFreeParams(my_params);
            for (size_t i(0); i < my_params.size(); i++) {
               params.push_back(my_params.at(i));
            }
         }
      }
   }

// Loop over LogLike components again, this time gathering up the
// common source type params (for the first component) and
// normalization parameters for subsequent ones.
   std::vector<optimizers::Parameter> my_params;
   for (it = m_components.begin(); it != m_components.end(); ++it) {
      const std::string & commonSrcName(it->first);
      const Source * my_source = 
         it->second->sources().find(commonSrcName)->second;
      if (it == m_components.begin()) {
         my_source->spectrum().getFreeParams(my_params);
         for (size_t i(0); i < my_params.size(); i++) {
            params.push_back(my_params.at(i));
         }
      } else {
         const optimizers::Parameter & normPar = 
            const_cast<optimizers::Function &>(my_source->spectrum()).normPar();
         if (normPar.isFree()) {
            params.push_back(normPar);
         }
      }
   }
}

void CompositeLikelihood::
setFreeParamValues(const std::vector<double> & values) {
   if (m_components.empty()) {
      throw std::runtime_error("getFreeParams: empty composite list");
   }

// Loop over LogLike components and set the free Parameters for all of
// the Sources, except for the common source types, which we handle
// at the end.
   std::vector<double>::const_iterator vals(values.begin());
   ComponentIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
      const std::string & commonSrcName(it->first);
      std::map<std::string, Source *>::const_iterator 
         src(it->second->sources().begin());
      for ( ; src != it->second->sources().end(); ++src) {
         if (src->first != commonSrcName) { 
            size_t npar(src->second->spectrum().getNumFreeParams());
            std::vector<double> my_values(npar);
            for (size_t i(0); i < npar; i++, ++vals) {
               my_values.push_back(*vals);
            }
            const_cast<optimizers::Function &>(src->second->spectrum())
               .setFreeParamValues(my_values);
         }
      }
   }

// Set the common source parameters for the first component, saving the
// values for the non-normPar parameters for the later components.
   it = m_components.begin();
   const std::string & commonSrcName(it->first);
   Source * my_source = const_cast<Source *>(it->second->sources()
                                             .find(commonSrcName)->second);
   optimizers::Function & my_spectrum = 
      const_cast<optimizers::Function &>(my_source->spectrum());

   std::vector<std::string> parNames;
   my_spectrum.getFreeParamNames(parNames);
   std::string normParName(my_spectrum.normPar().getName());

   std::map<std::string, double> parValues;
   for (size_t i(0); i < parNames.size(); i++, ++vals) {
      if (parNames.at(i) != normParName) {
         parValues[parNames.at(i)] = *vals;
      }
      my_spectrum.setParam(parNames.at(i), *vals);
   }
   ++it;

// Loop over remaining LogLike components, setting the common source
// type params for all components and the normalization parameters
// where needed.
   for ( ; it != m_components.end(); ++it) {
      const std::string & commonSrcName(it->first);
      Source * my_source = const_cast<Source *>(it->second->sources()
                                                .find(commonSrcName)->second);
      optimizers::Function & my_spectrum = 
         const_cast<optimizers::Function &>(my_source->spectrum());
      for (size_t i(0); i < parNames.size(); i++) {
         my_spectrum.setParam(parNames.at(i), parValues[parNames.at(i)]);
      }
      if (my_spectrum.normPar().isFree()) {
         my_spectrum.setParam(normParName, *vals);
         ++vals;
      }
   }
}

unsigned int CompositeLikelihood::getNumFreeParams() const {
   ComponentConstIterator_t it(m_components.begin());
   unsigned int npars(it->second->getNumFreeParams());
   return npars + m_components.size() - 1;
}

// void CompositeLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
//    std::vector<double> par_derivs;
//    std::vector<double> norm_derivs;
//    ComponentConstIterator_t it(m_components.begin());
//    for ( ; it != m_components.end(); ++it) {


//    size_t nsrcs(m_components.size());
// //
// // Get the derivatives wrt to the free parameters of the first component.
// // component values, then reset the normalization parameters.
// // 
//    ComponentConstIterator_t it(m_components.begin());
//    for ( ; it != m_components.end(); ++it) {
//       it->second->getFreeDerivs(derivs);
//    }
//    for (size_t i(nsrcs-1), it=m_components.end()-1; i > 0; --it, i--) {
//       it->second->normPar().setValue(values.at(i));
//    }
   
// }

} // namespace Likleihood
