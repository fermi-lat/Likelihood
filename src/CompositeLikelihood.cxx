/**
 * @file CompositeLikelihood.cxx
 * @brief Statistic object that comprises LogLike objects that have a
 * source (e.g. DM source) with common fit parameters intended to be
 * tied together.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/CompositeLikelihood.cxx,v 1.9 2010/05/17 21:17:49 jchiang Exp $
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Likelihood/CompositeLikelihood.h"

namespace Likelihood {

void CompositeLikelihood::addComponent(const std::string & srcName, 
                                       LogLike & component) {
   if (!component.hasSrcNamed(srcName)) {
      std::ostringstream message;
      message << "Log-likelihood component does not have source named "
              << srcName;
      throw std::runtime_error(message.str());
   }
   Source * my_source(component.sources().find(srcName)->second);
   if (m_components.empty()) {
      m_normParName = const_cast<optimizers::Function &>(my_source->spectrum())
         .normPar().getName();
      m_commonFuncName = my_source->spectrum().genericName();
   } else {
      if (m_commonFuncName != my_source->spectrum().genericName()) {
         throw std::runtime_error("Inconsistent common source type.");
      }
   }
   m_components[&component] = srcName;
}

double CompositeLikelihood::value() const {
   double my_value(0);
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
      my_value += it->first->value();
   }
   return my_value;
}

// void CompositeLikelihood::
// getIndices(std::vector<std::vector<size_t> > & indices) const {
//    size_t ncp(m_components.size());
//    ComponentConstIterator_t it(m_components.begin());
//    for ( ; it != m_components.end(); ++it) {
//       std::map<std::string, Source *>::const_iterator src
//          = it->first->sources().begin();
//       if (src != it->second) {
         
//       }
//    }
// }

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
      const std::string & commonSrcName(it->second);
      std::map<std::string, Source *>::const_iterator 
         src(it->first->sources().begin());
      for ( ; src != it->first->sources().end(); ++src) {
         if (src->first != commonSrcName) { 
            std::vector<optimizers::Parameter> my_params;
            src->second->spectrum().getFreeParams(my_params);
            for (size_t i(0); i < my_params.size(); i++) {
               params.push_back(my_params.at(i));
            }
         }
      }
   }

// Collect the common source type params for the first component, excluding its
// normalization parameter.
   it = m_components.begin();
   const std::string & commonSrcName(it->second);
   const Source * my_source(it->first->sources().find(commonSrcName)->second);
   std::vector<optimizers::Parameter> my_params;
   my_source->spectrum().getFreeParams(my_params);
   for (size_t i(0); i < my_params.size(); i++) {
      if (my_params.at(i).getName() != m_normParName) {
         params.push_back(my_params.at(i));
      }
   }

// Loop over components again and collect the normalization parameters
// for all of the common source type objects.
   for ( ; it != m_components.end(); ++it) {
      const Source * my_source(it->first->sources().find(it->second)->second);
      const optimizers::Parameter & normPar = 
         const_cast<optimizers::Function &>(my_source->spectrum()).normPar();
      if (normPar.isFree()) {
         params.push_back(normPar);
      }
   }
}

void CompositeLikelihood::
fetchParamValues(std::vector<double> &values, bool getFree) const {
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
      std::map<std::string, Source *>::const_iterator 
         src(it->first->sources().begin());
      for ( ; src != it->first->sources().end(); ++src) {
         if (src->first != it->second) { 
            size_t npar(src->second->spectrum().getNumFreeParams());
            std::vector<double> my_values;
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
   Source * my_source = const_cast<Source *>(it->first->sources()
                                             .find(it->second)->second);
   optimizers::Function & my_spectrum = 
      const_cast<optimizers::Function &>(my_source->spectrum());

   std::vector<std::string> parNames;
   my_spectrum.getFreeParamNames(parNames);

   std::map<std::string, double> parValues;
   for (size_t i(0); i < parNames.size(); i++) {
      if (parNames.at(i) != m_normParName) {
         parValues[parNames.at(i)] = *vals;
         my_spectrum.setParam(parNames.at(i), *vals);
         ++vals;
      }
   }


// Loop over remaining LogLike components, setting the common source
// type params for all components and the normalization parameters
// where needed.
   for ( ; it != m_components.end(); ++it) {
      my_source = const_cast<Source *>(it->first->sources()
                                       .find(it->second)->second);
      optimizers::Function & my_spec = 
         const_cast<optimizers::Function &>(my_source->spectrum());
      for (size_t i(0); i < parNames.size(); i++) {
         if (parNames.at(i) != m_normParName) {
            my_spec.setParam(parNames.at(i), parValues[parNames.at(i)]);
         }
      }
      if (my_spec.normPar().isFree()) {
         my_spec.setParam(m_normParName, *vals);
         ++vals;
      }
   }

   syncParams();
}

void CompositeLikelihood::syncParams() {
   m_parameter.clear();
   optimizers::Parameter saved_par;
   for (ComponentIterator_t it(m_components.begin()); 
        it != m_components.end(); ++it) {
      it->first->syncParams();
      const std::vector<optimizers::Parameter> & pars(it->first->parameters());
      for (size_t i(0); i < pars.size(); i++) {
	if(pars.at(i).getName()!=m_normParName){
	  m_parameter.push_back(pars.at(i));
	}
	else {
	  saved_par=pars.at(i);
	}
      }
   }
   //Finally sync the tied parameter, taken from the last component as it does not matter
   m_parameter.push_back(saved_par);
}

unsigned int CompositeLikelihood::getNumFreeParams() const {
   ComponentConstIterator_t it(m_components.begin());
   unsigned int npars(0);
   for ( ; it != m_components.end(); ++it) {
      std::map<std::string, Source *>::const_iterator src =
         it->first->sources().begin();
      for ( ; src != it->first->sources().end(); ++src) {
         if (it->second != src->first) {
            npars += src->second->spectrum().getNumFreeParams();
         }
      }
      const optimizers::Function & my_spectrum
         = it->first->sources().find(it->second)->second->spectrum();
      if (it == m_components.begin()) {
         std::vector<optimizers::Parameter> pars;
         my_spectrum.getFreeParams(pars);
         for (size_t i(0); i < pars.size(); i++) {
            if (pars.at(i).getName() != m_normParName) {
               npars++;
            }
         }
      }
      if (const_cast<optimizers::Function &>(my_spectrum).normPar().isFree()) {
         npars++;
      }
   }
   return npars;
}

void CompositeLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
// Build vector of component parameter names.  This list should have
// the same ordering as the derivatives wrt the free parameters.
   std::vector<double> freeDerivs;
   std::vector<std::pair<std::string, std::string> > componentPars;
   std::vector<int> commonSrcTypeFlag;

   ComponentConstIterator_t it(m_components.begin());
// Determine the number of spectral shape parameters (i.e., excluding
// the normalization) for the common source type object in the first
// component.
   const optimizers::Function & my_spectrum = 
      it->first->sources().find(it->second)->second->spectrum();
   std::vector<optimizers::Parameter> pars;
   my_spectrum.getFreeParams(pars);
   size_t npars(0);
   for (size_t i(0); i < pars.size(); i++) {
      if (pars.at(i).getName() != m_normParName) {
         npars++;
      }
   }

   for ( ; it != m_components.end(); ++it) {
      std::map<std::string, Source *>::const_iterator src 
         = it->first->sources().begin();
      for ( ; src != it->first->sources().end(); ++src) {
         std::vector<std::string> parNames;
         src->second->spectrum().getFreeParamNames(parNames);
         for (size_t i(0); i < parNames.size(); i++) {
            componentPars.push_back(std::make_pair(src->first, 
                                                   parNames.at(i)));
            if (it->second == src->first) {
               commonSrcTypeFlag.push_back(1);
            } else {
               commonSrcTypeFlag.push_back(0);
            }
         }
      }
      std::vector<double> my_derivs;
      it->first->getFreeDerivs(my_derivs);
      for (size_t i(0); i < my_derivs.size(); i++) {
         freeDerivs.push_back(my_derivs.at(i));
      }
   }

// Loop through elements and append derivatives not associated with
// the common source type. For the common source type objects, collect
// normalization derivatives, and sum up spectral derivatives.
   derivs.clear();
   std::vector<double> specParDerivs(npars, 0);
   std::vector<double> normDerivs;
   size_t ipar(0);
   for (size_t i(0); i < freeDerivs.size(); i++) {
      if (commonSrcTypeFlag.at(i) == 0) {
         derivs.push_back(freeDerivs.at(i));
      } else {
         if (componentPars.at(i).second == m_normParName) {
            normDerivs.push_back(freeDerivs.at(i));
         } else {
            specParDerivs.at(ipar % npars) += freeDerivs.at(i);
            ipar++;
         }
      }
   }

// Finally, append spectral shape derivatives and then normalization 
// derivatives.
   for (size_t i(0); i < specParDerivs.size(); i++) {
      derivs.push_back(specParDerivs.at(i));
   }
   for (size_t i(0); i < normDerivs.size(); i++) {
      derivs.push_back(normDerivs.at(i));
   }
}

} // namespace Likleihood
