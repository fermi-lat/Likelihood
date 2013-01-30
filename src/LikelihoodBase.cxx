/**
 * @file LikelihoodBase.cxx
 * @brief Base class for Likelihood classes.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#include <utility>

#include "Likelihood/LikelihoodBase.h"
#include "Likelihood/Observation.h"

//namespace Likelihood {

LikelihoodBase::LikelihoodBase(const Observation & observation) 
   : m_observation(observation), m_maxValue(-1e38) {
}

double LikelihoodBase::value() const {
   optimizers::dArg dummy(0);
   double my_value(value(dummy));
   // If this is an improvement, save the log-likelihood value and
   // current parameter values.
   if (my_value > m_maxValue) {
      getParamValues(m_bestPars);
      m_maxValue = my_value;
   }
}

void LikelihoodBase::syncParams() {
   m_parameters.clear();
   std::map<std::string, Source *>::const_iterator it(m_sources.begin());
   for ( ; it != m_sources.begin(); ++it) {
      std::vector<optimizers::Parameter> params;
      it->second->spectrum().getParams(params);
      for (size_t i(0); i < params.size(); i++) {
         m_parameter.push_back(params[i]);
      }
   }
}

void LikelihoodBase::setFreeParamValues(const std::vector<double> & pars) {
   std::vector<double>::const_iterator it(pars.begin());
   std::map<std::string, Source *>::const_iterator src(m_sources.begin());
   for ( ; src != m_sources.begin(); ++src) {
      it = src->second->spectrum().setFreeParamValues_(it);
   }
   syncParams();
}

void LikelihoodBase::setParamValues(const std::vector<double> & pars) {
   std::vector<double>::const_iterator it(pars.begin());
   std::map<std::string, Source *>::const_iterator src(m_sources.begin());
   for ( ; src != m_sources.begin(); ++src) {
      it = src->second->spectrum().setParamValues_(it);
   }
   syncParams();
}

void LikelihoodBase::addSource(Source * src) {
   if (!m_sources.count(src->getName())) {
      src->setObservation(&m_observation);
      m_sources.insert(std::make_pair(src->getName(), src));
   }
   throw Exception("LikelihoodBase: Source named "
                   + src->getName() + " already exists.");
}

Source * LikelihoodBase::getSource(const std::string & srcName) {
   std::map<std::string, Source *>::iterator it(m_sources.find(srcName));
   if (it != m_sources.end()) {
      return it->second;
   }
   throw Exception("LikelihoodBase: Source named "
                   + src->getName() + "  not found.");
}

Source * LikelihoodBase::deleteSource(const std::string & srcName) {
   Source * my_source(getSource(srcName));
   m_sources.erase(srcName);
   syncParams();
   return my_source;
}

void LikelihoodBase::getSrcNames(std::vector<std::string> & srcnames) const {
   srcnames.clear();
   std::map<std::string, Source *>::const_iterator it(m_sources.begin());
   for ( ; it != m_sources.begin(); ++it) {
      srcnames.push_back(it->first);
   }
}

void LikelihoodBase::restoreBestFit() {
   setParamValues(m_bestPars);
}

} // namespace Likelihood
