/**
 * @file LikelihoodBase.cxx
 * @brief Base class for Likelihood classes.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Attic/LikelihoodBase.cxx,v 1.1.2.1 2013/01/30 16:09:30 jchiang Exp $
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
   if(m_sources.count(src->getName())) {
      throw Exception("LikelihoodBase: Source named "
		      + src->getName() + " already exists.");
   }
  
   Source * my_src = src->clone();
   my_src->setObservation(&m_observation);

   if(my_src->fixedSpectrum()) {
      addFixedSource(my_src, src); 
   } else {
      addFreeSource(my_sec, src);
   }
}

Source * LikelihoodBase::getSource(const std::string & srcName) {
   std::map<std::string, Source *>::iterator it(m_sources.find(srcName));
   if (it != m_sources.end()) {
      return it->second;
   }
   throw Exception("LikelihoodBase: Source named "
                   + src->getName() + "  not found.");
}

Source * LikelihoodBase::removeSource(const std::string & srcName) {
   Source * my_source(getSource(srcName));
   m_sources.erase(srcName);

   if(m_fixed_source_param_values.count(my_source)) {
      removeFixedSource(my_source);
   } else {
      removeFreeSource(my_source);
   }

   syncParams();
   return my_source;
}

void LikelihoodBase::addFreeSource(Source* src, const Source* callers_src) {
   (void)(callers_src);
   std::vector<Source *>::iterator ifree = m_free_sources.begin();
   for(std::map<std::string, Source *>::iterator isrc = m_sources.begin();
       isrc != m_sources.end(); isrc++) {
      if(isrc->second == src) {
         break;
      } else if(isrc->second == *ifree) {
         ifree++;
      }
   }
   m_free_sources.insert(ifree, src);
   // In a concrete likelihood imlementation should go on and calculate
   // the response value for each bin (or event) and cache them
}

void removeFreeSource(Source* src) {
   for(std::vector<Source *>::iterator ifree = m_free_sources.begin();
       ifree != ifree.end(); ifree++) {
      if(*ifree == src) {
         m_free_sources.erase(ifree);
	 return;
      }
   }
   assert(0);
   // A concrete likelihood imlementation should remove cached response values
}

void LikelihoodBase::addFixedSource(Source* src, const Source* callers_src) {
   (void)(callers_src);
   src->spectrum()->getParamValues(m_fixed_source_param_values[src]);
   // In a concrete likelihood imlementation should go on and calculate
   // the contribution of this source to each bin (or event) and to npred.
}

void LikelihoodBase::removeFixedSource(Source* src) {
   assert(m_fixed_source_param_values.erase(src));
   rebuildFixedSourceCache();
}

void LikelihoodBase::rebuildFixedSourceCache() {
   // Clear the values cached for the fixed sources and then add them back.
   // Fixed sources that have become free or whose parameter values have
   // been changed are removed from the fixed-source list and added to the
   // free-source list (from which they will never return unless the source
   // is deleted and re-added)

   std::vector<Source*> fixed_sources;
   for(std::map<Source *, std::vecor<double> >::iterator isrc = 
	 m_fixed_source_param_values.begin(); 
       isrc!=m_fixed_source_param_values.end(); isrc++) {
      std::vector<double> params;
      isrc->first->getParamValues(params);
      if((isrc->fixedSpectrum()) && (params == isrc->second)) {
	 fixed_sources.push_back(isrc->first);
      } else {
	 addFreeSource(isrc->first, 0);
      }
   }

   m_fixed_source_param_values.clear();
   for(std::vector<KahanAccumulator>::iterator ibin = 
	 m_fixed_source_flux_in_bin.begin(); 
       ibin != m_fixed_source_flux_in_bin.end(); ibin++) {
      ibin->reset();
   }
   m_fixed_source_npred.reset();

   for(std::vector<Source*>::iterator isrc = fixed_sources.begin();
       isrc != fixed_sources.end(); isrc++) {
      addFixedSource(my_src, src); 
   }
}

bool LikelihoodBase::hasFixedSourceModelChanged() const {
   // We test ONLY that sources that are on the fixed list are still
   // fixed and that they still have the same parameter values. We
   // deliberately do not test for the inverse, i.e. that there are
   // now ADDITIONAL fixed sources. The only way to get on the
   // fixed-source list is through addSource.
   for(std::map<Source *, std::vecor<double> >::const_iterator isrc = 
	 m_fixed_source_param_values.begin(); 
       isrc!=m_fixed_source_param_values.end(); isrc++) {
      std::vector<double> params;
      isrc->first->getParamValues(params);
      if((!isrc->fixedSpectrum()) || (params != isrc->second)) {
	 return true;
      }
   }
   return false;
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
