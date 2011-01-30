/** 
 * @file LogLike.cxx
 * @brief LogLike class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/LogLike.cxx,v 1.76 2011/01/29 06:53:03 jchiang Exp $
 */

#include <cmath>
#include <ctime>

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "st_stream/StreamFormatter.h"

#include "Likelihood/Exception.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Npred.h"
#include "Likelihood/SrcArg.h"

namespace Likelihood {

LogLike::LogLike(const Observation & observation) 
  : SourceModel(observation), m_nevals(0), m_bestValueSoFar(-1e38),
    m_Npred(), m_accumulator(), m_npredValues(),    
    m_respCache() {
   const std::vector<Event> & events = m_observation.eventCont().events();
   m_respCache.clearAndResize(events.size());
   deleteAllSources();
}

double LogLike::value(optimizers::Arg&) const {
   std::clock_t start = std::clock();
   const std::vector<Event> & events = m_observation.eventCont().events();
   double my_value(0);
   
// The "data sum"
   for (size_t j = 0; j < events.size(); j++) {
      ResponseCache::EventRef rc_ref = m_respCache.getEventRef(j);
      double addend(logSourceModel(events.at(j), &rc_ref));
      my_value += addend;
      m_accumulator.add(addend);
   }

// The "model integral", a sum over Npred for each source
   if (m_useNewImp) {
      std::map<std::string, double>::const_iterator 
         npred(m_npredValues.begin());
      for ( ; npred != m_npredValues.end(); ++npred) {
         my_value -= npred->second;
         m_accumulator.add(-npred->second);
      }
   } else {
      std::map<std::string, Source *>::const_iterator srcIt(m_sources.begin());
      for ( ; srcIt != m_sources.end(); ++srcIt) {
         SrcArg sArg(srcIt->second);
	 double addend = m_Npred(sArg);
         my_value -= addend;
         m_accumulator.add(-addend);
      }
   }
   double my_total(m_accumulator.total());
   st_stream::StreamFormatter formatter("LogLike", "value", 4);
   formatter.info() << m_nevals << "  "
                    << my_value << "  "
                    << my_total << "  "
                    << std::clock() - start << std::endl;
   m_nevals++;
//    if (::getenv("BYPASS_ACCUMULATOR")) {
//       return my_value;
//    }

/// Add in contribution from priors.
   std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
   for ( ; par != m_parameter.end(); ++par) {
      if (par->isFree()) {
         my_total += par->log_prior_value();
      }
   }

   saveBestFit(my_total);
   return my_total;
}

double LogLike::logSourceModel(const Event & event,
			       ResponseCache::EventRef* srcRespCache) const {
   double my_value(0);
// This part was commented out in v15r8p2 (Feb 8, 2010), either for
// accuracy reasons or because there was some problem related to the
// handling of the free state of sources.
   // if (m_useNewImp) {
   //    for (size_t i = 0; i < m_freeSrcs.size(); i++) {
   //           const Source* source = m_freeSrcs.at(i);
   //           CachedResponse* cResp=0;
   //           if(srcRespCache)
   //            cResp = &srcRespCache->getCachedValue(source->getName());
   //        const_cast<Event &>(event).updateModelSum(*m_freeSrcs.at(i), cResp);
   //    }
   //    my_value = event.modelSum();
   // } else {
      std::map<std::string, Source *>::const_iterator 
         source(m_sources.begin());
      for ( ; source != m_sources.end(); ++source) {
         CachedResponse* cResp=0;
         if (srcRespCache) {
            cResp = &srcRespCache->getCachedValue(source->second->getName());
         }
// Event::modelSum() will be used for the per event source
// probabilities so we need to update the Event::m_modelSum value.
         if (std::count(m_freeSrcs.begin(), m_freeSrcs.end(), source->second)) {
            const_cast<Event &>(event).updateModelSum(*source->second, cResp);
         }
         double fluxDens(source->second->fluxDensity(event, cResp));
         fluxDens *= event.efficiency();
         my_value += fluxDens;
      }
//       if (my_value > 0 && 
//           std::fabs((my_value - event.modelSum())/my_value) > 1e-5) {
//          std::cout << event.getEnergy() << "  "
//                    << my_value << "  "
//                    << event.modelSum() << std::endl;
//       }
   // }
   if (my_value > 0) {
      return std::log(my_value);
   }
//    throw std::runtime_error("negative probability density for this event.");
   return 0;
}

void LogLike::getLogSourceModelDerivs(const Event & event,
                                      std::vector<double> & derivs,
				 ResponseCache::EventRef* srcRespCache) const {
   derivs.clear();
   derivs.reserve(getNumFreeParams());
   double srcSum = std::exp(logSourceModel(event,srcRespCache));

   std::map<std::string, Source *>::const_iterator source = m_sources.begin();
   for ( ; source != m_sources.end(); ++source) {
//       Source::FuncMap srcFuncs = source->second->getSrcFuncs();
//       Source::FuncMap::const_iterator func_it = srcFuncs.begin();
//       CachedResponse* cResp=0;
//       for (; func_it != srcFuncs.end(); func_it++) {
//          std::vector<std::string> paramNames;
//          (*func_it).second->getFreeParamNames(paramNames);
//          // Only set cResp for sources with at least 1 free param
//          if((cResp==0)&&(!paramNames.empty())&&(srcRespCache))
//             cResp = &srcRespCache->getCachedValue(source->second->getName());
//          for (size_t j = 0; j < paramNames.size(); j++) {
//             double fluxDensDeriv = 
// 	       source->second->fluxDensityDeriv(event, paramNames[j], cResp);
//             fluxDensDeriv *= event.efficiency();
//             derivs.push_back(fluxDensDeriv/srcSum);
//          }
//       }
      CachedResponse * cResp(0);
      std::vector<std::string> paramNames;
      source->second->spectrum().getFreeParamNames(paramNames);
      if ( (cResp == 0) && (!paramNames.empty()) && (srcRespCache) ) {
         cResp = &srcRespCache->getCachedValue(source->second->getName());
      }
      for (size_t j(0); j < paramNames.size(); j++) {
         double fluxDensDeriv = 
            source->second->fluxDensityDeriv(event, paramNames.at(j), cResp);
         fluxDensDeriv *= event.efficiency();
         derivs.push_back(fluxDensDeriv/srcSum);
      }
   }
}

void LogLike::getFreeDerivs(optimizers::Arg&,
                            std::vector<double> &freeDerivs) const {
// Retrieve the free derivatives for the log(SourceModel) part
   const std::vector<Event> & events = m_observation.eventCont().events();

   std::vector<double> logSrcModelDerivs(getNumFreeParams(), 0);
   for (size_t j = 0; j < events.size(); j++) {
      std::vector<double> derivs;
      ResponseCache::EventRef rc_ref = m_respCache.getEventRef(j);
      getLogSourceModelDerivs(events[j], derivs, &rc_ref);
      for (size_t i = 0; i < derivs.size(); i++) {
         logSrcModelDerivs[i] += derivs[i];
      }
   }

// The free derivatives for the Npred part must be appended 
// for each Source in m_sources.
   std::vector<double> NpredDerivs;
   NpredDerivs.reserve(getNumFreeParams());

   if (m_useNewImp) {
      std::vector<Source *>::const_iterator srcIt = m_freeSrcs.begin();
      for ( ; srcIt != m_freeSrcs.end(); ++srcIt) {
         SrcArg sArg(*srcIt);
         std::vector<double> derivs;
         m_Npred.getFreeDerivs(sArg, derivs);
         for (size_t i = 0; i < derivs.size(); i++) {
            NpredDerivs.push_back(derivs[i]);
         }
      }
   } else {
      std::map<std::string, Source *>::const_iterator srcIt 
         = m_sources.begin();
      for ( ; srcIt != m_sources.end(); ++srcIt) {
         SrcArg sArg(srcIt->second);
         std::vector<double> derivs;
         m_Npred.getFreeDerivs(sArg, derivs);
         for (size_t i = 0; i < derivs.size(); i++) {
            NpredDerivs.push_back(derivs[i]);
         }
      }
   }

   freeDerivs.reserve(NpredDerivs.size());
   freeDerivs.clear();
   for (size_t i = 0; i < NpredDerivs.size(); i++) {
      freeDerivs.push_back(logSrcModelDerivs.at(i) - NpredDerivs.at(i));
   }

   /// Derivatives from priors.
   size_t i(0);
   std::vector<optimizers::Parameter>::const_iterator par(m_parameter.begin());
   for ( ; par != m_parameter.end(); ++par) {
      if (par->isFree()) {
         freeDerivs[i] += par->log_prior_deriv();
         i++;
      }
   }
}

void LogLike::addSource(Source * src) {
   SourceModel::addSource(src);
   const std::vector<Event> & events = m_observation.eventCont().events();
   std::string srcName = src->getName();

   bool useCachedResp(false);
   if (m_useNewImp) {
      // New implementation: use the response cache if there are free parameters
//       Source::FuncMap & srcFuncs(src->getSrcFuncs());
//       for (Source::FuncMap::const_iterator func = srcFuncs.begin();
//            func != srcFuncs.end(); ++func) {
//          if (func->second->getNumFreeParams() > 0) {
//             useCachedResp = true;
//             break;
//          }
//       }
      if (src->spectrum().getNumFreeParams() > 0) {
         useCachedResp = true;
      }
   } else {
      // Old implementation: always use the response cache
      useCachedResp = true;
   }

   for (size_t j = 0; j < events.size(); j++) {
      CachedResponse* cResp = 0;
      if(useCachedResp)cResp = &m_respCache.getCachedValue(j,srcName);
      const_cast<std::vector<Event> &>(events).at(j).updateModelSum(*src, cResp);
   }
   SrcArg sArg(src);
   m_npredValues[src->getName()] = m_Npred(sArg);
   m_bestValueSoFar = -1e38;
}

Source * LogLike::deleteSource(const std::string & srcName) {
   const std::vector<Event> & events = m_observation.eventCont().events();
   for (size_t j = 0; j < events.size(); j++) {
      const_cast<std::vector<Event> &>(events).at(j).deleteSource(srcName);
   }
   m_respCache.deleteSource(srcName);
   m_npredValues.erase(srcName);
   m_bestValueSoFar = -1e38;
   return SourceModel::deleteSource(srcName);
}

void LogLike::getEvents(std::string event_file) {
   EventContainer & eventCont =
      const_cast<EventContainer &>(m_observation.eventCont());
   eventCont.getEvents(event_file);
   m_respCache.clearAndResize(eventCont.events().size());
}

void LogLike::computeEventResponses(double sr_radius) {
   std::vector<DiffuseSource *> diffuse_srcs;
   std::map<std::string, Source *>::iterator srcIt = m_sources.begin();
   for ( ; srcIt != m_sources.end(); ++srcIt) {
      if (srcIt->second->getType() == std::string("Diffuse")) {
         DiffuseSource *diffuse_src = 
            dynamic_cast<DiffuseSource *>(srcIt->second);
         diffuse_srcs.push_back(diffuse_src);
      }
   }
   if (diffuse_srcs.size() > 0) {
      EventContainer & eventCont =
         const_cast<EventContainer &>(m_observation.eventCont());
      eventCont.computeEventResponses(diffuse_srcs, sr_radius);
   }
}

void LogLike::syncParams() {
   SourceModel::syncParams();
   if (m_useNewImp) {
      for (size_t i = 0; i < m_freeSrcs.size(); i++) {
         SrcArg sArg(m_freeSrcs.at(i));
         m_npredValues[m_freeSrcs.at(i)->getName()] = m_Npred(sArg);
      }
   }
}

void LogLike::syncSrcParams(const std::string & srcName) {
   SourceModel::syncParams();
   std::map<std::string, Source *>::const_iterator source 
      = m_sources.find(srcName);
   if (source != m_sources.end()) {
      SrcArg sArg(source->second);
      m_npredValues[source->first] = m_Npred(sArg);
      const std::vector<Event> & events(m_observation.eventCont().events());
      for (size_t j(0); j < events.size(); j++) {
	 CachedResponse* cResp = &m_respCache.getCachedValue(j,srcName);
         const_cast<Event &>(events.at(j)).updateModelSum(*source->second,cResp);
      }
   }
}

double LogLike::NpredValue(const std::string & srcName) const {
   return const_cast<Source &>(source(srcName)).Npred();
}

void LogLike::saveBestFit(double logLikeValue) const {
   // This is called from value(...), so we must pass the current
   // log-likelihood value as a parameter, i.e., we cannot evaluate it
   // in this function.
   if (logLikeValue > m_bestValueSoFar) {
      getParamValues(m_bestFitParsSoFar);
      m_bestValueSoFar = logLikeValue;
   }
}

void LogLike::restoreBestFit() {
   setParamValues(m_bestFitParsSoFar);
   syncParams();
}

void LogLike::saveCurrentFit() {
   syncParams();
   m_bestValueSoFar = -1e38;
   optimizers::Arg dummy;
   saveBestFit(value(dummy));
}

void LogLike::addPrior(size_t index,
                       optimizers::Function & log_prior) {
   m_parameter[index].setPrior(log_prior);
}

} // namespace Likelihood
