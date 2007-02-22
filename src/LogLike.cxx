/** 
 * @file LogLike.cxx
 * @brief LogLike class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LogLike.cxx,v 1.57 2007/02/19 18:05:32 jchiang Exp $
 */

#include <cmath>
#include <ctime>

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "st_stream/StreamFormatter.h"

#include "Likelihood/Accumulator.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Npred.h"
#include "Likelihood/SrcArg.h"

namespace Likelihood {

LogLike::LogLike(const Observation & observation) 
   : SourceModel(observation), m_nevals(0) {
   deleteAllSources();
}

double LogLike::value(optimizers::Arg&) const {
   std::clock_t start = std::clock();
   const std::vector<Event> & events = m_observation.eventCont().events();
   double my_value(0);
   
// The "data sum"
   for (size_t j = 0; j < events.size(); j++) {
      my_value += logSourceModel(events.at(j));
      m_accumulator.add(logSourceModel(events.at(j)));
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
         my_value -= m_Npred(sArg);
         m_accumulator.add(-m_Npred(sArg));
      }
   }
   double my_total(m_accumulator.total());
   st_stream::StreamFormatter formatter("LogLike", "value", 4);
   formatter.info() << m_nevals << "  "
                    << my_value << "  "
                    << my_total << "  "
                    << std::clock() - start << std::endl;
   m_nevals++;
//   return my_value;
   return my_total;
}

double LogLike::logSourceModel(const Event & event) const {
   double my_value(0);
   if (m_useNewImp) {
      for (size_t i = 0; i < m_freeSrcs.size(); i++) {
         const_cast<Event &>(event).updateModelSum(*m_freeSrcs.at(i));
      }
      my_value = event.modelSum();
   } else {
      std::map<std::string, Source *>::const_iterator 
         source(m_sources.begin());
      for ( ; source != m_sources.end(); ++source) {
         double fluxDens(source->second->fluxDensity(event));
         my_value += fluxDens;
      }
   }
   if (my_value > 0) {
      return std::log(my_value);
   }
   return 0;
}

void LogLike::getLogSourceModelDerivs(const Event & event,
                                      std::vector<double> & derivs) const {
   derivs.clear();
   derivs.reserve(getNumFreeParams());
   double srcSum = std::exp(logSourceModel(event));

   std::map<std::string, Source *>::const_iterator source = m_sources.begin();
   for ( ; source != m_sources.end(); ++source) {
      Source::FuncMap srcFuncs = source->second->getSrcFuncs();
      Source::FuncMap::const_iterator func_it = srcFuncs.begin();
      for (; func_it != srcFuncs.end(); func_it++) {
         std::vector<std::string> paramNames;
         (*func_it).second->getFreeParamNames(paramNames);
         for (size_t j = 0; j < paramNames.size(); j++) {
            derivs.push_back(
               source->second->fluxDensityDeriv(event, paramNames[j])/srcSum
               );
         }
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
      getLogSourceModelDerivs(events[j], derivs);
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
      freeDerivs.push_back(logSrcModelDerivs[i] - NpredDerivs[i]);
   }
}

void LogLike::addSource(Source * src) {
   SourceModel::addSource(src);
   const std::vector<Event> & events = m_observation.eventCont().events();
   for (size_t j = 0; j < events.size(); j++) {
      const_cast<std::vector<Event> &>(events).at(j).updateModelSum(*src);
   }
   SrcArg sArg(src);
   m_npredValues[src->getName()] = m_Npred(sArg);
}

Source * LogLike::deleteSource(const std::string & srcName) {
   const std::vector<Event> & events = m_observation.eventCont().events();
   for (size_t j = 0; j < events.size(); j++) {
      const_cast<std::vector<Event> &>(events).at(j).deleteSource(srcName);
   }
   m_npredValues.erase(srcName);
   return SourceModel::deleteSource(srcName);
}

void LogLike::getEvents(std::string event_file) {
   EventContainer & eventCont =
      const_cast<EventContainer &>(m_observation.eventCont());
   eventCont.getEvents(event_file);
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
         const_cast<Event &>(events.at(j)).updateModelSum(*source->second);
      }
   }
}

} // namespace Likelihood
