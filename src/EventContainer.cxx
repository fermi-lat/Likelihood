/**
 * @file EventContainer.cxx
 * @brief Container for FT1 event data.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/EventContainer.cxx,v 1.25 2012/06/14 04:18:39 jchiang Exp $
 */

#include <cmath>

#include <algorithm>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "irfInterface/IEfficiencyFactor.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/EventContainer.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"

namespace Likelihood {

std::vector<std::string> EventContainer::s_FT1_columns;

EventContainer::EventContainer(const ResponseFunctions & respFuncs, 
                               const RoiCuts & roiCuts, const ScData & scData) 
   : m_respFuncs(respFuncs), m_roiCuts(roiCuts), m_scData(scData),
     m_formatter(new st_stream::StreamFormatter("EventContainer", "", 2)) {
   if (s_FT1_columns.size() == 0) {
      setFT1_columns();
   }
}

EventContainer::~EventContainer() {
   delete m_formatter;
}

void EventContainer::getEvents(std::string event_file) {

   facilities::Util::expandEnvVar(&event_file);

   unsigned int nTotal(0);
   unsigned int nReject(0);

   facilities::Util::expandEnvVar(&event_file);
   tip::Table * events = 
      tip::IFileSvc::instance().editTable(event_file, "events");

   int evclsver(0); // version of event class definition

   tip::Header & header(events->getHeader());
   try {
      header["EVCLSVER"].get(evclsver);
   } catch(tip::TipException) {
      // keyword missing so use default value
   }

   double ra;
   double dec;
   double energy;
   double time;
   double zenAngle;
   int eventClass;
   int conversionType;
   int eventType;

   double respValue;

   tip::Table::Iterator it = events->begin();
   tip::Table::Record & event = *it;

   std::vector<std::string> diffuseNames;
   DiffRespNames diffRespNames;
   bool haveOldDiffRespCols(false);
   try {
      int ndifrsp;
      header["NDIFRSP"].get(ndifrsp);
      get_diffuse_names(events, diffRespNames);
      diffuseNames = diffRespNames.colnames();
   } catch(tip::TipException) {
// Use old diffuse response column names.
      get_diffuse_names(events, diffuseNames);
      haveOldDiffRespCols = true;
   }

   for ( ; it != events->end(); ++it, nTotal++) {
      event["ra"].get(ra);
      event["dec"].get(dec);
      event["energy"].get(energy);
      event["time"].get(time);
      event["zenith_angle"].get(zenAngle);
      event["conversion_type"].get(conversionType);
      event["event_class"].get(eventClass);
      if (evclsver == 0) {
         eventType = conversionType;
      } else {
         eventType = conversionType + 2*eventClass;
      }
      const irfInterface::IEfficiencyFactor * eff_factor =
         m_respFuncs.respPtr(eventType)->efficiencyFactor();

      double efficiency(1);
      if (eff_factor) {
         efficiency = eff_factor->value(energy, m_scData.livetimefrac(time),
                                        time);
      }
      Event thisEvent(ra, dec, energy, time, m_scData.zAxis(time),
                      m_scData.xAxis(time), cos(zenAngle*M_PI/180.), 
                      m_respFuncs.useEdisp(), m_respFuncs.respName(),
                      eventType, efficiency);
      if (m_roiCuts.accept(thisEvent)) {
         m_events.push_back(thisEvent);
         for (std::vector<std::string>::iterator name = diffuseNames.begin();
              name != diffuseNames.end(); ++name) {
            std::string srcName = sourceName(*name);
            std::vector<double> gaussianParams;
            if (m_respFuncs.useEdisp()) {
               try {
                  event[*name].get(gaussianParams);
                  m_events.back().setDiffuseResponse(srcName, gaussianParams);
               } catch (tip::TipException & eObj) {
                  std::string message(eObj.what());
                  if (message.find("FitsColumn::getVector") ==
                      std::string::npos) {
                     throw;
                  }
               }
            } else {
               std::string colname;
               if (haveOldDiffRespCols) {
                  colname = *name;
               } else {
                  colname = diffRespNames.key(*name);
               }
               event[colname].get(respValue);
               m_events.back().setDiffuseResponse(*name, respValue);
            }
         }
      } else {
         nReject++;
      }
   }

   m_formatter->info(3) << "EventContainer::getEvents:\nOut of " 
                        << nTotal << " events in file "
                        << event_file << ",\n "
                        << nTotal - nReject << " were accepted, and "
                        << nReject << " were rejected.\n" << std::endl;

   delete events;
}

void EventContainer::computeEventResponses(Source & src, double sr_radius) {
                      
   DiffuseSource *diffuse_src = dynamic_cast<DiffuseSource *>(&src);
   std::vector<DiffuseSource *> srcs;
   srcs.push_back(diffuse_src);
   computeEventResponses(srcs, sr_radius);
//    m_formatter->warn() << "Computing Event responses for " << src.getName();
//    for (unsigned int i = 0; i < m_events.size(); i++) {
//       if ((i % (m_events.size()/20)) == 0) {
//          m_formatter->warn() << ".";
//       }
//       m_events[i].computeResponse(*diffuse_src, m_respFuncs, sr_radius);
//    }
//    m_formatter->warn() << "!" << std::endl;
}

void EventContainer::computeEventResponses(std::vector<DiffuseSource *> &srcs,
                                           double sr_radius) {
   (void)(sr_radius);
   if (m_events.size() == 0) {
      return;
   }
   std::vector<DiffuseSource *> new_srcs;
   m_events[0].getNewDiffuseSrcs(srcs, new_srcs);
   if (new_srcs.size() > 0) {
      m_formatter->info(2) << "Computing Event responses for "
                           << "the DiffuseSources:\n";
      for (size_t j = 0; j < new_srcs.size(); j++) {
         m_formatter->info(2) << new_srcs.at(j)->getName() << "\n";
      }
   }
   for (unsigned int i = 0; i < m_events.size(); i++) {
      if (m_events.size() > 20 && (i % (m_events.size()/20)) == 0) {
         m_formatter->info(3) << ".";
      }
// Use Gaussian quadrature calculation instead of default. 
/// @todo Implement an accurate and faster default calculation.
//       m_events[i].computeResponse(srcs, m_respFuncs, sr_radius);
      m_events[i].computeResponseGQ(srcs, m_respFuncs);
   }
   m_formatter->info(3) << "!" << std::endl;
}

std::vector<double> 
EventContainer::nobs(const std::vector<double> & ebounds,
                     const Source * src) const {
   std::vector<double> my_nobs(ebounds.size()-1, 0);
   for (size_t i(0); i < m_events.size(); i++) {
      const Event & event(m_events.at(i));
      double energy(event.getEnergy());
      if (energy < ebounds.front() || energy > ebounds.back()) {
         continue;
      }
      size_t k = std::upper_bound(ebounds.begin(), ebounds.end(), energy)
         - ebounds.begin() - 1;
      if (energy == ebounds.front()) {
         k = 0;
      }
      if (energy == ebounds.back()) {
         k = my_nobs.size() - 1;
      }
      if (src) {
         my_nobs.at(k) += (src->fluxDensity(event)*event.efficiency()
                           /event.modelSum());
      } else {
         my_nobs.at(k)++;
      }
   }
   return my_nobs;
}

void EventContainer::setFT1_columns() const {
   std::string colnames("energy ra dec l b theta phi zenith_angle "
                        + std::string("earth_azimuth_angle time event_id ")
                        + "recon_version calib_version "
                        + std::string("event_class conversion_type ")
                        + "livetime pulse_phase mc_src_id orbital_phase");
   facilities::Util::stringTokenize(colnames, " ", s_FT1_columns);
}

void EventContainer::
get_diffuse_names(tip::Table * events, 
                  std::vector<std::string> & names) const {
/// Read in the diffuse respnose columns assuming the old column naming
/// scheme.
   names.clear();
   const std::vector<std::string> & fields(events->getValidFields());
   std::string respName(m_respFuncs.respName());
   Event::toLower(respName);
   for (size_t i = 0; i < fields.size(); i++) {
      if (fields.at(i).find("__") != std::string::npos ||
          fields.at(i).find("::") != std::string::npos) {
         if (fields.at(i).find(respName) == 0) {
            names.push_back(fields.at(i));
         }
      }
   }
}
   
void EventContainer::
get_diffuse_names(tip::Table * events, 
                  DiffRespNames & diffRespNames) const {
   const tip::Header & header(events->getHeader());
   int nkeys;
   header["NDIFRSP"].get(nkeys);
   for (int i(0); i < nkeys; i++) {
      std::ostringstream keyname;
      keyname << "DIFRSP" << i;
      std::string colname;
      header[keyname.str()].get(colname);
      diffRespNames.addColumn(colname);
   }
}
   
std::string EventContainer::sourceName(const std::string & name) const {
// The column name for a diffuse response has the IRF name prepended.
// Strip the IRF name and use the underlying diffuse component name
// in setDiffuseResponse.
   std::vector<std::string> tokens;
   std::string::size_type pos;
   if ((pos = name.find("__")) != std::string::npos) {
      return name.substr(pos+2);
   } else if ((pos = name.find("::")) != std::string::npos) {
      return name.substr(pos+2);
   }
   return name;
}

} // namespace Likelihood
