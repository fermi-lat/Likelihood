/** 
 * @file LogLike.cxx
 * @brief LogLike class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LogLike.cxx,v 1.14 2004/04/03 06:14:14 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <stdexcept>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

// This include is temporary until tip provides column access.
#include "latResponse/../src/DC1.h"

#include "Likelihood/ScData.h"
#include "Likelihood/Npred.h"
#include "Likelihood/logSrcModel.h"
#include "Likelihood/EventArg.h"
#include "Likelihood/SrcArg.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/LogLike.h"

namespace Likelihood {

double LogLike::value(optimizers::Arg&) const {
// Compute the EML log-likelihood for a single-point source.

   double my_value = 0;
   
// The "data sum"
   for (unsigned int j = 0; j < m_events.size(); j++) {
      EventArg eArg(m_events[j]);
      my_value += m_logSrcModel(eArg);
   }

// The "model integral", a sum over Npred for each source
   for (unsigned int i = 0; i < getNumSrcs(); i++) {
      SrcArg sArg(s_sources[i]);
      my_value -= m_Npred(sArg);
   }
   return my_value;
}

void LogLike::getFreeDerivs(optimizers::Arg&,
                            std::vector<double> &freeDerivs) const {

// retrieve the free derivatives for the log(SourceModel) part
   m_logSrcModel.mySyncParams();
   std::vector<double> logSrcModelDerivs(m_logSrcModel.getNumFreeParams(), 0);
   for (unsigned int j = 0; j < m_events.size(); j++) {
      std::vector<double> derivs;
      EventArg eArg(m_events[j]);
      m_logSrcModel.getFreeDerivs(eArg, derivs);
      for (unsigned int i = 0; i < derivs.size(); i++)
         logSrcModelDerivs[i] += derivs[i];
   }

// the free derivatives for the Npred part must be appended 
// for each Source in s_sources
   std::vector<double> NpredDerivs;
   NpredDerivs.reserve(m_logSrcModel.getNumFreeParams());

   for (unsigned int i = 0; i < s_sources.size(); i++) {
      SrcArg sArg(s_sources[i]);
      std::vector<double> derivs;
      m_Npred.getFreeDerivs(sArg, derivs);
      for (unsigned int i = 0; i < derivs.size(); i++) 
         NpredDerivs.push_back(derivs[i]);
   }

   freeDerivs.reserve(NpredDerivs.size());
   freeDerivs.clear();
   for (unsigned int i = 0; i < NpredDerivs.size(); i++) 
      freeDerivs.push_back(logSrcModelDerivs[i] - NpredDerivs[i]);
}

void LogLike::computeEventResponses(Source &src, double sr_radius) {
   DiffuseSource *diffuse_src = dynamic_cast<DiffuseSource *>(&src);
   std::cerr << "Computing Event responses for " << src.getName();
   for (unsigned int i = 0; i < m_events.size(); i++) {
      if ((i % (m_events.size()/20)) == 0) std::cerr << ".";
      m_events[i].computeResponse(*diffuse_src, sr_radius);
   }
   std::cerr << "!" << std::endl;
}

void LogLike::computeEventResponses(std::vector<DiffuseSource *> &srcs, 
                                    double sr_radius) {
   std::cerr << "Computing Event responses for the DiffuseSources";
   for (unsigned int i = 0; i < m_events.size(); i++) {
      if ((i % (m_events.size()/20)) == 0) std::cerr << ".";
      m_events[i].computeResponse(srcs, sr_radius);
   }
   std::cerr << "!" << std::endl;
}

void LogLike::computeEventResponses(double sr_radius) {
   std::vector<DiffuseSource *> diffuse_srcs;
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (s_sources[i]->getType() == std::string("Diffuse")) {
         DiffuseSource *diffuse_src = 
            dynamic_cast<DiffuseSource *>(s_sources[i]);
         diffuse_srcs.push_back(diffuse_src);
      }
   }
   if (diffuse_srcs.size() > 0) computeEventResponses(diffuse_srcs, sr_radius);
}

#ifdef USE_FT1
void LogLike::getEvents(std::string event_file, int) {

   facilities::Util::expandEnvVar(&event_file);

   RoiCuts * roiCuts = RoiCuts::instance();
   ScData * scData = ScData::instance();

   if (!scData) {
      std::cout << "LogLike::getEvents: "
                << "The spacecraft data must be read in first."
                << std::endl;
      assert(scData);
   }

   if (!roiCuts) {
      std::cout << "LogLike::getEvents: "
                << "The region-of-interest data must be read in first."
                << std::endl;
      assert(roiCuts);
   }

   unsigned int nTotal(0);
   unsigned int nReject(0);

   facilities::Util::expandEnvVar(&event_file);
   tip::Table * events = 
      tip::IFileSvc::instance().editTable(event_file, "events");

   double ra;
   double dec;
   double energy;
   double time;
   double raSCZ;
   double decSCZ;
   double zenAngle;
   int convLayer;
   int eventType;

   tip::Table::Iterator it = events->begin();
   tip::Table::Record & event = *it;
   for ( ; it != events->end(); ++it, nTotal++) {
      event["ra"].get(ra);
      event["dec"].get(dec);
      event["energy"].get(energy);
      event["time"].get(time);
      raSCZ = scData->zAxis(time).ra();
      decSCZ = scData->zAxis(time).dec();
      event["zenith_angle"].get(zenAngle);
      event["conversion_layer"].get(convLayer);
      if (convLayer < 12) { // Front
         eventType = 0;
      } else {
         eventType = 1;
      }
      Event thisEvent(ra, dec, energy, time, raSCZ, decSCZ, 
                      cos(zenAngle*M_PI/180.), eventType);
      if (roiCuts->accept(thisEvent)) {
         m_events.push_back(thisEvent);
      } else {
            nReject++;
      }
   }

   std::cout << "LogLike::getEvents:\nOut of " 
             << nTotal << " events in file "
             << event_file << ",\n "
             << nTotal - nReject << " were accepted, and "
             << nReject << " were rejected.\n" << std::endl;

   delete events;
}
#else // USE_FT1
void LogLike::getEvents(std::string event_file, int hdu) {

   facilities::Util::expandEnvVar(&event_file);

   readEventData(event_file, hdu);

   typedef std::pair<long, std::vector<double> > tableColumn;
   tableColumn ra = getEventColumn("RA");
   tableColumn dec = getEventColumn("DEC");
   tableColumn energy = getEventColumn("energy");
   tableColumn time = getEventColumn("time");
   tableColumn sc_x = getEventColumn("SC_x");
   tableColumn sc_y = getEventColumn("SC_y");
   tableColumn sc_z = getEventColumn("SC_z");
   tableColumn zenangle = getEventColumn("zenith_angle");

// get pointer to RoiCuts
   RoiCuts *roi_cuts = RoiCuts::instance();

   unsigned int nReject = 0;

   long nevents = m_events.size();
   m_events.reserve(ra.first + nevents);
   for (int i = 0; i < ra.first; i++) {
// compute sc_ra and sc_dec from direction cosines 
      double sc_ra = atan2(sc_y.second[i], sc_x.second[i])*180./M_PI;
      double sc_dec = asin(sc_z.second[i])*180./M_PI;

      Event thisEvent(ra.second[i], dec.second[i], energy.second[i],
                      time.second[i], sc_ra, sc_dec,
                      cos(zenangle.second[i]*M_PI/180.));

      if (roi_cuts->accept(thisEvent)) {
         m_events.push_back(thisEvent);
      } else {
         nReject++;
      }
   }

   std::cout << "LogLike::getEvents:\nOut of " 
             << ra.first << " events in file "
             << event_file << ",\n "
             << m_events.size() - nevents << " were accepted, and "
             << nReject << " were rejected.\n" << std::endl;
}

void LogLike::readEventData(const std::string &eventFile, int hdu) {
   m_eventFile = eventFile;
   m_eventHdu = hdu;

   std::string extName;
   latResponse::DC1::getFitsHduName(eventFile, hdu, extName);

   std::vector<std::string> columnNames;
   latResponse::DC1::getFitsColNames(eventFile, hdu, columnNames);

   std::vector<double> my_column;
   for (unsigned int i = 0; i < columnNames.size(); i++) {
      latResponse::DC1::getTableVector(eventFile, extName, columnNames[i],
                                       my_column);
      m_eventColumns[columnNames[i]] = my_column;
   }
}

std::pair<long, std::vector<double> >
LogLike::getEventColumn(const std::string &colname) const {
   std::pair<long, std::vector<double> > my_column;

   if (m_eventColumns.count(colname)) {
      my_column = std::make_pair(m_eventColumns[colname].size(), 
                                 m_eventColumns[colname]);
   } else {
      std::ostringstream errorMessage;
      errorMessage << "LogLike::getColumn:\n"
                   << "Column " << colname 
                   << " was not found in event data.\n"
                   << "Valid names are \n" << colnames << "\n";
      throw std::runtime_error(errorMessage.str());
   }
   return my_column;
}
#endif // USE_FT1

} // namespace Likelihood
