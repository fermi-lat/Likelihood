/** 
 * @file LogLike.cxx
 * @brief LogLike class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LogLike.cxx,v 1.31 2004/10/11 01:34:59 jchiang Exp $
 */

#include <cmath>
#include <cassert>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/EventArg.h"
#include "Likelihood/logSrcModel.h"
#include "Likelihood/Npred.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SrcArg.h"

#include "Likelihood/LogLike.h"

#include "Verbosity.h"

namespace Likelihood {

std::vector<std::string> LogLike::s_FT1_columns;

double LogLike::value(optimizers::Arg&) const {
   double my_value = 0;
   
// The "data sum"
   for (unsigned int j = 0; j < m_events.size(); j++) {
      EventArg eArg(m_events[j]);
      my_value += m_logSrcModel(eArg);
   }

// The "model integral", a sum over Npred for each source
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      SrcArg sArg(srcIt->second);
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

   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      SrcArg sArg(srcIt->second);
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
   if (print_output()) {
      std::cerr << "Computing Event responses for " << src.getName();
   }
   for (unsigned int i = 0; i < m_events.size(); i++) {
      if (print_output() && (i % (m_events.size()/20)) == 0) std::cerr << ".";
      m_events[i].computeResponse(*diffuse_src, sr_radius);
   }
   if (print_output()) std::cerr << "!" << std::endl;
}

void LogLike::computeEventResponses(std::vector<DiffuseSource *> &srcs, 
                                    double sr_radius) {
   if (print_output()) {
      std::cerr << "Computing Event responses for the DiffuseSources";
   }
   for (unsigned int i = 0; i < m_events.size(); i++) {
      if (print_output() && m_events.size() > 20 &&
          (i % (m_events.size()/20)) == 0) std::cerr << ".";
      m_events[i].computeResponse(srcs, sr_radius);
//       if (i < 10) {
//          std::ostringstream filename;
//          filename << "diffuse_response_" << i << ".dat";
//          m_events[i].writeDiffuseResponses(filename.str());
//       }
   }
   if (print_output()) std::cerr << "!" << std::endl;
// // Write out the diffuse responses.
//    std::ofstream outfile("diffuse_responses.dat");
//    for (unsigned int i = 0; i < m_events.size(); i++) {
//       for (unsigned int j = 0; j < srcs.size(); j++) {
//          outfile << m_events[i].diffuseResponse(1., srcs[j]->getName())
//                  << "  ";
//       }
//       outfile << "\n";
//    }
//    outfile.close();
}

void LogLike::computeEventResponses(double sr_radius) {
   std::vector<DiffuseSource *> diffuse_srcs;
   std::map<std::string, Source *>::iterator srcIt = s_sources.begin();
   for ( ; srcIt != s_sources.end(); ++srcIt) {
      if (srcIt->second->getType() == std::string("Diffuse")) {
         DiffuseSource *diffuse_src = 
            dynamic_cast<DiffuseSource *>(srcIt->second);
         diffuse_srcs.push_back(diffuse_src);
      }
   }
   if (diffuse_srcs.size() > 0) {
      computeEventResponses(diffuse_srcs, sr_radius);
   }
}

void LogLike::getEvents(std::string event_file) {

   facilities::Util::expandEnvVar(&event_file);

   RoiCuts * roiCuts = RoiCuts::instance();
   ScData * scData = ScData::instance();

   if (!scData) {
      std::string message = std::string("LogLike::getEvents: ")
         + "The spacecraft data must be read in first.";
      throw std::runtime_error(message);
   }

   if (!roiCuts) {
      std::string message = std::string("LogLike::getEvents: ")
         + "The region-of-interest data must be read in first.";
      throw std::runtime_error(message);
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
   double respValue;

   tip::Table::Iterator it = events->begin();
   tip::Table::Record & event = *it;

   std::vector<std::string> diffuseNames;
   get_diffuse_names(events, diffuseNames);

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
         for (std::vector<std::string>::iterator name = diffuseNames.begin();
              name != diffuseNames.end(); ++name) {
            std::vector<double> gaussianParams;
            if (ResponseFunctions::useEdisp()) {
               try {
                  event[*name].get(gaussianParams);
                  m_events.back().setDiffuseResponse(*name, gaussianParams);
               } catch (tip::TipException & eObj) {
                  std::string message(eObj.what());
                  if (message.find_first_of("FitsColumn::getVector") ==
                      std::string::npos) {
                     throw;
                  }
               }
            } else {
               try {
                  event[*name].get(gaussianParams);
                  m_events.back().setDiffuseResponse(*name, gaussianParams[0]);
               } catch (tip::TipException &eObj) {
                  std::string message(eObj.what());
                  if (message.find_first_of("FitsColumn::getVector") !=
                      std::string::npos) {
                     event[*name].get(respValue);
                     m_events.back().setDiffuseResponse(*name, respValue);
                  } else {
                     throw;
                  }
               }
            }
         }
      } else {
         nReject++;
      }
   }

   if (print_output()) {
      std::cerr << "LogLike::getEvents:\nOut of " 
                << nTotal << " events in file "
                << event_file << ",\n "
                << nTotal - nReject << " were accepted, and "
                << nReject << " were rejected.\n" << std::endl;
   }

   delete events;
}

void LogLike::setFT1_columns() {
   std::string colnames("energy ra dec theta phi zenith_angle "
                        + std::string("earth_azimuth_angle time event_id ")
                        + "recon_version calib_version imgoodcalprob "
                        + std::string("imvertexprob imcoreprob impsferrpred ")
                        + "calenergysum caltotrln imgammaprob "
                        + std::string("conversion_point conversion_layer"));
   facilities::Util::stringTokenize(colnames, " ", s_FT1_columns);
}

void LogLike::get_diffuse_names(tip::Table * events, 
                                std::vector<std::string> & names) {
   names.clear();
   const std::vector<std::string> & fields = events->getValidFields();
   for (unsigned int i = 0; i < fields.size(); i++) {
      if (!std::count(s_FT1_columns.begin(), s_FT1_columns.end(), fields[i])) {
         names.push_back(fields[i]);
      }
   }
}

} // namespace Likelihood
