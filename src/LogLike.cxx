/** 
 * @file LogLike.cxx
 * @brief LogLike class implementation
 * @author J. Chiang
 *
 * $Header$
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "Likelihood/LogLike.h"
#include "Likelihood/Npred.h"
#include "Likelihood/logSrcModel.h"
#include "Likelihood/EventArg.h"
#include "Likelihood/SrcArg.h"
#include "Likelihood/DiffuseSource.h"

namespace Likelihood {

/* compute the EML log-likelihood for a single-point source */

double LogLike::value(optimizers::Arg&) const {

   double my_value = 0;
   
// the "data sum"
   for (unsigned int j = 0; j < m_events.size(); j++) {
      EventArg eArg(m_events[j]);
      my_value += m_logSrcModel(eArg);
   }

// the "model integral", a sum over Npred for each source
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

void LogLike::getEvents(const std::string &event_file, int hdu) {

   std::string colnames = "RA DEC energy time SC_x SC_y SC_z zenith_angle";

   readEventData(event_file, colnames, hdu);

   typedef std::pair<long, double*> tableColumn;
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
             << event_file << ", "
             << m_events.size() - nevents << " were accepted, and "
             << nReject << " were rejected.\n" << std::endl;
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

void LogLike::computeEventResponses(std::vector<DiffuseSource> &srcs, 
                                          double sr_radius) {
   std::cerr << "Computing Event responses for the DiffuseSources";
   for (unsigned int i = 0; i < m_events.size(); i++) {
      if ((i % (m_events.size()/20)) == 0) std::cerr << ".";
      m_events[i].computeResponse(srcs, sr_radius);
   }
   std::cerr << "!" << std::endl;
}

void LogLike::computeEventResponses(double sr_radius) {
   std::vector<DiffuseSource> diffuse_srcs;
   for (unsigned int i = 0; i < s_sources.size(); i++) {
      if (s_sources[i]->getType() == std::string("Diffuse")) {
         DiffuseSource *diffuse_src = 
            dynamic_cast<DiffuseSource *>(s_sources[i]);
         diffuse_srcs.push_back(*diffuse_src);
      }
   }
   computeEventResponses(diffuse_srcs, sr_radius);
}

// Methods from the old Likelihood::Statistic class:
//
// read in the event data
void LogLike::readEventData(const std::string &eventFile, 
                                  const std::string &colnames, int hdu) {
   m_eventFile = eventFile;
   m_eventHdu = hdu;
//   m_eventData.add_columns(colnames);

   delete m_eventData;
   m_eventData = new Table();

   std::vector<std::string> columnNames;
   m_eventData->read_FITS_colnames(m_eventFile, m_eventHdu, columnNames);
   m_eventData->add_columns(columnNames);
   std::cerr << "Columns in " << m_eventFile 
             << ", HDU " << m_eventHdu 
             << ": \n";
   for (unsigned int i = 0; i < columnNames.size(); i++) {
      std::cerr << columnNames[i] << "  ";
   }
   std::cerr << std::endl;
   m_eventData->read_FITS_table(m_eventFile, m_eventHdu);

}

//! return pointer to data columns
std::pair<long, double*> 
LogLike::getColumn(const Table &tableData, 
                         const std::string &colname) const
   throw(optimizers::Exception) {
   std::pair<long, double*> my_column(0, 0);
   std::string colnames;
   
// loop over column names, return the matching one
   for (int i = 0; i < tableData.npar(); i++) {
      if (colname == std::string(tableData[i].colname)) {
         my_column.first = tableData[i].dim;
         my_column.second = tableData[i].val;
         return my_column;
      }
      colnames += " "; colnames += tableData[i].colname;
   }
   std::ostringstream errorMessage;
   errorMessage << "LogLike::getColumn:\n"
                << "Column " << colname << " was not found in event data.\n"
                << "Valid names are \n" << colnames << "\n";
   throw Exception(errorMessage.str());
}

} // namespace Likelihood
