/** @file logLike_ptsrc.cxx
 * @brief logLike_ptsrc class implementation
 * @author J. Chiang
 *
 * $Header$
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "logLike_ptsrc.h"
#include "Likelihood/Npred.h"
#include "Likelihood/logSrcModel.h"
#include "Likelihood/EventArg.h"
#include "Likelihood/SrcArg.h"

namespace Likelihood {

/* compute the EML log-likelihood for a single-point source */

double logLike_ptsrc::value(const std::vector<double> &paramVec) {
   setFreeParamValues(paramVec);
   
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

void logLike_ptsrc::getFreeDerivs(std::vector<double> &freeDerivs) {

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
      freeDerivs[i] = logSrcModelDerivs[i] - NpredDerivs[i];
}

void logLike_ptsrc::getEvents(const std::string &event_file, int hdu) {

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

   m_events.reserve(ra.first);
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

   std::cout << "logLike_ptsrc::getEvents:\nOut of " 
             << ra.first << " events, "
             << m_events.size() << " were accepted, and "
             << nReject << " were rejected.\n" << std::endl;
}

} // namespace Likelihood
