/** @file logLike_ptsrc.cxx
 * @brief logLike_ptsrc class implementation
 *
 * $Header:
 */

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "logLike_ptsrc.h"

namespace Likelihood {

/* compute the EML log-likelihood for a single-point source */

double logLike_ptsrc::value(const std::vector<double> &paramVec) {
   setParamValues(paramVec);
   
   double my_value = 0;
   
// the "data sum"
   for (unsigned int j = 0; j < m_events.size(); j++) {
      double src_sum = 0.;
      src_sum += evaluate_at(m_events[j]);
      if (src_sum > 0) my_value += log(src_sum);
   }

// the "model integral", a sum over Npred for each source
   for (unsigned int i = 0; i < getNumSrcs(); i++) {
      my_value -= m_sources[i]->Npred();
   }
   
   return my_value;
}

double logLike_ptsrc::evaluate_at(const Event &evt) const {
   double my_value = 0;
   for (unsigned int i = 0; i < getNumSrcs(); i++) {
      my_value += (*m_sources[i]).fluxDensity(evt);
   }
   return my_value;
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
