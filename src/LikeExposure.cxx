/**
 * @file LikeExposure.cxx
 * @brief Implementation of Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LikeExposure.cxx,v 1.11 2005/03/07 05:55:37 jchiang Exp $
 */

#include <algorithm>
#include <iostream>

#include "facilities/Util.h"

#include "tip/Table.h"

#include "Likelihood/LikeExposure.h"
#include "Likelihood/RoiCuts.h"

namespace Likelihood {

LikeExposure::
LikeExposure(double skybin, double costhetabin, 
             const std::vector< std::pair<double, double> > & timeCuts,
             const std::vector< std::pair<double, double> > & gtis)
   : map_tools::Exposure(skybin, costhetabin), m_timeCuts(timeCuts),
     m_gtis(gtis) {}

void LikeExposure::load(tip::Table * scData, bool verbose) {
   
   double ra, dec, start, stop, livetime;

   tip::Table::Iterator it = scData->begin();
   tip::Table::Record & row = *it;
   long nrows = scData->getNumRecords();

   double maxTime(3.2e8);  // No random access iterator in tip, so we
                           // have to assume some large value.
   for (unsigned int i=0; i < m_timeCuts.size(); i++) {
      if (m_timeCuts.at(i).second < maxTime) {
         maxTime = m_timeCuts.at(i).second;
      }
   }
// We assume that the time intervals are 30 sec long, even though the
// time_candle source in flux is hacked to behave incorrectly such
// that we cannot extract the time interval size from the data; plus
// there is again no random access iterator available to allow us to
// infer the time interval size a priori.  In any case, nrows is just
// used for calculating the progress bar.
   if (static_cast<long>(maxTime/30.) < nrows) {
      nrows = static_cast<long>(maxTime/30.);
   }
   
   for (long irow = 0; it != scData->end(); ++it, ++irow) {
      if (verbose && (irow % (nrows/20)) == 0 ) std::cerr << "."; 
      row["livetime"].get(livetime);
      row["start"].get(start);
      row["stop"].get(stop);
      if (start > maxTime) {
         break;
      }
      double deltat = livetime > 0 ? livetime : stop-start;
      double fraction;
      if (acceptInterval(start, stop, m_timeCuts, m_gtis, fraction)) {
         row["ra_scz"].get(ra);
         row["dec_scz"].get(dec);
         fill(astro::SkyDir(ra, dec), deltat*fraction);
      }
   }
   if (verbose) std::cerr << "!" << std::endl;
}

bool LikeExposure::
acceptInterval(double start, double stop, 
               const std::vector< std::pair<double, double> > & timeCuts,
               const std::vector< std::pair<double, double> > & gtis,
               double & fraction) {
                                  
   std::pair<double, double> candidateInterval(start, stop);

   std::vector< std::pair<double, double> >::const_iterator it 
      = timeCuts.begin();

   for ( ; it != timeCuts.end(); ++it) {
      if (!overlaps(*it, candidateInterval)) {
         fraction = 0;
         return false;
      }
   }
   
   double total(0);
   double maxTotal(candidateInterval.second - candidateInterval.first);
   for (it = gtis.begin(); it != gtis.end(); ++it) {
      total += overlap(*it, candidateInterval);
   }
   if (total > maxTotal || gtis.size() == 0) {
      total = maxTotal;
   }
   fraction = total/(stop - start);
   return true;
}

bool LikeExposure::
overlaps(const std::pair<double, double> & interval1,
         std::pair<double, double> & interval2) {
   double start = std::max(interval1.first, interval2.first);
   double stop = std::min(interval1.second, interval2.second);
   if (start < stop) {
      interval2.first = start;
      interval2.second = stop;
      return true;
   }
   return false;
}

double LikeExposure::
overlap(const std::pair<double, double> & interval1,
        const std::pair<double, double> & interval2) {
   double start = std::max(interval1.first, interval2.first);
   double stop = std::min(interval1.second, interval2.second);
   if (start < stop) {
      return stop - start;
   }
   return 0;
}

} // namespace Likelihood
