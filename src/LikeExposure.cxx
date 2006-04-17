/**
 * @file LikeExposure.cxx
 * @brief Implementation of Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LikeExposure.cxx,v 1.18 2006/04/01 00:04:03 jchiang Exp $
 */

#include <algorithm>
#include <iostream>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "tip/Table.h"

#include "Likelihood/LikeExposure.h"
#include "Likelihood/RoiCuts.h"

namespace {
   bool compareFirst(const std::pair<double, double> & a, 
                     const std::pair<double, double> & b) {
      return a.first < b.first;
   }
   bool compareSecond(const std::pair<double, double> & a, 
                      const std::pair<double, double> & b) {
      return a.second < b.second;
   }
}

namespace Likelihood {

LikeExposure::
LikeExposure(double skybin, double costhetabin, 
             const std::vector< std::pair<double, double> > & timeCuts,
             const std::vector< std::pair<double, double> > & gtis)
   : map_tools::Exposure(skybin, costhetabin), m_timeCuts(timeCuts),
     m_gtis(gtis) {}

void LikeExposure::load(const tip::Table * scData, bool verbose) {
   
   double ra, dec, start, stop, livetime;

   tip::Table::ConstIterator it = scData->begin();
   tip::ConstTableRecord & row = *it;
   long nrows = scData->getNumRecords();

   double maxTime(0);
   for (unsigned int i=0; i < m_timeCuts.size(); i++) {
      if (m_timeCuts.at(i).second > maxTime) {
         maxTime = m_timeCuts.at(i).second;
      }
   }
   for (unsigned int i=0; i < m_gtis.size(); i++) {
      if (m_gtis.at(i).second > maxTime) {
         maxTime = m_gtis.at(i).second;
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
   
   st_stream::StreamFormatter formatter("LikeExposure", "load", 2);
   for (long irow = 0; it != scData->end(); ++it, ++irow) {
      if (verbose && (irow % (nrows/20)) == 0 ) {
         formatter.info() << "."; 
      }
      row["livetime"].get(livetime);
      row["start"].get(start);
      row["stop"].get(stop);
      if (start > maxTime) {
         break;
      }
      double deltat = livetime;
      double fraction;
      if (acceptInterval(start, stop, m_timeCuts, m_gtis, fraction)) {
         row["ra_scz"].get(ra);
         row["dec_scz"].get(dec);
         fill(astro::SkyDir(ra, dec), deltat*fraction);
      }
   }
   if (verbose) {
      formatter.info() << "!" << std::endl;
   }
}

bool LikeExposure::
acceptInterval(double start, double stop, 
               const std::vector< std::pair<double, double> > & timeCuts,
               const std::vector< std::pair<double, double> > & gtis,
               double & fraction) {
                                  
   std::pair<double, double> candidateInterval(start, stop);

   typedef std::vector< std::pair<double, double> > IntervalCont_t;
   IntervalCont_t::const_iterator it;

   for (it = timeCuts.begin(); it != timeCuts.end(); ++it) {
      if (!overlaps(*it, candidateInterval)) {
         fraction = 0;
         return false;
      }
   }
   
   double total(0);
   double maxTotal(candidateInterval.second - candidateInterval.first);
   
   IntervalCont_t::const_iterator gti_first = 
      std::upper_bound(gtis.begin(), gtis.end(), candidateInterval,
                       ::compareFirst);
   if (gti_first != gtis.begin()) {
      --gti_first;
   }

   IntervalCont_t::const_iterator gti_last = 
      std::upper_bound(gti_first, gtis.end(), candidateInterval,
                       ::compareSecond);

   if (gti_last != gtis.end()) {
      ++gti_last;
   }

   for (it = gti_first; it != gti_last; ++it) {
      double dt(overlap(*it, candidateInterval));
      total += dt;
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
