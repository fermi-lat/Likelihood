/**
 * @file LikeExposure.cxx
 * @brief Implementation of Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LikeExposure.cxx,v 1.23 2007/12/03 20:33:48 jchiang Exp $
 */

#include <algorithm>
#include <iostream>
#include <stdexcept>

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
     m_gtis(gtis), m_numIntervals(0) {
   if (!gtis.empty()) {
      for (size_t i = 0; i < gtis.size(); i++) {
         if (i == 0 || gtis.at(i).first < m_tmin) {
            m_tmin = gtis.at(i).first;
         }
         if (i == 0 || gtis.at(i).second > m_tmax) {
            m_tmax = gtis.at(i).second;
         }
      }
   } else {
      throw std::runtime_error("LikeExposure::LikeExposure: GTIs are empty.\n"
                               "Cannot proceed with livetime calculation.");
   }
}

void LikeExposure::load(const tip::Table * scData, bool verbose) {
   st_stream::StreamFormatter formatter("LikeExposure", "load", 2);
   
   double ra, dec, start, stop, livetime;

   tip::Table::ConstIterator it = scData->end();
   tip::ConstTableRecord & row = *it;

   --it;
   long nrows(scData->getNumRecords());
   if (nrows == 0) {
      return;
   }
   for ( ; it != scData->begin(); --it, nrows--) {
      row["stop"].get(stop);
      if (stop < m_tmax) {
         break;
      }
   }

   double last_start;
   it = scData->begin();
   for ( ; it != scData->end(); ++it, nrows--) {
      last_start = start;
      row["start"].get(start);
      if (start > m_tmin) {
         break;
      }
   }
// Reset to the FT2 interval start time that precedes the
// user-selected interval.
   start = last_start; 

   long istep(nrows/20);
   if (istep == 0) {
      istep = 1;
   }

   --it;

   for (long irow = 0; it != scData->end() && start < m_tmax; ++it, ++irow) {
      if (verbose && (irow % istep) == 0 ) {
         formatter.warn() << "."; 
      }
      row["livetime"].get(livetime);
      row["start"].get(start);
      row["stop"].get(stop);
      double deltat = livetime;
      double fraction;
      if (acceptInterval(start, stop, m_timeCuts, m_gtis, fraction)) {
         row["ra_scz"].get(ra);
         row["dec_scz"].get(dec);
         fill(astro::SkyDir(ra, dec), deltat*fraction);
         m_numIntervals++;
      }
   }
   if (verbose) {
      formatter.warn() << "!" << std::endl;
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
