/**
 * @file LikeExposure.h
 * @brief Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LikeExposure.h,v 1.7 2005/03/07 05:18:29 jchiang Exp $
 */

#ifndef Likelihood_LikeExposure_h
#define Likelihood_LikeExposure_h

#include <utility>

#include "map_tools/Exposure.h"

namespace tip {
   class Table;
}

namespace Likelihood {

/**
 * @class LikeExposure
 *
 * @brief Class to aid in computing an exposure time hypercube that
 * includes the ROI time range cuts and GTIs.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LikeExposure.h,v 1.7 2005/03/07 05:18:29 jchiang Exp $
 */

class LikeExposure : public map_tools::Exposure {

public:

   LikeExposure(double skybin, double costhetabin,
                const std::vector< std::pair<double, double> > & timeCuts,
                const std::vector< std::pair<double, double> > & gtis);

   void load(tip::Table * tuple, bool verbose=true);

   /// @param start MET start time of interval (seconds)
   /// @param stop MET stop time of interval (seconds)
   /// @param timeCuts Time range cuts
   /// @param gtis Good Time Intervals
   /// @param fraction Fraction of the interval to use in exposure
   ///        calculation.  This is a return value.
   static bool 
   acceptInterval(double start, double stop, 
                  const std::vector< std::pair<double, double> > & timeCuts,
                  const std::vector< std::pair<double, double> > & gtis,
                  double & fraction);

   static double overlap(const std::pair<double, double> & interval1,
                         const std::pair<double, double> & interval2);

private:

   const std::vector< std::pair<double, double> > & m_timeCuts;
   const std::vector< std::pair<double, double> > & m_gtis;

   static bool overlaps(const std::pair<double, double> & interval1,
                        std::pair<double, double> & interval2);

};

} // namespace Likelihood

#endif // Likelihood_LikeExposure_h
