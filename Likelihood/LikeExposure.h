/**
 * @file LikeExposure.h
 * @brief Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LikeExposure.h,v 1.4 2004/09/17 01:19:22 jchiang Exp $
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
 * includes the ROI time-interval cuts, SAA passages, etc..
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LikeExposure.h,v 1.4 2004/09/17 01:19:22 jchiang Exp $
 */

class LikeExposure : public map_tools::Exposure {

public:

   LikeExposure() {}

   LikeExposure(double skybin, double costhetabin);

   void load(tip::Table * tuple, bool verbose=true);

private:

   std::vector< std::pair<double, double> > m_timeCuts;

   /// @param start MET start time of interval (seconds)
   /// @param stop MET stop time of interval (seconds)
   /// @param latGeo Ground point latitude (degrees)
   /// @param lonGeo Ground point longitude (degrees)
   /// @param fraction Fraction of the interval to use in exposure
   ///        calculation.  This is a return value.
   bool acceptInterval(double start, double stop, 
                       double latGeo, double lonGeo, double & fraction);

};

} // namespace Likelihood

#endif // Likelihood_LikeExposure_h
