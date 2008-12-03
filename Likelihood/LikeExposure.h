/**
 * @file LikeExposure.h
 * @brief Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/LikeExposure.h,v 1.14 2008/07/05 21:44:27 jchiang Exp $
 */

#ifndef Likelihood_LikeExposure_h
#define Likelihood_LikeExposure_h

#include <utility>

#include "tip/tip_types.h"
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
 */

class LikeExposure : public map_tools::Exposure {

public:

   LikeExposure(double skybin, double costhetabin,
                const std::vector< std::pair<double, double> > & timeCuts,
                const std::vector< std::pair<double, double> > & gtis,
                double zenmax=180.);

   void load(const tip::Table * tuple, bool verbose=true);

   tip::Index_t numIntervals() const {
//   size_t numIntervals() const {
      return m_numIntervals;
   }

   /// @brief Normally one would re-implement the
   /// map_tools::Exposure::write(...) member function from the base
   /// class, but it is not virtual, so we add this method instead to
   /// avoid possible confusion if these classes are used
   /// polymorphically.
   void writeFile(const std::string & outfile) const;

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

   double m_costhetabin;

   const std::vector< std::pair<double, double> > & m_timeCuts;
   const std::vector< std::pair<double, double> > & m_gtis;

   /// Minimum time to be considered given GTIs (MET s)
   double m_tmin;

   /// Maximum time to be considered given GTIs (MET s)
   double m_tmax;

   /// Number of FT2 intervals that have been loaded.
//   size_t m_numIntervals;
   tip::Index_t m_numIntervals;

   static bool overlaps(const std::pair<double, double> & interval1,
                        std::pair<double, double> & interval2);

   void writeFilename(const std::string & outfile) const;

   void writeLivetimes(const std::string & outfile) const;

   void writeCosbins(const std::string & outfile) const;

   void setCosbinsFieldFormat(const std::string & outfile) const;

   void fitsReportError(int status, const std::string & routine) const;

   void computeCosbins(std::vector<double> & mubounds) const;

};

} // namespace Likelihood

#endif // Likelihood_LikeExposure_h
