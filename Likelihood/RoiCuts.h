/** @file RoiCuts.h
 * @brief Declaration for RoiCuts class
 * @author J. Chiang
 * 
 * $Header$
 */

#ifndef RoiCuts_h
#define RoiCuts_h

#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include "astro/SkyDir.h"
#include "Likelihood/Event.h"
#include "Likelihood/ScData.h"

namespace Likelihood {

/** 
 * @class RoiCuts
 *
 * @brief NTuple Singleton class to represent Region-of-interest cuts.
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class RoiCuts {

public:

   ~RoiCuts(){}

   static RoiCuts * instance();

   //! access to the cuts
   void getTimeCuts(std::vector< std::pair<double, double> > &tLimVec) const
      {tLimVec = s_tLimVec;}

   std::pair<double, double> getEnergyCuts() const
      {return std::make_pair(s_eMin, s_eMax);}

   std::pair<astro::SkyDir, double> getExtractionRegion() const
      {return std::make_pair(s_roiCenter, s_roiRadius);}

   double getMuZenMax() {return s_muZenMax;}

   //! methods to allow cuts to be specified
   // prompt user (or use PIL for defaults)
   static void setCuts(double ra = 193.98, double dec = -5.82, 
                       double roi_radius = 50);
//   static void setCuts(std::string RoiFile);    // read from (XML?) file

   //! apply these cuts to an Event
   bool accept(const Event &);

protected:

   RoiCuts(){}

private:

   static RoiCuts * s_instance;

   //! cuts on photon "MET" arrival times in seconds; 
   //! this vector of pairs specify time intervals for event acceptance;
   //! the *intersection* of these intervals will be used
   typedef std::pair<double, double> timeInterval; // this will be generalized
   static std::vector<timeInterval> s_tLimVec;

   //! minimum and maximum energies in MeV,
   static double s_eMin;
   static double s_eMax;

   //! center of the sky extraction region
   static astro::SkyDir s_roiCenter;

   //! acceptance cone half angle in degrees
   static double s_roiRadius;

   //! cosine of the maximum Zenith angle
   static double s_muZenMax;

};

} // namespace Likelihood
#endif // RoiCuts_h
