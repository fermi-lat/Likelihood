/**
 * @file LikeExposure.cxx
 * @brief Implementation of Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header$
 */

#include "facilities/Util.h"

#include "astro/EarthCoordinate.h"

#include "tuple/ITable.h"

#include "Likelihood/RoiCuts.h"
#include "Likelihood/LikeExposure.h"

namespace Likelihood {

LikeExposure::LikeExposure(double skybin, double costhetabin, 
                           const std::string & roiFile) 
   : Exposure(skybin, costhetabin) {

   RoiCuts::setCuts(roiFile);
   RoiCuts::instance()->getTimeCuts(m_timeCuts);
}

void LikeExposure::load(tuple::ITable & scData) { 
   
   const double & ra = scData.selectColumn("ra_scz");
   const double & dec = scData.selectColumn("dec_scz");
   const double & start = scData.selectColumn("start");
   const double & stop = scData.selectColumn("stop");
   const double & livetime = scData.selectColumn("livetime");
   const double & latGeo = scData.selectColumn("lat_geo");
   const double & lonGeo = scData.selectColumn("lon_geo");

   for (tuple::Iterator it = scData.begin(); it != scData.end(); ++it) {
      double deltat = livetime > 0 ? livetime : stop-start;
      double fraction;
      if (acceptInterval(start, stop, latGeo, lonGeo, fraction)) {
         deltat *= fraction;
         add(astro::SkyDir(ra, dec), deltat);
      }
   }
}

bool LikeExposure::acceptInterval(double start, double stop, 
                                  double latGeo, double lonGeo, 
                                  double & fraction) {
                                  
   std::vector< std::pair<double, double> >::const_iterator it 
      = m_timeCuts.begin();

// At some point, we need to implement fraction properly.
   fraction = 1.;

   bool accept(true);
   for ( ; it != m_timeCuts.end(); it++) {
      if (stop < it->first || start > it->second) {
         accept = false;
         break;
      }
   }

   (void)(latGeo); (void)(lonGeo);
//    astro::EarthCoordinate earthCoord(latGeo, lonGeo);
//    if (!earthCoord.insideSAA) {
//       accept = false;
//    }
      
   return accept;
}

} // namespace Likelihood
