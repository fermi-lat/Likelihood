/**
 * @file LikeExposure.cxx
 * @brief Implementation of Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LikeExposure.cxx,v 1.1 2004/03/11 05:19:36 jchiang Exp $
 */

#include "facilities/Util.h"

#include "astro/EarthCoordinate.h"

#include "tip/Table.h"

#include "Likelihood/RoiCuts.h"
#include "Likelihood/LikeExposure.h"

namespace Likelihood {

LikeExposure::LikeExposure(double skybin, double costhetabin, 
                           const std::string & roiFile) 
   : Exposure(skybin, costhetabin) {

   RoiCuts::setCuts(roiFile);
   RoiCuts::instance()->getTimeCuts(m_timeCuts);
}

void LikeExposure::load(tip::Table * scData) {
   
   double ra, dec, start, stop, livetime, latGeo, lonGeo;

   tip::Table::Iterator it = scData->begin();
   tip::Table::Record & row = *it;
   for ( ; it != scData->end(); ++it) {
      row["livetime"].get(livetime);
      row["start"].get(start);
      row["stop"].get(stop);
      row["latGeo"].get(latGeo);
      row["lonGeo"].get(lonGeo);
      double deltat = livetime > 0 ? livetime : stop-start;
      double fraction;
      if (acceptInterval(start, stop, latGeo, lonGeo, fraction)) {
         deltat *= fraction;
         row["ra"].get(ra);
         row["dec"].get(dec);
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
