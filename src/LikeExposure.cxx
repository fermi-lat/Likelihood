/**
 * @file LikeExposure.cxx
 * @brief Implementation of Exposure class for use by the Likelihood tool.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/LikeExposure.cxx,v 1.4 2004/04/12 23:22:49 jchiang Exp $
 */

#include <iostream>

#include "facilities/Util.h"

#include "astro/EarthCoordinate.h"

#include "tip/Table.h"

#include "Likelihood/RoiCuts.h"
#include "Likelihood/LikeExposure.h"

namespace Likelihood {

LikeExposure::LikeExposure(double skybin, double costhetabin, 
                           const std::string & roiFile) 
   : map_tools::Exposure(skybin, costhetabin) {

   RoiCuts::setCuts(roiFile);
   RoiCuts::instance()->getTimeCuts(m_timeCuts);
}

void LikeExposure::load(tip::Table * scData) {
   
   double ra, dec, start, stop, livetime, latGeo, lonGeo;

   tip::Table::Iterator it = scData->begin();
   tip::Table::Record & row = *it;
   long nrows = scData->getNumRecords();
   for (long irow = 0; it != scData->end(); ++it, ++irow) {
      if ( (irow % (nrows/20)) == 0 ) std::cerr << "."; 
      row["livetime"].get(livetime);
      row["start"].get(start);
      row["stop"].get(stop);
      row["lat_Geo"].get(latGeo);
      row["lon_Geo"].get(lonGeo);
      double deltat = livetime > 0 ? livetime : stop-start;
      double fraction;
      if (acceptInterval(start, stop, latGeo, lonGeo, fraction)) {
         deltat *= fraction;
         row["ra_scz"].get(ra);
         row["dec_scz"].get(dec);
         add(astro::SkyDir(ra, dec), deltat);
      }
   }
   std::cerr << "!" << std::endl;
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
