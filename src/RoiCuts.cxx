/** @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Regoin-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header$
 */

#include "Likelihood/RoiCuts.h"

namespace Likelihood {

// definitions of static data
std::vector<RoiCuts::timeInterval> RoiCuts::s_tLimVec;
double RoiCuts::s_eMin;
double RoiCuts::s_eMax;
astro::SkyDir RoiCuts::s_roiCenter;
double RoiCuts::s_roiRadius;
double RoiCuts::s_muZenMax;
RoiCuts * RoiCuts::s_instance = 0;

void RoiCuts::setCuts(double ra, double dec, double roi_radius) {
// get everything for now....
   s_tLimVec.push_back(std::make_pair(0., HUGE));

// default min and max energies in MeV
   s_eMin = 31.623;
   s_eMax = 3.1622e5;

// // prompt the user for sky extraction region info
//    double ra;
//    double dec;
//    std::cout << "Enter ra and dec (in degrees) of ROI center: ";
//    std::cin >> ra >> dec;
    s_roiCenter = astro::SkyDir(ra, dec);

//    std::cout << "Enter acceptance cone radius (in degrees): ";
//    std::cin >> s_roiRadius;
    s_roiRadius = roi_radius;

// accept everything until effect of Zenith angle cuts can be computed
   s_muZenMax = -1.;
}

bool RoiCuts::accept(const Event &event) {
   bool acceptEvent = true;

   for (unsigned int i = 0; i < s_tLimVec.size(); i++) {
      if (event.getArrTime() < s_tLimVec[i].first ||
          event.getArrTime() > s_tLimVec[i].second) acceptEvent = false;
   }

   if (event.getEnergy() < s_eMin || event.getEnergy() > s_eMax) 
      acceptEvent = false;

   double dist = event.getSeparation(s_roiCenter)*180./M_PI;
   if (dist > s_roiRadius) {
//      std::cerr << dist << std::endl;
      acceptEvent = false;
   }

   if (event.getMuZenith() < s_muZenMax) acceptEvent = false;

   return acceptEvent;
}

RoiCuts * RoiCuts::instance() {
   if (s_instance == 0) {
      s_instance = new RoiCuts();
   }
   return s_instance;
}

} // namespace Likelihood
