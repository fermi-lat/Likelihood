#include "Likelihood/RoiCuts.h"

namespace Likelihood {

// definitions of static data
std::vector<RoiCuts::timeInterval> RoiCuts::m_tLimVec;
double RoiCuts::m_eMin;
double RoiCuts::m_eMax;
astro::SkyDir RoiCuts::m_roiCenter;
double RoiCuts::m_roiRadius;
double RoiCuts::m_muZenMax;
RoiCuts * RoiCuts::s_instance = 0;

void RoiCuts::setCuts(double ra, double dec, double roi_radius) {
// get everything for now....
    m_tLimVec.push_back(std::make_pair(0., HUGE));

// default min and max energies in MeV
   m_eMin = 31.623;
   m_eMax = 3.1622e5;

// // prompt the user for sky extraction region info
//    double ra;
//    double dec;
//    std::cout << "Enter ra and dec (in degrees) of ROI center: ";
//    std::cin >> ra >> dec;
    m_roiCenter = astro::SkyDir(ra, dec);

//    std::cout << "Enter acceptance cone radius (in degrees): ";
//    std::cin >> m_roiRadius;
    m_roiRadius = roi_radius;

// accept everything until effect of Zenith angle cuts can be computed
   m_muZenMax = -1.;
}

bool RoiCuts::accept(const Event &event) {
   bool acceptEvent = true;

   for (unsigned int i = 0; i < m_tLimVec.size(); i++) {
      if (event.getArrTime() < m_tLimVec[i].first ||
          event.getArrTime() > m_tLimVec[i].second) acceptEvent = false;
   }

   if (event.getEnergy() < m_eMin || event.getEnergy() > m_eMax) 
      acceptEvent = false;

   double dist = event.getSeparation(m_roiCenter)*180./M_PI;
   if (dist > m_roiRadius) {
//      std::cerr << dist << std::endl;
      acceptEvent = false;
   }

   if (event.getMuZenith() < m_muZenMax) acceptEvent = false;

   return acceptEvent;
}

RoiCuts * RoiCuts::instance() {
   if (s_instance == 0) {
      s_instance = new RoiCuts();
   }
   return s_instance;
}

} // namespace Likelihood
