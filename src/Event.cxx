/** @file Event.cxx
 * @brief Event class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Event.cxx,v 1.6 2003/03/25 23:22:03 jchiang Exp $
 */

#include "Likelihood/Event.h"

namespace Likelihood {

Event::Event(double ra, double dec, double energy, 
             double time, double sc_ra, double sc_dec, double muZenith) {
   
   m_appDir = astro::SkyDir(ra, dec);
   m_energy = energy;
   m_arrTime = time;
   m_scDir = astro::SkyDir(sc_ra, sc_dec);
   m_muZenith = muZenith;

//   m_computeResponse();
}

Event::Event(const Event &event) {
   m_appDir = event.m_appDir;
   m_energy = event.m_energy;
   m_arrTime = event.m_arrTime;
   m_scDir = event.m_scDir;
   m_respEg = event.m_respEg;
   m_respGal = event.m_respGal;
   m_respDiffuseSrcs = event.m_respDiffuseSrcs;
}

double Event::diffuseResponse(double energy, 
                              std::string diffuseComponent) const {
   if (m_respDiffuseSrcs.count(diffuseComponent)) {

// Since the energy resolution is presently assumed to be infinite,
// simply return the (second member of the pair of the) first (and
// only) element of the diffuse_response vector.

// Avoid operator[], so use iterator to respect const-ness of map
      std::map<std::string, diffuse_response>::const_iterator it 
         = m_respDiffuseSrcs.begin();
      for (; it != m_respDiffuseSrcs.end(); it++) {
         if (it->first == diffuseComponent) {
            return (*it).second[0].second;
         }
      }
   } else {
      std::cerr << "Event::diffuseResponse: Diffuse component " 
                << diffuseComponent 
                << " does not have an associated diffuse response."
                << std::endl;
      assert(false);
      return 0;
   }
   return 0;
}

} // namespace Likelihood
