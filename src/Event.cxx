/** @file Event.cxx
 * @brief Event class implementation
 *
 * $Header:
 */

#include "../Likelihood/Event.h"

namespace Likelihood {

Event::Event(double ra, double dec, double energy, 
             double time, double sc_ra, double sc_dec, double muZenith) {
   
//! astro::SkyDir needs a setDir(ra, dec) method
   astro::SkyDir appDir(ra, dec);
   m_appDir = appDir;
   m_energy = energy;
   m_arrTime = time;
   astro::SkyDir scDir(sc_ra, sc_dec);
   m_scDir = scDir;
   m_muZenith = muZenith;

//   m_computeResponse();
}

//! copy constructor
Event::Event(const Event &event) {
   m_appDir = event.m_appDir;
   m_energy = event.m_energy;
   m_arrTime = event.m_arrTime;
   m_scDir = event.m_scDir;
   m_respEg = event.m_respEg;
   m_respGal = event.m_respGal;
   m_respDiffuseSrcs = event.m_respDiffuseSrcs;
}

//  void Event::setDir(double ra, double dec) {
//     astro::SkyDir appDir(ra, dec);
//     m_appDir = appDir;
//  }   

} // namespace Likelihood
