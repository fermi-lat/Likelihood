/** @file Event.cxx
 * @brief Event class implementation
 *
 * $Header:
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

} // namespace Likelihood
