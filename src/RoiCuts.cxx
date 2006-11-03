/** 
 * @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Region-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RoiCuts.cxx,v 1.48 2006/09/20 23:43:38 jchiang Exp $
 */

#include <cstdlib>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "facilities/Util.h"

#include "tip/Header.h"

#include "dataSubselector/Gti.h"

#include "Likelihood/Exception.h"
#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"

namespace Likelihood {

void RoiCuts::addTimeInterval(double tmin, double tmax) {
   m_timeRangeCuts.push_back(new dataSubselector::RangeCut("TIME", "s",
                                                           tmin, tmax));
}

void RoiCuts::
getTimeCuts(std::vector< std::pair<double, double> > & timeCuts) const {
   timeCuts.clear();
   for (size_t i = 0; i < m_timeRangeCuts.size(); i++) {
      timeCuts.push_back(std::make_pair(m_timeRangeCuts.at(i)->minVal(),
                                        m_timeRangeCuts.at(i)->maxVal()));
   }
}

void RoiCuts::getGtis(std::vector< std::pair<double, double> > & gtis) const {
   gtis.clear();
// There should be at most one GTI extension.
   if (!m_gtiCuts.empty()) {
      const dataSubselector::Gti & gti(m_gtiCuts.at(0)->gti());
      evtbin::Gti::ConstIterator it(gti.begin());
      for ( ; it != gti.end(); ++it) {
         gtis.push_back(std::make_pair(it->first, it->second));
      }
   }
}

void RoiCuts::setCuts(double ra, double dec, double roi_radius,
                      double emin, double emax,
                      double tmin, double tmax,
                      double muZenMax, bool reset_tlims) {
   if (reset_tlims) {
      for (size_t i = 0; i < m_timeRangeCuts.size(); i++) {
         delete m_timeRangeCuts.at(i);
      }
      m_timeRangeCuts.clear();
      addTimeInterval(tmin, tmax);
   }
   m_minTime = tmin;
   m_maxTime = tmax;

   m_eMin = emin;
   m_eMax = emax;
   makeEnergyVector();

   m_roiCone = irfInterface::AcceptanceCone(astro::SkyDir(ra, dec),
                                            roi_radius);
   m_muZenMax = muZenMax;
}

void RoiCuts::readCuts(const std::string & eventFile, 
                       const std::string & ext, bool strict) {
   m_cuts = new dataSubselector::Cuts(eventFile, ext, false);
   sortCuts(strict);
   setRoiData();
}

void RoiCuts::readCuts(const std::vector<std::string> & eventFiles, 
                       const std::string & ext, bool strict) {
//   m_cuts = new dataSubselector::Cuts(eventFiles, ext, false, false, true);
   m_cuts = new dataSubselector::Cuts(eventFiles, ext, false, false, false);
   sortCuts(strict);
   setRoiData();
}

void RoiCuts::writeDssKeywords(tip::Header & header) const {
   if (m_cuts) {
      m_cuts->writeDssKeywords(header);
   }
}

void RoiCuts::writeDssTimeKeywords(tip::Header & header) const {
   if (m_cuts) {
      m_cuts->writeDssTimeKeywords(header);
   }
}

void RoiCuts::writeGtiExtension(const std::string & filename) const {
   if (m_cuts) {
      m_cuts->writeGtiExtension(filename);
   }
}

void RoiCuts::makeEnergyVector(int nee) {
   m_energies.clear();
   m_energies.reserve(nee);
   double estep = std::log(m_eMax/m_eMin)/(nee - 1.);
   for (int i = 0; i < nee; i++) {
      m_energies.push_back(m_eMin*std::exp(estep*i));
   }
}

void RoiCuts::setRoiData() {
   double ra(0);
   double dec(0);
   double radius(180.);
   if (m_skyConeCut) {
      ra = m_skyConeCut->ra();
      dec = m_skyConeCut->dec();
      radius = m_skyConeCut->radius();
   }
   double emin(20.);
   double emax(2e5);
   if (m_energyCut) {
      emin = m_energyCut->minVal();
      emax = m_energyCut->maxVal();
   }
   setCuts(ra, dec, radius, emin, emax);
}

void RoiCuts::sortCuts(bool strict) {
   typedef dataSubselector::CutBase CutBase;
   typedef dataSubselector::RangeCut RangeCut;
   typedef dataSubselector::GtiCut GtiCut;
   typedef dataSubselector::SkyConeCut SkyConeCut;

   unsigned int nenergy(0), ncone(0), ntime(0);
   for (unsigned int i = 0; i < m_cuts->size(); i++) {
      CutBase & cut = const_cast<CutBase &>(m_cuts->operator[](i));
      if (cut.type() == "range") {
         RangeCut & rangeCut = dynamic_cast<RangeCut &>(cut);
         std::string colname = rangeCut.colname();
         if (colname == "ENERGY") {
            if (nenergy == 0 || rangeCut.supercedes(*m_energyCut)) {
               nenergy++;
               m_energyCut = &rangeCut;
            }
         } else if (colname == "TIME") {
            ntime++;
            m_timeRangeCuts.push_back(&rangeCut);
         }
      }
      if (cut.type() == "SkyCone") {
         ncone++;
         m_skyConeCut = reinterpret_cast<SkyConeCut *>(&cut);
      }
      if (cut.type() == "GTI") {
         ntime++;
         m_gtiCuts.push_back(reinterpret_cast<GtiCut *>(&cut));
      }
   }
/// @todo Sort out the correct way to handle multiple, and perhaps
/// inconsistent, energy range and SkyCone cuts.
   if (strict && (nenergy != 1 || ncone != 1 || ntime == 0)) {
//   if (strict && (ncone != 1 || ntime == 0)) {
      std::ostringstream message;
      message << "RoiCuts::sortCuts:\n"
              << "There should be exactly "
              << "one acceptance cone cut,\n"
              << "and at least one time range and/or GTI cut.\n"
              << "The event file contains the following DSS selections:\n\n";
      bool suppressGtis(true);
      m_cuts->writeCuts(message, suppressGtis);
      throw std::runtime_error(message.str());
   }
}

bool RoiCuts::accept(const Event &event) const {
   bool acceptEvent(false);

   std::map<std::string, double> thisEvent;
   thisEvent["TIME"] = event.getArrTime();

   if (m_gtiCuts.size() == 0) {
      acceptEvent = true;
   } else {
// Accept the event if it appears in any of the GTIs.
      for (size_t i = 0; i < m_gtiCuts.size(); i++) {
         acceptEvent = acceptEvent || m_gtiCuts.at(i)->accept(thisEvent);
      }
   }

// Require the event to lie within all time range cuts.
   for (size_t i = 0; i < m_timeRangeCuts.size(); i++) {
      acceptEvent = m_timeRangeCuts.at(i)->accept(thisEvent);
   }

   if (event.getEnergy() < m_eMin || event.getEnergy() > m_eMax) { 
      acceptEvent = false;
   }

   double dist = event.getSeparation(m_roiCone.center())*180./M_PI;
   if (dist > m_roiCone.radius()) {
      acceptEvent = false;
   }

   if (event.getMuZenith() < m_muZenMax) {
      acceptEvent = false;
   }

   return acceptEvent;
}

} // namespace Likelihood
