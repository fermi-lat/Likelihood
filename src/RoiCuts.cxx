/** 
 * @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Region-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RoiCuts.cxx,v 1.36 2005/03/22 00:18:13 jchiang Exp $
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
   m_timeCuts.push_back(std::make_pair(tmin, tmax));
}

void RoiCuts::addGoodTimeInterval(double tmin, double tmax) {
   m_gtis.push_back(std::make_pair(tmin, tmax));
}

void RoiCuts::setCuts(double ra, double dec, double roi_radius,
                      double emin, double emax,
                      double tmin, double tmax,
                      double muZenMax) {
   m_timeCuts.clear();
   addTimeInterval(tmin, tmax);
   m_minTime = tmin;
   m_maxTime = tmax;

   m_gtis.clear();
    
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

void RoiCuts::writeDssKeywords(tip::Header & header) const {
   if (m_cuts) {
      m_cuts->writeDssKeywords(header);
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
   m_timeCuts.clear();
   for (unsigned int i = 0; i < m_timeRangeCuts.size(); i++) {
      addTimeInterval(m_timeRangeCuts.at(i)->minVal(),
                      m_timeRangeCuts.at(i)->maxVal());
      if (i == 0 || m_timeRangeCuts.at(i)->minVal() < m_minTime) {
         m_minTime = m_timeRangeCuts.at(i)->minVal();
      }
      if (i == 0 || m_timeRangeCuts.at(i)->maxVal() > m_maxTime) {
         m_maxTime = m_timeRangeCuts.at(i)->maxVal();
      }
   }
   for (unsigned int i = 0; i < m_gtiCuts.size(); i++) {
      const dataSubselector::Gti & gti = m_gtiCuts.at(i)->gti();
      evtbin::Gti::ConstIterator it;
      for (it = gti.begin(); it != gti.end(); ++it) {
         addGoodTimeInterval(it->first, it->second);
      }
   }
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
//   if (strict && (nenergy != 1 || ncone != 1 || ntime == 0)) {
   if (strict && (ncone != 1 || ntime == 0)) {
      std::ostringstream message;
      message << "RoiCuts::sortCuts:\n"
//              << "There should be exactly one energy range cut, "
              << "There should be exactly "
              << "one acceptance cone cut,\n"
              << "and at least one time range and/or GTI cut.\n"
              << "The event file contains the following DSS selections:\n\n";
      m_cuts->writeCuts(message);
      throw std::runtime_error(message.str());
   }
}

bool RoiCuts::accept(const Event &event) const {
//    if (m_cuts) {
//       std::map<std::string, double> params;
//       params["TIME"] = event.getArrTime();
//       params["ENERGY"] = event.getEnergy();
//       params["RA"] = event.getDir().ra();
//       params["DEC"] = event.getDir().dec();
//       return m_cuts->accept(params);
//    } else {
      bool acceptEvent(false);

      if (m_gtis.size() == 0) {
         acceptEvent = true;
      } else {
// Accept the event if it appears in any of the GTIs.
         for (unsigned int i = 0; i < m_gtis.size(); i++) {
            if (m_gtis.at(i).first <= event.getArrTime() &&
                event.getArrTime() <= m_gtis.at(i).second) {
               acceptEvent = true;
               break;
            }
         }
      }

// Require the event to lie within all time range cuts.
      for (unsigned int i = 0; i < m_timeCuts.size(); i++) {
         if (event.getArrTime() < m_timeCuts[i].first ||
             event.getArrTime() > m_timeCuts[i].second) {
            acceptEvent = false;
         }
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
//    }
//    return false;
}

} // namespace Likelihood
