/** 
 * @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Region-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RoiCuts.cxx,v 1.32 2005/03/01 07:17:07 jchiang Exp $
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

std::vector<RoiCuts::timeInterval> RoiCuts::s_timeCuts;
double RoiCuts::s_eMin(20.);
double RoiCuts::s_eMax(2e5);
irfInterface::AcceptanceCone RoiCuts::s_roiCone;
double RoiCuts::s_muZenMax;

dataSubselector::Cuts * RoiCuts::s_cuts(0);

RoiCuts * RoiCuts::s_instance(0);

RoiCuts * RoiCuts::instance() {
   if (s_instance == 0) {
      s_instance = new RoiCuts();
   }
   return s_instance;
}

void RoiCuts::addTimeInterval(double tmin, double tmax) {
   s_timeCuts.push_back(std::make_pair(tmin, tmax));
}

void RoiCuts::setCuts(double ra, double dec, double roi_radius,
                      double emin, double emax,
                      double tmin, double tmax,
                      double muZenMax) {
   s_timeCuts.clear();
   addTimeInterval(tmin, tmax);
    
// min and max energies in MeV.
   s_eMin = emin;
   s_eMax = emax;

   makeEnergyVector();
        
   s_roiCone = irfInterface::AcceptanceCone(astro::SkyDir(ra, dec),
                                            roi_radius);
   s_muZenMax = muZenMax;
}

void RoiCuts::readCuts(const std::string & eventFile, 
                       const std::string & ext, bool strict) {
   s_cuts = new dataSubselector::Cuts(eventFile, ext, false);
   s_instance->sortCuts(strict);
   s_instance->setRoiData();
}

void RoiCuts::writeDssKeywords(tip::Header & header) const {
   if (s_cuts) {
      s_cuts->writeDssKeywords(header);
   }
}

void RoiCuts::writeGtiExtension(const std::string & filename) const {
   if (s_cuts) {
      s_cuts->writeGtiExtension(filename);
   }
}

void RoiCuts::makeEnergyVector(int nee) {
   m_energies.clear();
   m_energies.reserve(nee);
   double estep = std::log(s_eMax/s_eMin)/(nee - 1.);
   for (int i = 0; i < nee; i++) {
      m_energies.push_back(s_eMin*std::exp(estep*i));
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
   s_timeCuts.clear();
   for (unsigned int i = 0; i < m_timeCuts.size(); i++) {
      addTimeInterval(m_timeCuts.at(i)->minVal(), m_timeCuts.at(i)->maxVal());
   }
   for (unsigned int i = 0; i < m_gtiCuts.size(); i++) {
      const dataSubselector::Gti & gti = m_gtiCuts.at(i)->gti();
      evtbin::Gti::ConstIterator it;
      for (it = gti.begin(); it != gti.end(); ++it) {
         addTimeInterval(it->first, it->second);
      }
   }
}

void RoiCuts::sortCuts(bool strict) {
   typedef dataSubselector::CutBase CutBase;
   typedef dataSubselector::RangeCut RangeCut;
   typedef dataSubselector::GtiCut GtiCut;
   typedef dataSubselector::SkyConeCut SkyConeCut;

   unsigned int nenergy(0), ncone(0), ntime(0);
   for (unsigned int i = 0; i < s_cuts->size(); i++) {
      CutBase & cut = const_cast<CutBase &>(s_cuts->operator[](i));
      if (cut.type() == "range") {
         RangeCut & rangeCut = dynamic_cast<RangeCut &>(cut);
         std::string colname = rangeCut.colname();
         if (colname == "ENERGY") {
            nenergy++;
            m_energyCut = &rangeCut;
         } else if (colname == "TIME") {
            ntime++;
            m_timeCuts.push_back(&rangeCut);
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
   if (strict && (nenergy != 1 || ncone != 1 || ntime == 0)) {
      std::ostringstream message;
      message << "RoiCuts::sortCuts:\n"
              << "There should be exactly one energy range cut, "
              << "one acceptance cone cut,\n"
              << "and at least one time range and/or GTI cut.\n"
              << "The event file contains the following DSS selections:\n\n";
      s_cuts->writeCuts(message);
      throw std::runtime_error(message.str());
   }
}

bool RoiCuts::accept(const Event &event) const {
   if (s_cuts) {
      std::map<std::string, double> params;
      params["TIME"] = event.getArrTime();
      params["ENERGY"] = event.getEnergy();
      params["RA"] = event.getDir().ra();
      params["DEC"] = event.getDir().dec();
      return s_cuts->accept(params);
   } else {
      bool acceptEvent = true;

/// @todo treat RangeCuts in time differently from GTIs.      
      for (unsigned int i = 0; i < s_timeCuts.size(); i++) {
         if (event.getArrTime() < s_timeCuts[i].first ||
             event.getArrTime() > s_timeCuts[i].second) {
            acceptEvent = false;
         }
      }

      if (event.getEnergy() < s_eMin || event.getEnergy() > s_eMax) { 
         acceptEvent = false;
      }

      double dist = event.getSeparation(s_roiCone.center())*180./M_PI;
      if (dist > s_roiCone.radius()) {
         acceptEvent = false;
      }

      if (event.getMuZenith() < s_muZenMax){
         acceptEvent = false;
      }

      return acceptEvent;
   }
   return false;
}

} // namespace Likelihood
