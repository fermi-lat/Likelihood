/** 
 * @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Region-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RoiCuts.cxx,v 1.26 2004/12/08 21:44:11 jchiang Exp $
 */

#include <cstdlib>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/dom/DOM.hpp>

#include "xml/Dom.h"
#include "xml/XmlParser.h"

#include "facilities/Util.h"

#include "tip/Header.h"

#include "optimizers/Dom.h"

#include "dataSubselector/Gti.h"

#include "Likelihood/Exception.h"
#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"

namespace Likelihood {

XERCES_CPP_NAMESPACE_USE

std::vector<RoiCuts::timeInterval> RoiCuts::s_tLimVec;
double RoiCuts::s_eMin;
double RoiCuts::s_eMax;
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
   s_tLimVec.push_back(std::make_pair(tmin, tmax));
}

void RoiCuts::setCuts(double ra, double dec, double roi_radius,
                      double emin, double emax,
                      double tmin, double tmax,
                      double muZenMax) {
   s_tLimVec.clear();
   addTimeInterval(tmin, tmax);
    
// min and max energies in MeV.
   s_eMin = emin;
   s_eMax = emax;
        
   s_roiCone = irfInterface::AcceptanceCone(astro::SkyDir(ra, dec),
                                            roi_radius);
   s_muZenMax = muZenMax;
}

void RoiCuts::setCuts(std::string xmlFile) {
        
   facilities::Util::expandEnvVar(&xmlFile);
        
   xml::XmlParser * parser = new xml::XmlParser();
        
   DOMDocument * doc = parser->parse(xmlFile.c_str());
        
   if (doc == 0) { // xml file not parsed successfully
      std::string errorMessage = "RoiCuts::setCuts:\nInput xml file, "
         + xmlFile + " not parsed successfully.";
      throw Exception(errorMessage);
   }

// Direct Xerces API call.        
   DOMElement * roi = doc->getDocumentElement();
   if (!xml::Dom::checkTagName(roi, "Region-of-Interest")) {
      throw Exception(std::string("RoiCuts::setCuts:\n")
                      + "Region-of-Interest root element not found in "
                      + xmlFile);
   }
        
// Read in time intervals.
   s_tLimVec.clear();
   std::vector<DOMElement *> times;
   xml::Dom::getChildrenByTagName(roi, "timeInterval", times);
   
   std::vector<DOMElement *>::const_iterator timeIt = times.begin();
   for ( ; timeIt != times.end(); timeIt++) {
      double start 
         = std::atof(xml::Dom::getAttribute(*timeIt, "start").c_str());
      double stop = std::atof(xml::Dom::getAttribute(*timeIt, "stop").c_str());
      if (xml::Dom::getAttribute(*timeIt, "unit") == "days") { 
// convert to seconds
         start *= 8.64e4;
         stop *= 8.64e4;
      }
      addTimeInterval(start, stop);
   }
   
// Energy interval.
   std::vector<DOMElement *> child;
   xml::Dom::getChildrenByTagName(roi, "energies", child);
   s_eMin = std::atof(xml::Dom::getAttribute(child[0], "emin").c_str());
   s_eMax = std::atof(xml::Dom::getAttribute(child[0], "emax").c_str());
   if (xml::Dom::getAttribute(child[0], "unit") == "GeV") {
// Convert to MeV.
      s_eMin *= 1e3;
      s_eMax *= 1e3;
   }

   xml::Dom::getChildrenByTagName(roi, "acceptanceCone", child);
   astro::SkyDir roiCenter;
   double lon 
      = std::atof(xml::Dom::getAttribute(child[0], "longitude").c_str());
   double lat 
      = std::atof(xml::Dom::getAttribute(child[0], "latitude").c_str());
   if (xml::Dom::getAttribute(child[0], "coordsys") 
       == std::string("Galactic")) {
      roiCenter = astro::SkyDir(lon, lat, astro::SkyDir::GALACTIC);
   } else {
      roiCenter = astro::SkyDir(lon, lat, astro::SkyDir::EQUATORIAL);
   }
   double roiRadius 
      = atof(xml::Dom::getAttribute(child[0], "radius").c_str());

   s_roiCone = irfInterface::AcceptanceCone(roiCenter, roiRadius);

// Do not apply zenith angle cut for now.
   s_muZenMax = -1.;

   delete parser;
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

void RoiCuts::writeGtiExtension(const std::string & filename) {
   if (s_cuts) {
      s_cuts->writeGtiExtension(filename);
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
   s_tLimVec.clear();
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

void RoiCuts::writeXml(std::string xmlFile, const std::string &roiTitle) {
   facilities::Util::expandEnvVar(&xmlFile);
   std::ofstream outFile(xmlFile.c_str());
   writeXml(outFile, roiTitle, true);
}

void RoiCuts::writeXml(std::ostream & ostr, const std::string & roiTitle,
                       bool pretty) {
   DOMElement * roiElt = rootDomElement(roiTitle);
   if (pretty) {
      ostr << "<?xml version='1.0' standalone='no'?>\n"
           << "<!DOCTYPE Region-of-Interest SYSTEM "
           << "\"$(LIKELIHOODROOT)/xml/RoiCuts.dtd\" >\n";
      xml::Dom::prettyPrintElement(roiElt, ostr, "");
   } else {
      ostr << "<?xml version='1.0' standalone='no'?>"
           << "<!DOCTYPE Region-of-Interest SYSTEM "
           << "\"$(LIKELIHOODROOT)/xml/RoiCuts.dtd\" >";
      xml::Dom::printElement(roiElt, ostr);
   }
   roiElt->release();
}

bool RoiCuts::accept(const Event &event) {
   if (s_cuts) {
      std::map<std::string, double> params;
      params["TIME"] = event.getArrTime();
      params["ENERGY"] = event.getEnergy();
      params["RA"] = event.getDir().ra();
      params["DEC"] = event.getDir().dec();
      return s_cuts->accept(params);
   } else {
      bool acceptEvent = true;
      
      for (unsigned int i = 0; i < s_tLimVec.size(); i++) {
         if (event.getArrTime() < s_tLimVec[i].first ||
             event.getArrTime() > s_tLimVec[i].second) {
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

DOMElement * RoiCuts::rootDomElement(const std::string &roiTitle) {

   xml::XmlParser * parser = new xml::XmlParser();

   DOMDocument * doc = optimizers::Dom::createDocument();

   DOMElement * roiElt = optimizers::Dom::createElement(doc,
                                                        "Region-of-Interest");
   xml::Dom::addAttribute(roiElt, "title", roiTitle.c_str());

// Loop over time intervals
   std::vector<timeInterval>::iterator tintIt = s_tLimVec.begin();
   for ( ; tintIt != s_tLimVec.end(); tintIt++) {
      DOMElement * tintElt = optimizers::Dom::createElement(doc,
                                                            "timeInterval");
      xml::Dom::addAttribute(tintElt, std::string("start"), tintIt->first);
      xml::Dom::addAttribute(tintElt, std::string("stop"), tintIt->second);
      xml::Dom::addAttribute(tintElt, "unit", "seconds");
      optimizers::Dom::appendChild(roiElt, tintElt);
   }

   DOMElement * energElt = optimizers::Dom::createElement(doc, "energies");
   xml::Dom::addAttribute(energElt, std::string("emin"), s_eMin);
   xml::Dom::addAttribute(energElt, std::string("emax"), s_eMax);
   xml::Dom::addAttribute(energElt, "unit", "MeV");
   optimizers::Dom::appendChild(roiElt, energElt);

   DOMElement * coneElt = optimizers::Dom::createElement(doc,
                                                         "acceptanceCone");
   xml::Dom::addAttribute(coneElt, std::string("longitude"), 
                          s_roiCone.center().ra());
   xml::Dom::addAttribute(coneElt, std::string("latitude"), 
                          s_roiCone.center().dec());
   xml::Dom::addAttribute(coneElt, "coordsys", "J2000");
   xml::Dom::addAttribute(coneElt, std::string("radius"), 
                          s_roiCone.radius());

   optimizers::Dom::appendChild(roiElt, coneElt);
   delete parser;

   return roiElt;
}

} // namespace Likelihood
