/** 
 * @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Region-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RoiCuts.cxx,v 1.19 2004/07/19 14:16:58 jchiang Exp $
 */

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>

#include "xml/Dom.h"
#include "xml/XmlParser.h"

#include "facilities/Util.h"

#include "Likelihood/Exception.h"
#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"

namespace {
   DomDocument * createDocument() {
      DomDocument * doc = new DOM_Document();
      *doc = DOM_Document::createDocument();
      return doc;
   }
   DomElement * createElement(DomDocument * doc, const std::string & name) {
      DomElement * elt = new DOM_Element();
      *elt = doc->createElement(name.c_str());
      return elt;
   }
}

namespace Likelihood {

// Definitions of static data.
std::vector<RoiCuts::timeInterval> RoiCuts::s_tLimVec;
double RoiCuts::s_eMin;
double RoiCuts::s_eMax;
irfInterface::AcceptanceCone RoiCuts::s_roiCone;
double RoiCuts::s_muZenMax;
RoiCuts * RoiCuts::s_instance = 0;

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
        
   xml::XmlParser *parser = new xml::XmlParser();
        
   DomDocument doc = parser->parse(xmlFile.c_str());
        
   if (doc == 0) { // xml file not parsed successfully
      std::string errorMessage = "RoiCuts::setCuts:\nInput xml file, "
         + xmlFile + " not parsed successfully.";
      throw Exception(errorMessage);
   }

// Direct Xerces API call.        
   DomElement roi = doc.getDocumentElement();
   if (!xml::Dom::checkTagName(roi, "Region-of-Interest")) {
      throw Exception(std::string("RoiCuts::setCuts:\n")
                      + "Region-of-Interest root element not found in "
                      + xmlFile);
   }
        
// Read in time intervals.
   s_tLimVec.clear();
   std::vector<DomElement> times;
   xml::Dom::getChildrenByTagName(roi, "timeInterval", times);
   
   std::vector<DomElement>::const_iterator timeIt = times.begin();
   for ( ; timeIt != times.end(); timeIt++) {
      double start = ::atof(xml::Dom::getAttribute(*timeIt, "start").c_str());
      double stop = ::atof(xml::Dom::getAttribute(*timeIt, "stop").c_str());
      if (xml::Dom::getAttribute(*timeIt, "unit") == "days") { 
// convert to seconds
         start *= 8.64e4;
         stop *= 8.64e4;
      }
      addTimeInterval(start, stop);
   }
   
// Energy interval.
   std::vector<DomElement> child;
   xml::Dom::getChildrenByTagName(roi, "energies", child);
   s_eMin = ::atof(xml::Dom::getAttribute(child[0], "emin").c_str());
   s_eMax = ::atof(xml::Dom::getAttribute(child[0], "emax").c_str());
   if (xml::Dom::getAttribute(child[0], "unit") == "GeV") {
// Convert to MeV.
      s_eMin *= 1e3;
      s_eMax *= 1e3;
   }

   xml::Dom::getChildrenByTagName(roi, "acceptanceCone", child);
   astro::SkyDir roiCenter;
   double lon = ::atof(xml::Dom::getAttribute(child[0], "longitude").c_str());
   double lat = ::atof(xml::Dom::getAttribute(child[0], "latitude").c_str());
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

DomElement * RoiCuts::rootDomElement(const std::string &roiTitle) {

   xml::XmlParser *parser = new xml::XmlParser();

   DomDocument * doc = ::createDocument();

   DomElement * roiElt = ::createElement(doc, "Region-of-Interest");
   xml::Dom::addAttribute(*roiElt, "title", roiTitle.c_str());

// Loop over time intervals
   std::vector<timeInterval>::iterator tintIt = s_tLimVec.begin();
   for ( ; tintIt != s_tLimVec.end(); tintIt++) {
      DomElement * tintElt = ::createElement(doc, "timeInterval");
      xml::Dom::addAttribute(*tintElt, std::string("start"), tintIt->first);
      xml::Dom::addAttribute(*tintElt, std::string("stop"), tintIt->second);
      xml::Dom::addAttribute(*tintElt, "unit", "seconds");
      roiElt->appendChild(*tintElt);
   }

   DomElement * energElt = ::createElement(doc, "energies");
   xml::Dom::addAttribute(*energElt, std::string("emin"), s_eMin);
   xml::Dom::addAttribute(*energElt, std::string("emax"), s_eMax);
   xml::Dom::addAttribute(*energElt, "unit", "MeV");
   roiElt->appendChild(*energElt);

   DomElement * coneElt = ::createElement(doc, "acceptanceCone");
   xml::Dom::addAttribute(*coneElt, std::string("longitude"), 
                          s_roiCone.center().ra());
   xml::Dom::addAttribute(*coneElt, std::string("latitude"), 
                          s_roiCone.center().dec());
   xml::Dom::addAttribute(*coneElt, "coordsys", "J2000");
   xml::Dom::addAttribute(*coneElt, std::string("radius"), 
                          s_roiCone.radius());

   roiElt->appendChild(*coneElt);
   delete parser;

   return roiElt;
}

void RoiCuts::writeXml(std::string xmlFile, const std::string &roiTitle) {
   facilities::Util::expandEnvVar(&xmlFile);
   std::ofstream outFile(xmlFile.c_str());
   writeXml(outFile, roiTitle, true);
}

void RoiCuts::writeXml(std::ostream & ostr, const std::string & roiTitle,
                       bool pretty) {
   DomElement * roiElt = rootDomElement(roiTitle);
   if (pretty) {
      ostr << "<?xml version='1.0' standalone='no'?>\n"
           << "<!DOCTYPE Region-of-Interest SYSTEM "
           << "\"$(LIKELIHOODROOT)/xml/RoiCuts.dtd\" >\n";
      xml::Dom::prettyPrintElement(*roiElt, ostr, "");
   } else {
      ostr << "<?xml version='1.0' standalone='no'?>"
           << "<!DOCTYPE Region-of-Interest SYSTEM "
           << "\"$(LIKELIHOODROOT)/xml/RoiCuts.dtd\" >";
      xml::Dom::printElement(*roiElt, ostr);
   }
}

bool RoiCuts::accept(const Event &event) {
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

RoiCuts * RoiCuts::instance() {
   if (s_instance == 0) {
      s_instance = new RoiCuts();
   }
   return s_instance;
}

} // namespace Likelihood
