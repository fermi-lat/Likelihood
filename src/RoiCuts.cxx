/** 
 * @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Region-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RoiCuts.cxx,v 1.13 2003/11/14 00:13:20 petry Exp $
 */

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include <xercesc/dom/DOM_Document.hpp>
#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>
#include <xercesc/dom/DOM_DOMException.hpp>

#include "facilities/Util.h"

#include "optimizers/Dom.h"

#include "Likelihood/Exception.h"
#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"

namespace Likelihood {

// Definitions of static data.
std::vector<RoiCuts::timeInterval> RoiCuts::s_tLimVec;
double RoiCuts::s_eMin;
double RoiCuts::s_eMax;
latResponse::AcceptanceCone RoiCuts::s_roiCone;
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
        
   s_roiCone = latResponse::AcceptanceCone(astro::SkyDir(ra, dec),
                                           roi_radius);
   s_muZenMax = muZenMax;
}
    

void RoiCuts::setCuts(std::string xmlFile) {
        
// Expand any environment variables in the xmlFile name.
   facilities::Util::expandEnvVar(&xmlFile);
        
   xml::XmlParser *parser = new xml::XmlParser();
        
   DOM_Document doc = parser->parse(xmlFile.c_str());
        
   if (doc == 0) { // xml file not parsed successfully
      std::ostringstream errorMessage;
      errorMessage << "RoiCuts::setCuts: "
                   << "Input xml file, " << xmlFile 
                   << " not parsed successfully.\n";
      throw Exception(errorMessage.str());
   }
        
   DOM_Element roi = doc.getDocumentElement();
   optimizers::Dom::checkTag(roi, "Region-of-Interest",
                             "RoiCuts::setCuts");
        
// Read in time intervals.
   s_tLimVec.clear();
   std::vector<DOM_Element> times;
   optimizers::Dom::getElements(roi, "timeInterval", times);
   
   std::vector<DOM_Element>::const_iterator timeIt = times.begin();
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
   std::vector<DOM_Element> child;
   optimizers::Dom::getElements(roi, "energies", child);
   s_eMin = ::atof(xml::Dom::getAttribute(child[0], "emin").c_str());
   s_eMax = ::atof(xml::Dom::getAttribute(child[0], "emax").c_str());
   if (xml::Dom::getAttribute(child[0], "unit") == "GeV") {
// convert to MeV
      s_eMin *= 1e3;
      s_eMax *= 1e3;
   }

   optimizers::Dom::getElements(roi, "acceptanceCone", child);
   astro::SkyDir roiCenter;
   double lon = ::atof(xml::Dom::getAttribute(child[0], "longitude").c_str());
   double lat = ::atof(xml::Dom::getAttribute(child[0], "latitude").c_str());
   if (xml::Dom::getAttribute(child[0], "coordsys") 
       == std::string("Galactic")) {
      roiCenter = astro::SkyDir(lon, lat, astro::SkyDir::GALACTIC);
   } else {
      roiCenter = astro::SkyDir(lon, lat, astro::SkyDir::CELESTIAL);
   }
   double roiRadius 
      = atof(xml::Dom::getAttribute(child[0], "radius").c_str());

   s_roiCone = latResponse::AcceptanceCone(roiCenter, roiRadius);

// Do not apply zenith angle cut for now.
   s_muZenMax = -1.;

   delete parser;
}

void RoiCuts::writeXml(std::string xmlFile, const std::string &roiTitle) {

   xml::XmlParser *parser = new xml::XmlParser();

   DOM_Document doc = DOM_Document::createDocument();

   DOM_Element roiElt = doc.createElement("Region-of-Interest");
   roiElt.setAttribute("title", roiTitle.c_str());

// Loop over time intervals
   std::vector<timeInterval>::iterator tintIt = s_tLimVec.begin();
   for ( ; tintIt != s_tLimVec.end(); tintIt++) {
      DOM_Element tintElt = doc.createElement("timeInterval");
      std::ostringstream startTime;
      startTime << tintIt->first;
      tintElt.setAttribute("start", startTime.str().c_str());
      std::ostringstream stopTime;
      stopTime << tintIt->second;
      tintElt.setAttribute("stop", stopTime.str().c_str());
      tintElt.setAttribute("unit", "seconds");
      roiElt.appendChild(tintElt);
   }

   DOM_Element energElt = doc.createElement("energies");
   std::ostringstream eMin;
   eMin << s_eMin;
   energElt.setAttribute("emin", eMin.str().c_str());
   std::ostringstream eMax;
   eMax << s_eMax;
   energElt.setAttribute("emax", eMax.str().c_str());
   energElt.setAttribute("unit", "MeV");
   roiElt.appendChild(energElt);

   DOM_Element coneElt = doc.createElement("acceptanceCone");
   std::ostringstream ra;
   ra << s_roiCone.center().ra();
   coneElt.setAttribute("ra", ra.str().c_str());
   std::ostringstream dec;
   dec << s_roiCone.center().dec();
   coneElt.setAttribute("dec", dec.str().c_str());
   coneElt.setAttribute("coordsys", "Equatorial");
   std::ostringstream roi_radius;
   roi_radius << s_roiCone.radius();
   coneElt.setAttribute("radius", roi_radius.str().c_str()); 
   roiElt.appendChild(coneElt);
            
// Expand any environment variables in the xmlFile name.
   facilities::Util::expandEnvVar(&xmlFile);

   std::ofstream outFile(xmlFile.c_str());
//    outFile << "<?xml version='1.0' standalone='no'?>\n"
//            << "<!DOCTYPE Region-of-Interest SYSTEM \"RoiCuts.dtd\" >\n";
   xml::Dom::prettyPrintElement(roiElt, outFile, "");

   delete parser;
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
