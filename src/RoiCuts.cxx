/** 
 * @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Region-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RoiCuts.cxx,v 1.7 2003/05/29 20:10:46 jchiang Exp $
 */

#include <sstream>

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>

#include "optimizers/Dom.h"

#include "Likelihood/Exception.h"
#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"

namespace Likelihood {

// definitions of static data
std::vector<RoiCuts::timeInterval> RoiCuts::s_tLimVec;
double RoiCuts::s_eMin;
double RoiCuts::s_eMax;
astro::SkyDir RoiCuts::s_roiCenter;
double RoiCuts::s_roiRadius;
double RoiCuts::s_muZenMax;
RoiCuts * RoiCuts::s_instance = 0;

void RoiCuts::setCuts(double ra, double dec, double roi_radius) {
// get everything for now....
   s_tLimVec.push_back(std::make_pair(0., HUGE));

// default min and max energies in MeV
   s_eMin = 31.623;
   s_eMax = 3.1622e5;

// // prompt the user for sky extraction region info
//    double ra;
//    double dec;
//    std::cout << "Enter ra and dec (in degrees) of ROI center: ";
//    std::cin >> ra >> dec;
    s_roiCenter = astro::SkyDir(ra, dec);

//    std::cout << "Enter acceptance cone radius (in degrees): ";
//    std::cin >> s_roiRadius;
    s_roiRadius = roi_radius;

// accept everything until effect of Zenith angle cuts can be computed
   s_muZenMax = -1.;
}

void RoiCuts::setCuts(const std::string &xmlFile) {
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
      s_tLimVec.push_back(std::make_pair(start, stop));
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
   double lon = ::atof(xml::Dom::getAttribute(child[0], "longitude").c_str());
   double lat = ::atof(xml::Dom::getAttribute(child[0], "latitude").c_str());
   if (xml::Dom::getAttribute(child[0], "unit") == "Galactic") {
      s_roiCenter = astro::SkyDir(lon, lat, astro::SkyDir::GALACTIC);
   } else {
      s_roiCenter = astro::SkyDir(lon, lat, astro::SkyDir::CELESTIAL);
   }
   s_roiRadius = ::atof(xml::Dom::getAttribute(child[0], "radius").c_str());

// Do not apply zenith angle cut for now.
   s_muZenMax = -1.;
}

bool RoiCuts::accept(const Event &event) {
   bool acceptEvent = true;

   for (unsigned int i = 0; i < s_tLimVec.size(); i++) {
      if (event.getArrTime() < s_tLimVec[i].first ||
          event.getArrTime() > s_tLimVec[i].second) acceptEvent = false;
   }

   if (event.getEnergy() < s_eMin || event.getEnergy() > s_eMax) 
      acceptEvent = false;

   double dist = event.getSeparation(s_roiCenter)*180./M_PI;
   if (dist > s_roiRadius) {
      acceptEvent = false;
   }

   if (event.getMuZenith() < s_muZenMax) acceptEvent = false;

   return acceptEvent;
}

RoiCuts * RoiCuts::instance() {
   if (s_instance == 0) {
      s_instance = new RoiCuts();
   }
   return s_instance;
}

} // namespace Likelihood
