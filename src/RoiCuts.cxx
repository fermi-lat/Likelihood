/** 
 * @file RoiCuts.cxx
 * @brief Implementation for RoiCuts, a Singleton class that contains
 * the Region-of-Interest cuts.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/RoiCuts.cxx,v 1.12 2003/11/07 02:27:10 jchiang Exp $
 */

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>

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
                          double coszmax) {
        s_tLimVec.clear();
        addTimeInterval(tmin, tmax);
    
//  min and max energies in MeV.
        s_eMin = emin;
        s_eMax = emax;
        
        s_roiCone = latResponse::AcceptanceCone(astro::SkyDir(ra, dec),
                                                roi_radius);
        s_muZenMax = coszmax;
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
        if (xml::Dom::getAttribute(child[0], "coordsys") 
            == std::string("Galactic")) {
            double lon = ::atof(xml::Dom::getAttribute(child[0], "longitude").c_str());
            double lat = ::atof(xml::Dom::getAttribute(child[0], "latitude").c_str());
            roiCenter = astro::SkyDir(lon, lat, astro::SkyDir::GALACTIC);
        } else {
            double lon = ::atof(xml::Dom::getAttribute(child[0], "ra").c_str());
            double lat = ::atof(xml::Dom::getAttribute(child[0], "dec").c_str());
            roiCenter = astro::SkyDir(lon, lat, astro::SkyDir::CELESTIAL);
        }
        double roiRadius 
            = atof(xml::Dom::getAttribute(child[0], "radius").c_str());

        s_roiCone = latResponse::AcceptanceCone(roiCenter, roiRadius);

// Do not apply zenith angle cut for now.
        s_muZenMax = -1.;

        delete parser;
    }


    void RoiCuts::writeXml(std::string xmlFile, std::string title){

        char auxString[80];

        xml::XmlParser *parser = new xml::XmlParser();

        DOM_Document doc = DOM_Document::createDocument();

        DOM_Element roiElt = doc.createElement("Region-of-Interest");
        roiElt.setAttribute("title", title.c_str());

// Loop over time intervals
        std::vector<timeInterval>::iterator tintIt = s_tLimVec.begin();
        for ( ; tintIt != s_tLimVec.end(); tintIt++) {
            DOM_Element tintElt = doc.createElement("timeInterval");
            sprintf(auxString, "%.8g", tintIt->first);
            tintElt.setAttribute("start", auxString);
            sprintf(auxString, "%.8g", tintIt->second);
            tintElt.setAttribute("stop", auxString);
            tintElt.setAttribute("unit", "seconds");
            roiElt.appendChild(tintElt);
        }

        DOM_Element energElt = doc.createElement("energies");
        sprintf(auxString, "%.8g", s_eMin);
        energElt.setAttribute("emin", auxString);
        sprintf(auxString, "%.8g", s_eMax);
        energElt.setAttribute("emax", auxString);
        energElt.setAttribute("unit", "MeV");
        roiElt.appendChild(energElt);

        DOM_Element coneElt = doc.createElement("acceptanceCone");
        sprintf(auxString, "%.8g", s_roiCone.center().ra());
        coneElt.setAttribute("ra", auxString);
        sprintf(auxString, "%.8g", s_roiCone.center().dec());
        coneElt.setAttribute("dec", auxString);
        coneElt.setAttribute("coordsys", "Equatorial");
        sprintf(auxString, "%.8g", s_roiCone.radius());
        coneElt.setAttribute("radius", auxString); 
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
//                cout << "Event is outside time range " << i << ".\n";
            }
        }

        if (event.getEnergy() < s_eMin || event.getEnergy() > s_eMax){ 
            acceptEvent = false;
//            cout << "Event is outside energy range.\n";
        }

        double dist = event.getSeparation(s_roiCone.center())*180./M_PI;
        if (dist > s_roiCone.radius()) {
            acceptEvent = false;
//            cout << "Event is outside cone.\n";
        }

        if (event.getMuZenith() < s_muZenMax){
            acceptEvent = false;
//            cout << "Event is outside ZA range.\n";
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
