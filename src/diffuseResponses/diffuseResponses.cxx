/**
 * @file diffuseResponses.cxx
 * @brief Adds diffuse response information for extragalactic and Galactic
 * diffuse emission.  
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/diffuseResponses/diffuseResponses.cxx,v 1.17 2004/12/30 16:45:44 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iostream>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SourceModel.h"
#include "Verbosity.h"

using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
using namespace Likelihood;

/**
 * @class diffuseResponses
 * @brief FTOOL to add diffuse response information to an FT1 file for
 * extragalactic and Galactic diffuse emission.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/diffuseResponses/diffuseResponses.cxx,v 1.17 2004/12/30 16:45:44 jchiang Exp $
 */

class diffuseResponses : public st_app::StApp {

public:

   diffuseResponses();

   virtual ~diffuseResponses() throw() {}

   virtual void run();

private:

   AppHelpers * m_helper;  //blech.
   SourceModel * m_srcModel;
   double m_srRadius;
   st_app::AppParGroup & m_pars;

   std::vector<Event> m_events;
   std::vector<DiffuseSource *> m_srcs;
   std::vector<std::string> m_srcNames;

   void promptForParameters();
   void readDiffuseNames(std::vector<std::string> & srcNames);
   bool haveDiffuseColumns();
   void buildSourceModel();
   void readEventData();
   void computeEventResponses();
   void writeEventResponses();
   void getDiffuseSources();
   void setGaussianParams(const Event & event, const std::string & name,
                          tip::Table::Vector<double> & params);
};

st_app::StAppFactory<diffuseResponses> myAppFactory;

diffuseResponses::diffuseResponses() 
   : st_app::StApp(), m_helper(0), m_srcModel(0), m_srRadius(30.),
     m_pars(st_app::StApp::getParGroup("diffuseResponses")) {}


void diffuseResponses::run() {
   promptForParameters();
   Likelihood::Verbosity::instance(m_pars["chatter"]);
   if (m_pars["clobber"] || !haveDiffuseColumns()) {
      m_helper = new AppHelpers(m_pars);
      m_helper->readScData();
      m_srcModel = new SourceModel(true);
      ResponseFunctions::setEdispFlag(m_pars["use_energy_dispersion"]);
      buildSourceModel();
      readEventData();
      computeEventResponses();
      writeEventResponses();
   } else {
      if (Likelihood::print_output()) {
         std::cout << "Diffuse columns have already been computed...exiting."
                   << std::endl;
      }
   }
}

void diffuseResponses::promptForParameters() {
   m_pars.Prompt();
   m_pars.Save();
}

void diffuseResponses::readDiffuseNames(std::vector<std::string> & srcNames) {
   srcNames.clear();
   std::string xmlFile = m_pars["source_model_file"];
   xmlBase::XmlParser * parser = new xmlBase::XmlParser();
   DOMDocument * doc = parser->parse(xmlFile.c_str());
   DOMElement * source_library = doc->getDocumentElement();
   std::vector<DOMElement *> srcs;
   xmlBase::Dom::getChildrenByTagName(source_library, "source", srcs);
   for (unsigned int i = 0; i < srcs.size(); i++) {
      if (xmlBase::Dom::getAttribute(srcs[i], "type") == "DiffuseSource") {
         srcNames.push_back(xmlBase::Dom::getAttribute(srcs[i], "name"));
      }
   }
}

bool diffuseResponses::haveDiffuseColumns() {
   std::auto_ptr<const tip::Table> 
      events(tip::IFileSvc::instance().readTable(m_pars["evfile"],
                                                 m_pars["evtable"]));
   const std::vector<std::string> & colNames = events->getValidFields();
   std::vector<std::string> srcNames;
   readDiffuseNames(srcNames);
   for (std::vector<std::string>::iterator name = srcNames.begin();
        name != srcNames.end(); ++name) {
      Event::toLower(*name);
      if (std::find(colNames.begin(), colNames.end(), *name) 
          == colNames.end()) {
         return false;
      }
   }
   return true;
}

void diffuseResponses::buildSourceModel() {
   std::string sourceModel = m_pars["source_model_file"];
   st_facilities::Util::file_ok(sourceModel);
   m_srcModel->readXml(sourceModel, m_helper->funcFactory(), false);
}

void diffuseResponses::readEventData() {
   std::string eventFile = m_pars["evfile"];
   facilities::Util::expandEnvVar(&eventFile);
   st_facilities::Util::file_ok(eventFile);
   const tip::Table * events 
      = tip::IFileSvc::instance().readTable(eventFile, m_pars["evtable"]);

   double ra;
   double dec;
   double energy;
   double time;
   double zenAngle;
   int convLayer;
   int eventType;

   ScData * scData = ScData::instance();

   tip::Table::ConstIterator it = events->begin();
   tip::ConstTableRecord & event = *it;
   for ( ; it != events->end(); ++it) {
      event["ra"].get(ra);
      event["dec"].get(dec);
      event["energy"].get(energy);
      event["time"].get(time);
      event["zenith_angle"].get(zenAngle);
      event["conversion_layer"].get(convLayer);
      if (convLayer < 12) { // Front
         eventType = 0;
      } else {
         eventType = 1;
      }
      Event thisEvent(ra, dec, energy, time, scData->zAxis(time), 
                      scData->xAxis(time), cos(zenAngle*M_PI/180.), eventType);
      m_events.push_back(thisEvent);
   }
   delete events;
}

void diffuseResponses::computeEventResponses() {
   getDiffuseSources();
   std::vector<Event>::iterator it = m_events.begin();
   for (int i = 0; it != m_events.end(); ++it, i++) {
      if (Likelihood::print_output() && (i % (m_events.size()/20)) == 0) {
         std::cerr << ".";
      }
      it->computeResponse(m_srcs, m_srRadius);
   }
   if (Likelihood::print_output()) {
      std::cerr << "!" << std::endl;
   }
}

void diffuseResponses::writeEventResponses() {
   if (m_srcNames.size() > 0) {
      std::string eventFile = m_pars["evfile"];
      facilities::Util::expandEnvVar(&eventFile);
      st_facilities::Util::file_ok(eventFile);
      tip::Table * events 
         = tip::IFileSvc::instance().editTable(eventFile, m_pars["evtable"]);
      if (static_cast<unsigned int>(events->getNumRecords()) 
          != m_events.size()) {
         throw("LogLike::writeEventResponses:\nNumber of records in " 
               + eventFile + " does not match number of events.");
      }
      for (unsigned int i = 0; i < m_srcNames.size(); i++) {
         try {
            if (ResponseFunctions::useEdisp()) {
// Add a 3 dim vector containing the Gaussian parameters describing
// the energy response.
               events->appendField(m_srcNames[i], "3D");
            } else {
// Infinite energy response, so just add the single value.
               events->appendField(m_srcNames[i], "1D");
            }
         } catch (tip::TipException &eObj) {
            if (Likelihood::print_output()) {
               std::cout << eObj.what() << "\n"
                         << "Using existing column." << std::endl;
            }
         }
      }
      tip::Table::Iterator it = events->begin();
      tip::Table::Record & row = *it;
      for (int j = 0 ; it != events->end(); j++, ++it) {
         std::vector<std::string>::iterator name = m_srcNames.begin();
         for ( ; name != m_srcNames.end(); ++name) {
            if (ResponseFunctions::useEdisp()) {
               tip::Table::Vector<double> respParams = row[*name];
               setGaussianParams(m_events[j], *name, respParams);
            } else {
// Assume infinite energy resolution.
               row[*name].set(m_events[j].diffuseResponse(1., *name));
            }
         }
      }
      delete events;
   }
}

void diffuseResponses::setGaussianParams(const Event & event,
                                         const std::string & name,
                                         tip::Table::Vector<double> & params) {
   double norm, mean, sigma;
   event.computeGaussianParams(name, norm, mean, sigma);
   params[0] = norm;
   params[1] = mean;
   params[2] = sigma;
}

void diffuseResponses::getDiffuseSources() {
   m_srcs.clear();
   m_srcNames.clear();
   std::vector<std::string> srcNames;
   m_srcModel->getSrcNames(srcNames);
   for (unsigned int i = 0; i < srcNames.size(); i++) {
      Source * my_src = m_srcModel->getSource(srcNames[i]);
      if (my_src->getType() == std::string("Diffuse")) {
         m_srcs.push_back(dynamic_cast<DiffuseSource *>(my_src));
         m_srcNames.push_back(srcNames[i]);
      }
   }
}

