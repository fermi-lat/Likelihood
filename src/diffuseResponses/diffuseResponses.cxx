/**
 * @file diffuseResponses.cxx
 * @brief Adds diffuse response information for extragalactic and Galactic
 * diffuse emission.  Assumes infinite energy resolution.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/diffuseResponses/diffuseResponses.cxx,v 1.1 2004/06/05 00:28:00 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SourceModel.h"
#include "Likelihood/Util.h"

using namespace Likelihood;

/**
 * @class diffuseResponses
 * @brief FTOOL to add diffuse response information to an FT1 file for
 * extragalactic and Galactic diffuse emission.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/diffuseResponses/diffuseResponses.cxx,v 1.1 2004/06/05 00:28:00 jchiang Exp $
 */

class diffuseResponses : public st_app::StApp {

public:

   diffuseResponses();

   virtual ~diffuseResponses() throw() {
      try {
      } catch (std::exception & eObj) {
         std::cout << eObj.what() << std::endl;
      } catch (...) {
      }
   }

   virtual void run();

private:

   AppHelpers * m_helper;  //blech.
   SourceModel * m_srcModel;
   double m_srRadius;
   st_app::AppParGroup & m_pars;

   std::vector<Event> m_events;
   std::vector<std::string> m_srcNames;

   void setRoi();
   void buildSourceModel();
   void readEventData();
   void computeEventResponses();
   void writeEventResponses();
   void getDiffuseSources(std::vector<DiffuseSource *> &srcs);

};

st_app::StAppFactory<diffuseResponses> myAppFactory;

diffuseResponses::diffuseResponses() 
   : st_app::StApp(), m_helper(0), m_srcModel(0), m_srRadius(30.),
     m_pars(st_app::StApp::getParGroup("diffuseResponses")) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      m_helper = new AppHelpers(m_pars);
      m_srcModel = new SourceModel(true);
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
      std::exit(1);
   } catch (...) {
      std::cerr << "Caught unknown exception in diffuseResponses constructor." 
                << std::endl;
      std::exit(1);
   }
}

void diffuseResponses::run() {
   setRoi();
   buildSourceModel();
   readEventData();
   computeEventResponses();
   writeEventResponses();
}

void diffuseResponses::setRoi() {
   RoiCuts * roiCuts = RoiCuts::instance();
   roiCuts->setCuts();
}

void diffuseResponses::buildSourceModel() {
   std::string sourceModel = m_pars["Source_model_file"];
   Util::file_ok(sourceModel);
   m_srcModel->readXml(sourceModel, m_helper->funcFactory(), false);
}

void diffuseResponses::readEventData() {
   std::string eventFile = m_pars["event_file"];
   facilities::Util::expandEnvVar(&eventFile);
   Util::file_ok(eventFile);
   tip::Table * events 
      = tip::IFileSvc::instance().editTable(eventFile, "events");

   double ra;
   double dec;
   double energy;
   double time;
   double raSCZ;
   double decSCZ;
   double zenAngle;
   int convLayer;
   int eventType;

   ScData * scData = ScData::instance();

   tip::Table::Iterator it = events->begin();
   tip::Table::Record & event = *it;
   for ( ; it != events->end(); ++it) {
      event["ra"].get(ra);
      event["dec"].get(dec);
      event["energy"].get(energy);
      event["time"].get(time);
      raSCZ = scData->zAxis(time).ra();
      decSCZ = scData->zAxis(time).dec();
      event["zenith_angle"].get(zenAngle);
      event["conversion_layer"].get(convLayer);
      if (convLayer < 12) { // Front
         eventType = 0;
      } else {
         eventType = 1;
      }
      Event thisEvent(ra, dec, energy, time, raSCZ, decSCZ, 
                      cos(zenAngle*M_PI/180.), eventType);
      m_events.push_back(thisEvent);
   }
   delete events;
}

void diffuseResponses::computeEventResponses() {
   std::vector<DiffuseSource *> srcs;
   getDiffuseSources(srcs);
   std::vector<Event>::iterator it = m_events.begin();
   for (int i = 0; it != m_events.end(); ++it, i++) {
      if ((i % (m_events.size()/20)) == 0) std::cerr << ".";
      it->computeResponse(srcs, m_srRadius);
   }
   std::cerr << "!" << std::endl;
}

void diffuseResponses::writeEventResponses() {
   if (m_srcNames.size() > 0) {
      std::string eventFile = m_pars["event_file"];
      facilities::Util::expandEnvVar(&eventFile);
      Util::file_ok(eventFile);
      tip::Table * events 
         = tip::IFileSvc::instance().editTable(eventFile, "events");
      if (static_cast<unsigned int>(events->getNumRecords()) 
          != m_events.size()) {
         throw("LogLike::writeEventResponses:\nNumber of records in " 
               + eventFile + " does not match number of events.");
      }
      for (unsigned int i = 0; i < m_srcNames.size(); i++) {
         try {
            events->appendField(m_srcNames[i], "1D");
         } catch (tip::TipException &eObj) {
            std::cout << eObj.what() << "\n"
                      << "Using existing column." << std::endl;
         }
      }
      tip::Table::Iterator it = events->begin();
      tip::Table::Record & row = *it;
      for (int j = 0 ; it != events->end(); j++, ++it) {
         std::vector<std::string>::iterator name = m_srcNames.begin();
         for ( ; name != m_srcNames.end(); ++name) {
// For now, assume infinite energy resolution.
            row[*name].set(m_events[j].diffuseResponse(1., *name));
         }
      }
      delete events;
   }
}

void diffuseResponses::getDiffuseSources(std::vector<DiffuseSource *> &srcs) {
   srcs.clear();
   m_srcModel->getSrcNames(m_srcNames);
   for (unsigned int i = 0; i < m_srcNames.size(); i++) {
      Source * my_src = m_srcModel->getSource(m_srcNames[i]);
      if (my_src->getType() == std::string("Diffuse")) {
         srcs.push_back(dynamic_cast<DiffuseSource *>(my_src));
      }
   }
}

