/**
 * @file diffuseResponses.cxx
 * @brief Adds diffuse response information for extragalactic and Galactic
 * diffuse emission.  Assumes infinite energy resolution.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/likelihood.cxx,v 1.9 2004/06/02 05:27:25 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Util.h"

using namespace Likelihood;

/**
 * @class diffuseResponses
 * @brief FTOOL to add diffuse response information to an FT1 file for
 * extragalactic and Galactic diffuse emission.
 *
 * @author J. Chiang
 *
 * $Header$
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
   LogLike * m_logLike;
   st_app::AppParGroup & m_pars;

   void buildSourceModel();
   void readEventData();
   void computeResponses();
};

st_app::StAppFactory<diffuseResponses> myAppFactory;

diffuseResponses::diffuseResponses() 
   : st_app::StApp(), m_helper(0), m_logLike(0), 
     m_pars(st_app::StApp::getParGroup("diffuseResponses")) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      m_helper = new AppHelpers(m_pars);
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
   buildSourceModel();
   readEventData();
   computeReponses();
}

void diffuseResponses::buildSourceModel() {
   std::string sourceModel = m_pars["Source_model_file"];
   if (m_logLike->getNumSrcs() == 0) {
// Read in the Source model for the first time.
      Util::file_ok(sourceModel);
      m_logLike->readXml(sourceModel, m_helper->funcFactory());
      m_logLike->computeEventResponses();
   } else {
// Re-read the Source model from the xml file, allowing only for 
// Parameter adjustments.
      Util::file_ok(sourceModel);
      m_logLike->reReadXml(sourceModel);
   }
}

void diffuseResponses::readEventData() {
   std::string eventFile = m_pars["event_file"];
   long eventFileHdu = m_pars["event_file_hdu"];
   Util::file_ok(eventFile);
   m_logLike->getEvents(*evIt, eventFileHdu);
}

void diffuseResponses::computeResponses() {
   m_logLike->computeEventResponses();
   m_logLike->writeEventResponses(m_pars["event_file"]);
}
