/**
 * @file gtexpcube2.cxx
 * @brief Application for creating binned exposure maps.
 * @author J. Chiang
w *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtexpcube2/gtexpcube2.cxx,v 1.45 2010/06/02 22:04:16 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/Observation.h"

using namespace Likelihood;

class ExpCube : public st_app::StApp {
public:
   ExpCube();
   virtual ~ExpCube() throw() {
      try {
         delete m_helper;
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;
   double m_srRadius;
   void promptForParameters();
   static std::string s_cvs_id;
};

st_app::StAppFactory<ExpCube> myAppFactory("gtexpcube2");

std::string ExpCube::s_cvs_id("$Name: $");

ExpCube::ExpCube() : st_app::StApp(), m_helper(0), 
                     m_pars(st_app::StApp::getParGroup("gtexpcube2")) {
   setVersion(s_cvs_id);
}

void ExpCube::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void ExpCube::run() {
   promptForParameters();
   m_helper = new AppHelpers(&m_pars, "BINNED");
   CountsMap * cmap(0);
   if (m_pars["cmfile"] != "none") {
      cmap = new CountsMap(m_pars["cmfile"]);
   }
   BinnedExposure bexpmap(cmap, m_helper->observation());
   bexpmap.writeOutput(m_pars["outfile"]);

   delete cmap;
}

void ExpCube::promptForParameters() {
   m_pars.Prompt();
   m_pars.Save();
}
