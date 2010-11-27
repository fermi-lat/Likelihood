/**
 * @file gtexpcube2.cxx
 * @brief Application for creating binned exposure maps.
 * @author J. Chiang
w *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtexpcube2/gtexpcube2.cxx,v 1.1 2010/11/24 05:11:27 jchiang Exp $
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
   void generateEnergies(bool useEbounds, std::vector<double> & energies) const;
   static std::string s_cvs_id;
};

st_app::StAppFactory<ExpCube> myAppFactory("gtexpcube2");

std::string ExpCube::s_cvs_id("$Name:  $");

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
   m_helper->checkOutputFile();
   std::string ltcube_file = m_pars["infile"];
   ExposureCube & ltcube = 
      const_cast<ExposureCube &>(m_helper->observation().expCube());
   ltcube.readExposureCube(ltcube_file);

   bool useEbounds(true);
   if (m_pars["bincalc"] == "CENTER") {
      useEbounds = false;
   }
   if (m_pars["cmap"] != "none") {
      CountsMap cmap(m_pars["cmap"]);
      BinnedExposure bexpmap(cmap, m_helper->observation(), useEbounds);
      bexpmap.writeOutput(m_pars["outfile"]);
      return;
   }
   std::vector<double> energies;
   generateEnergies(useEbounds, energies);
   BinnedExposure bexpmap(energies, m_pars["proj"], m_pars["coordsys"],
                          m_helper->observation());
   bexpmap.writeOutput(m_pars["outfile"]);
}

void ExpCube::generateEnergies(bool useEbounds,
                               std::vector<double> & energies) const {
   double emin = m_pars["emin"];
   double emax = m_pars["emax"];
   size_t enumbins = m_pars["enumbins"];
   double estep = std::log(emax/emin)/(enumbins - 1);
   energies.clear();
   for (size_t k(0); k < enumbins; k++) {
      energies.push_back(emin*std::exp(estep*k));
   }
   if (!useEbounds) {
      for (size_t k(0); k < enumbins - 1; k++) {
         energies[k] = std::sqrt(energies[k]*energies[k+1]);
      }
      energies.pop_back();
   }
}

void ExpCube::promptForParameters() {
   m_pars.Prompt("infile");
   m_pars.Prompt("cmap");
   m_pars.Prompt("outfile");
   m_pars.Prompt("irfs");
   if (m_pars["cmap"] == "none") {
      m_pars.Prompt("emin");
      m_pars.Prompt("emax");
      m_pars.Prompt("enumbins");
      m_pars.Prompt("coordsys");
      m_pars.Prompt("proj");
   }
   m_pars.Save();
}
