/**
 * @file gtexpcube2.cxx
 * @brief Application for creating binned exposure maps.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtexpcube2/gtexpcube2.cxx,v 1.8 2010/11/30 07:51:13 jchiang Exp $
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

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Header.h"

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
   void set_phi_status();
   void generateEnergies(std::vector<double> & energies) const;
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

   std::string cmap_file = m_pars["cmap"];
   std::string ltcube_file = m_pars["infile"];

   bool useEbounds(true);
   std::string bincalc = m_pars["bincalc"];
   if (bincalc == "CENTER") {
      useEbounds = false;
   }
   set_phi_status();

   if (cmap_file != "none") {
      m_helper = new AppHelpers(&m_pars, "BINNED");
   } else {
      m_helper = new AppHelpers(&m_pars, "NONE");
   }
   m_helper->checkOutputFile();
   ExposureCube & ltcube = 
      const_cast<ExposureCube &>(m_helper->observation().expCube());
   ltcube.readExposureCube(ltcube_file);

   if (cmap_file != "none") {
// Create map to match counts map.
      m_helper->checkTimeCuts(cmap_file, "", ltcube_file, "Exposure");
      CountsMap cmap(cmap_file);
      BinnedExposure bexpmap(cmap, m_helper->observation(), useEbounds);
      bexpmap.writeOutput(m_pars["outfile"]);
      return;
   }

// Create map for user-defined geometry.
   std::vector<double> energies;
   generateEnergies(energies);
   if (!useEbounds) {
      for (size_t k(0); k < energies.size() - 1; k++) {
         energies[k] = std::sqrt(energies[k]*energies[k+1]);
      }
      energies.pop_back();
   }
   BinnedExposure bexpmap(energies, m_helper->observation(), &m_pars);
   bexpmap.writeOutput(m_pars["outfile"]);
}

void ExpCube::generateEnergies(std::vector<double> & energies) const {
   double emin = m_pars["emin"];
   double emax = m_pars["emax"];
   int enumbins = m_pars["enumbins"];
   double estep = std::log(emax/emin)/enumbins;
   energies.clear();
   for (size_t k(0); k < enumbins + 1; k++) {
      energies.push_back(emin*std::exp(estep*k));
   }
}

void ExpCube::promptForParameters() {
   m_pars.Prompt("infile");
   m_pars.Prompt("cmap");
   m_pars.Prompt("outfile");
   m_pars.Prompt("irfs");
   std::string cmap = m_pars["cmap"];
   if (cmap == "none") {
      m_pars.Prompt("nxpix");
      m_pars.Prompt("nypix");
      m_pars.Prompt("binsz");
      m_pars.Prompt("coordsys");
      m_pars.Prompt("xref");
      m_pars.Prompt("yref");
      m_pars.Prompt("axisrot");
      m_pars.Prompt("proj");
      m_pars.Prompt("emin");
      m_pars.Prompt("emax");
      m_pars.Prompt("enumbins");
   }
   m_pars.Save();
}

void ExpCube::set_phi_status() {
   // If indicated, turn off phi-dependence for all IRFs.
   bool ignorephi = m_pars["ignorephi"];
   if (ignorephi) {
      const Observation & observation(m_helper->observation());
      std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
         = observation.respFuncs().begin();
      for ( ; respIt != observation.respFuncs().end(); ++respIt) {
         respIt->second->aeff()->setPhiDependence(false);
      }
   }
}
