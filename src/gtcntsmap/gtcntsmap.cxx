/**
 * @file gtcntsmap.cxx
 * @brief Creates counts maps for use by binned likelihood.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtcntsmap/gtcntsmap.cxx,v 1.16 2006/04/18 05:43:45 jchiang Exp $
 */

#include <cstdlib>

#include <iostream>
#include <map>
#include <stdexcept>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/RoiCuts.h"

using namespace Likelihood;

class gtcntsmap : public st_app::StApp {
public:
   gtcntsmap();
   virtual ~gtcntsmap() throw() {}
   virtual void run();
   virtual void banner() const;
private:
   st_app::AppParGroup & m_pars;
   dataSubselector::Cuts * m_cuts;

   void promptForParameters();
   void writeDssCuts() const;
   void logArray(double xmin, double xmax, unsigned int nx,
                 std::vector<double> & xx) const;
   void checkEnergies(double emin, double emax) const;

   static std::string s_cvs_id;
};

st_app::StAppFactory<gtcntsmap> myAppFactory("gtcntsmap");

gtcntsmap::gtcntsmap() : st_app::StApp(), 
   m_pars(st_app::StApp::getParGroup("gtcntsmap")), m_cuts(0) {
   setVersion(s_cvs_id);
   m_pars.setSwitch("use_lb");
   m_pars.setCase("use_lb", "yes", "glon");
   m_pars.setCase("use_lb", "yes", "glat");
   m_pars.setCase("use_lb", "no", "ra");
   m_pars.setCase("use_lb", "no", "dec");
}

std::string gtcntsmap::s_cvs_id("$Name:  $");

void gtcntsmap::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void gtcntsmap::run() {
   promptForParameters();
   AppHelpers::checkOutputFile(m_pars["clobber"], m_pars["outfile"]);
                               
   std::string event_file = m_pars["evfile"];
   std::string evtable = m_pars["evtable"];
   std::vector<std::string> eventFiles;
   st_facilities::Util::resolve_fits_files(event_file, eventFiles);
   bool compareGtis(false);
   bool relyOnStreams(false);
   bool skipEventClassCuts(true);
   for (unsigned int i = 1; i < eventFiles.size(); i++) {
      AppHelpers::checkCuts(eventFiles[0], evtable,
                            eventFiles[i], evtable,
                            compareGtis, relyOnStreams, 
                            skipEventClassCuts);
   }

   m_cuts = new dataSubselector::Cuts(eventFiles, evtable);

   std::string sc_file = m_pars["scfile"];
   std::string sc_table = m_pars["sctable"];

   std::vector<std::string> scDataFiles;
   st_facilities::Util::resolve_fits_files(sc_file, scDataFiles);

   double emin = m_pars["emin"];
   double emax = m_pars["emax"];
   checkEnergies(emin, emax);
   long nenergies = m_pars["nenergies"];

   std::vector<double> energies;
   logArray(emin, emax, nenergies, energies);

   double ra = m_pars["ra"];
   double dec = m_pars["dec"];
   bool use_lb = m_pars["use_lb"];
   if (use_lb) {
      ra = m_pars["glon"];
      dec = m_pars["glat"];
   }
   long nra = m_pars["nra"];
   long ndec = m_pars["ndec"];
   double pixel_size = m_pars["pixel_size"];
   std::string projection = m_pars["proj"];
   CountsMap cmap(eventFiles[0], evtable, scDataFiles[0], sc_table,
                  ra, dec, projection, nra, ndec, pixel_size, 0, use_lb, 
                  "RA", "DEC", energies);
   
   for (unsigned int i = 0; i < eventFiles.size(); i++) {
      const tip::Table * events 
         = tip::IFileSvc::instance().readTable(eventFiles[i], "events");
      cmap.binInput(events->begin(), events->end());
      delete events;
   }
   std::string output_file = m_pars["outfile"];
   cmap.writeOutput("gtcntsmap", output_file);
   writeDssCuts();
   delete m_cuts;
}

void gtcntsmap::promptForParameters() {
   m_pars.Prompt("evfile");
   m_pars.Prompt("scfile");
   m_pars.Prompt("outfile");
   m_pars.Prompt("emin");
   m_pars.Prompt("emax");
   m_pars.Prompt("nenergies");
   m_pars.Prompt("use_lb");
   bool use_lb = m_pars["use_lb"];
   if (use_lb) {
      m_pars.Prompt("glon");
      m_pars.Prompt("glat");
   } else {
      m_pars.Prompt("ra");
      m_pars.Prompt("dec");
   }
   m_pars.Prompt("nra");
   m_pars.Prompt("ndec");
   m_pars.Prompt("pixel_size");
   m_pars.Prompt("proj");
   m_pars.Save();
}

void gtcntsmap::writeDssCuts() const {
   std::string output_file = m_pars["outfile"];
   std::auto_ptr<tip::Image> 
      image(tip::IFileSvc::instance().editImage(output_file, ""));
   m_cuts->writeDssKeywords(image->getHeader());
   m_cuts->writeGtiExtension(output_file);
}

void gtcntsmap::logArray(double xmin, double xmax, unsigned int nx,
                         std::vector<double> & xx) const {
   double xstep = log(xmax/xmin)/(nx - 1.);
   for (unsigned int i = 0; i < nx; i++) {
      xx.push_back(xmin*exp(i*xstep));
   }
}

void gtcntsmap::checkEnergies(double emin, double emax) const {
   std::map<std::string, double> pars1, pars2;
   pars1["ENERGY"] = emin;
   pars2["ENERGY"] = emax;
   if (!m_cuts->accept(pars1) || !m_cuts->accept(pars2)) {
      std::ostringstream message;
      message << "The requested energies, " << emin << ", " << emax
              << " do not lie within the event file DSS selections:\n\n";
      m_cuts->writeCuts(message);
      throw std::runtime_error(message.str());
   }
}
