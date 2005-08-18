/**
 * @file gtcntsmap.cxx
 * @brief Creates counts maps for use by binned likelihood.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtcntsmap/gtcntsmap.cxx,v 1.8 2005/05/17 13:44:14 jchiang Exp $
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

#include "Verbosity.h"

using namespace Likelihood;

class gtcntsmap : public st_app::StApp {
public:
   gtcntsmap();
   virtual ~gtcntsmap() throw() {}
   virtual void run();
   virtual void banner() const {}
private:
   st_app::AppParGroup & m_pars;
   dataSubselector::Cuts * m_cuts;

   void writeDssCuts() const;
   void logArray(double xmin, double xmax, unsigned int nx,
                 std::vector<double> & xx) const;
   void checkEnergies(double emin, double emax) const;
};

st_app::StAppFactory<gtcntsmap> myAppFactory("gtcntsmap");

gtcntsmap::gtcntsmap() : st_app::StApp(), 
   m_pars(st_app::StApp::getParGroup("gtcntsmap")), m_cuts(0) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      Likelihood::Verbosity::instance(m_pars["chatter"]);
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
      std::exit(1);
   } catch (...) {
      std::cerr << "Caught unknown exception in gtcntsmap constructor." 
                << std::endl;
      std::exit(1);
   }
}

void gtcntsmap::run() {
   AppHelpers::checkOutputFile(m_pars["clobber"], m_pars["outfile"]);
                               
   std::string event_file = m_pars["evfile"];
   std::vector<std::string> eventFiles;
   st_facilities::Util::resolve_fits_files(event_file, eventFiles);
   bool compareGtis(false);
   for (unsigned int i = 1; i < eventFiles.size(); i++) {
      AppHelpers::checkCuts(eventFiles[0], "EVENTS",
                            eventFiles[i], "EVENTS",
                            compareGtis);
   }

   m_cuts = new dataSubselector::Cuts(eventFiles[0]);

   std::string sc_file = m_pars["scfile"];
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
   long nra = m_pars["nra"];
   long ndec = m_pars["ndec"];
   double pixel_size = m_pars["pixel_size"];
   bool use_lb = m_pars["use_lb"];
   CountsMap cmap(eventFiles[0], scDataFiles[0], ra, dec, "CAR",
                  nra, ndec, pixel_size, 0, use_lb, "RA", "DEC", 
                  energies);
   
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
