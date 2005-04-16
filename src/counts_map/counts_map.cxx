#include <cstdlib>

#include <iostream>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

#include "Likelihood/CountsMap.h"

using namespace Likelihood;

class gtcntsmap : public st_app::StApp {
public:
   gtcntsmap();
   virtual ~gtcntsmap() throw() {
      try {
         delete m_helper;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
         std::cerr << "gtcntsmap::~gtcntsmap: "
                   << "unknown exception encountered."
                   << std::endl;
      }
   }
   virtual void run();
private:
   st_app::AppParGroup & m_pars;
   void logArray(double xmin, double xmax, unsigned int nx,
                 std::vector<double> & xx) const;
};

st_app::StAppFactory<gtcntsmap> myAppFactory;

gtcntsmap::gtcntsmap() : public st_app::StApp(), 
   m_pars(st_app::StApp::getParGroup("gtcntsmap")) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      m_helper = new AppHelpers(&m_pars);
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
   std::string event_file = m_pars["event_file"];
   std::vector<std::string> eventFiles;
   st_facilities::Util::resolve_fits_files(event_file, eventFiles);

   std::string sc_file = m_pars["spacecraft_file"];
   std::vector<std::string> scDataFiles;
   st_facilities::Util::resolve_fits_files(sc_file, scDataFiles);

   double emin = m_pars["emin"];
   double emax = m_pars["emax"];
   unsigned long nenergies = m_pars["nenergies"];

   std::vector<double> energies;
   logArray(emin, emax, nenergies, energies);

   double ra = m_pars["ra"];
   double dec = m_pars["dec"];
   unsigned long nra = m_pars["nra"];
   unsigned long ndec = m_pars["ndec"];
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
   std::string output_file = m_pars["output_file_name"];
   cmap.writeOutput("gtcntsmap", output_file);
}

void gtcntsmap::logArray(double xmin, double xmax, unsigned int nx,
                         std::vector<double> & xx) const {
   double xstep = log(xmax/xmin)/(nx - 1.);
   for (unsigned int i = 0; i < nx; i++) {
      xx.push_back(xmin*exp(i*xstep));
   }
}
