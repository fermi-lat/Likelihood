/**
 * @file gtexpcube2.cxx
 * @brief Application for creating binned exposure maps.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtexpcube2/gtexpcube2.cxx,v 1.21 2015/01/16 05:36:28 jchiang Exp $
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

#include "dataSubselector/Cuts.h"
#include "dataSubselector/Gti.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedExposure.h"
#include "Likelihood/BinnedHealpixExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/Observation.h"

#include "evtbin/HealpixMap.h"

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
   void copyGtis() const;
   void copyHeaderKeywords() const;
   void copyDssKeywords(tip::Header & header) const;
   bool have_user_sky_geom() const;
   bool have_user_energies() const;
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

   if (cmap_file != "none") {
      m_helper = new AppHelpers(&m_pars, "BINNED");
   } else {
      m_helper = new AppHelpers(&m_pars, "NONE");
   }
   set_phi_status();
   m_helper->checkOutputFile();
   ExposureCube & ltcube = 
      const_cast<ExposureCube &>(m_helper->observation().expCube());
   ltcube.readExposureCube(ltcube_file);

   if (cmap_file != "none") {
      m_helper->checkTimeCuts(cmap_file, "", ltcube_file, "Exposure");

      // Handle HEALPIX case.
      const tip::Image * image(tip::IFileSvc::instance().readImage(cmap_file,
                                                                   ""));
      int image_size(image->getImageDimensions().size());
      delete image;
      if (image_size == 0) {
         const tip::Table * table = 
            tip::IFileSvc::instance().readTable(cmap_file, "SKYMAP");
         std::string hpx("");
         const tip::Header & header(table->getHeader());
         header["PIXTYPE"].get(hpx);
         delete table;
         if (hpx == "HEALPIX") {
            evtbin::HealpixMap cmap(cmap_file);
            BinnedHealpixExposure bexpmap(cmap, m_helper->observation(), 
                                          useEbounds, &m_pars);
            bexpmap.writeOutput(m_pars["outfile"]);
         } else {
            throw std::runtime_error("Unrecognized extension for the "
                                     "input count map");
         }
         copyHeaderKeywords();
         copyGtis();
         return;
      }
   }
      
   // Conventional map geometries, with possible overrides by user.
   bool user_energies(have_user_energies());
   bool user_sky_geom(have_user_sky_geom());
   if (user_energies || user_sky_geom) {
      // Create a local copy of the m_pars object.
      st_app::AppParGroup & pars(m_pars);
      // Create map for user-defined geometry.
      std::vector<double> energies;
      if (user_energies) {
         generateEnergies(energies);
         if (!useEbounds) {
            for (size_t k(0); k < energies.size() - 1; k++) {
               energies[k] = std::sqrt(energies[k]*energies[k+1]);
            }
            energies.pop_back();
         }
      } else {
         // Harvest energies from the counts map.
         const tip::Table * 
            ebounds(tip::IFileSvc::instance().readTable(cmap_file, "EBOUNDS"));
         tip::Table::ConstIterator it(ebounds->begin());
         tip::ConstTableRecord & row(*it);
         double emin, emax;
         for ( ; it != ebounds->end(); ++it) {
            row["E_MIN"].get(emin);
            row["E_MAX"].get(emax);
            energies.push_back(emin/1e3);
         }
         energies.push_back(emax/1e3);
         delete ebounds;
      }
      if (!user_sky_geom) {
         // Harvest sky geometry from counts map.
         CountsMap cmap(cmap_file);
         pars["nxpix"] = cmap.naxis1();
         pars["nypix"] = cmap.naxis2();
         pars["binsz"] = std::fabs(cmap.cdelt1());
         pars["xref"] = cmap.crval1();
         pars["yref"] = cmap.crval2();
         pars["axisrot"] = cmap.crota2();
         pars["proj"] = cmap.proj_name();
         if (cmap.isGalactic()) {
            pars["coordsys"] = "GAL";
         } else {
            pars["coordsys"] = "CEL";
         }
      }
      BinnedExposure bexpmap(energies, m_helper->observation(), &pars);
      bexpmap.writeOutput(m_pars["outfile"]);
   } else {
      CountsMap cmap(cmap_file);
      BinnedExposure bexpmap(cmap, m_helper->observation(), 
                             useEbounds, &m_pars);
      bexpmap.writeOutput(m_pars["outfile"]);
   }
   copyHeaderKeywords();
   copyGtis();
}

bool ExpCube::have_user_sky_geom() const {
   // Check for INDEFs in m_pars.
   try {
      unsigned int nxpix = m_pars["nxpix"];
      unsigned int nypix = m_pars["nypix"];
      double binsz = m_pars["binsz"];
      double xref = m_pars["xref"];
      double yref = m_pars["yref"];
   } catch (hoops::Hexception &) {
      // At least one of the parameters needed to specify the map
      // geometry on the sky has an INDEF value, so use input cmap sky
      // geometry.
      return false;
   }
   return true;
}

bool ExpCube::have_user_energies() const {
   // Check for INDEFs in m_pars.
   std::string ebinalg = m_pars["ebinalg"];
   if (ebinalg == "FILE") {
      std::string ebinfile = m_pars["ebinfile"];
      if (ebinfile != "NONE") {
         return true;
      }
      return false;
   }
   try {
      double emin = m_pars["emin"];
      double emax = m_pars["emax"];
      unsigned int enumbins = m_pars["enumbins"];
   } catch (hoops::Hexception &) {
      // At least one of the energy bounds parameters has a value of
      // INDEF, so use cmap energy bounds.
      return false;
   }
   return true;
}

void ExpCube::generateEnergies(std::vector<double> & energies) const {
   energies.clear();
   std::string ebinalg = m_pars["ebinalg"];
   if (ebinalg == "FILE") {
      std::string ebinfile = m_pars["ebinfile"];
      const tip::Table * energybins = 
         tip::IFileSvc::instance().readTable(ebinfile, "ENERGYBINS");
      tip::Table::ConstIterator it = energybins->begin();
      tip::ConstTableRecord & row = *it;
      double energy;
      double emax;
      for ( ; it != energybins->end(); ++it) {
         row["E_MIN"].get(energy);
         // Note that energies in gtbindef output are in units of keV.
         energies.push_back(energy/1e3);
         row["E_MAX"].get(emax);
      }
      energies.push_back(emax/1e3);
      delete energybins;
   } else {
      double emin = m_pars["emin"];
      double emax = m_pars["emax"];
      int enumbins = m_pars["enumbins"];
      double estep = std::log(emax/emin)/enumbins;
      for (size_t k(0); k < enumbins + 1; k++) {
         energies.push_back(emin*std::exp(estep*k));
      }
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
      std::string ebinalg = m_pars["ebinalg"];
      if (ebinalg == "FILE") {
         m_pars.Prompt("ebinfile");
      } else {
         m_pars.Prompt("emin");
         m_pars.Prompt("emax");
         m_pars.Prompt("enumbins");
      }
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

void ExpCube::copyGtis() const {
   std::string infile = m_pars["infile"];
   dataSubselector::Gti gti(infile);
   std::string outfile = m_pars["outfile"];
   gti.writeExtension(outfile);
}

void ExpCube::copyHeaderKeywords() const {
   std::string infile = m_pars["infile"];
   std::string inext("EXPOSURE");
   const tip::Table * intab = 
     tip::IFileSvc::instance().readTable(infile, inext);
   const tip::Header & inheader(intab->getHeader());
   
   std::string outfile = m_pars["outfile"];
   std::string outext("PRIMARY");
   tip::Image * outimg = tip::IFileSvc::instance().editImage(outfile, outext);
   tip::Header & outheader(outimg->getHeader());

   // Unfortunately TIP does not provide access to header comments, so
   // we cannot copy them using it. Could use cfitsio directly but
   // prefer TIP as it is standard here.

#define COPYKEYWORD(type, name)			\
   try {					\
     type x;					\
     inheader[name].get(x);			\
     outheader[name].set(x);			\
   } catch(...) {				\
   }
   
   COPYKEYWORD(std::string, "DATE-OBS");
   COPYKEYWORD(std::string, "DATE-END");
   COPYKEYWORD(double,      "TSTART");
   COPYKEYWORD(double,      "TSTOP");
   COPYKEYWORD(double,      "MJDREFI");
   COPYKEYWORD(double,      "MJDREFF");
   COPYKEYWORD(std::string, "TIMEUNIT");
   COPYKEYWORD(double,      "TIMEZERO");
   COPYKEYWORD(std::string, "TIMESYS");
   COPYKEYWORD(std::string, "TIMEREF");
   COPYKEYWORD(bool,        "CLOCKAPP");
   COPYKEYWORD(bool,        "GPS_OUT");

#undef COPYKEYWORD

   copyDssKeywords(outheader);

   delete outimg;
   delete intab;
}

void ExpCube::copyDssKeywords(tip::Header & header) const {
   std::string cmap_file = m_pars["cmap"];
   dataSubselector::Cuts * irfs_cuts(0);
   if (cmap_file == "none") {
      // No DSS keywords to copy from counts map file, so
      // use livetime cube, in order to get GTIs.
      std::string ltcube = m_pars["infile"];
      irfs_cuts = new dataSubselector::Cuts(ltcube, "EXPOSURE", false,
                                            false, false);
   } else {
      // Copy DSS keywords from input cmap file, ensuring that the
      // irfs version is set.
      irfs_cuts = new dataSubselector::Cuts(cmap_file, "PRIMARY", false,
                                            false, false);
   }
   std::string irfs = m_pars["irfs"];
   if (irfs != "CALDB") {
      // Update the BitMaskCuts to the user's selection via the irfs
      // and evtype options.
      m_helper->setBitMaskCuts(*irfs_cuts);
   }
   irfs_cuts->addVersionCut("IRF_VERSION", m_helper->irfsName());
   irfs_cuts->writeDssKeywords(header);
   delete irfs_cuts;
}
