/**
 * @file gtexpcube2.cxx
 * @brief Application for creating binned exposure maps.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtexpcube2/gtexpcube2.cxx,v 1.22 2015/01/16 17:27:29 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <stdexcept>

#include "astro/ProjBase.h"

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
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/Observation.h"

#include "evtbin/HealpixMap.h"

using namespace Likelihood;

class ExpCube : public st_app::StApp {

private:
   static bool check_cmap_null(const std::string& cmap_name);
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

   astro::ProjBase::Method m_projMethod;
   bool m_useEbounds;

   void run_wcs(const std::string* cmap_file);    
   void run_hpx(const std::string* cmap_file);
 
   void promptForParameters();
   void set_phi_status();
   void generateEnergies(std::vector<double> & energies,
			 const std::string* cmap_file) const;

   void copyGtis() const;
   void copyHeaderKeywords() const;
   void copyDssKeywords(tip::Header & header) const;
   bool have_user_sky_geom_wcs() const;
   bool have_user_sky_geom_hpx() const;
   bool have_user_energies() const;
   static std::string s_cvs_id;
};

st_app::StAppFactory<ExpCube> myAppFactory("gtexpcube2");

std::string ExpCube::s_cvs_id("$Name:  $");

bool ExpCube::check_cmap_null(const std::string& cmap_name) {
   if ( cmap_name == "none" ||
	cmap_name == "hpx" ) {
      return true;
   }
   return false;
}

ExpCube::ExpCube() : st_app::StApp(), m_helper(0), 
		     m_projMethod(astro::ProjBase::UNKNOWN),
		     m_useEbounds(true),
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
     m_useEbounds = false;
   }

   bool cmap_is_null = check_cmap_null(cmap_file);
   
   if ( !cmap_is_null ) {
      m_helper = new AppHelpers(&m_pars, "NONE");
   } else {
      m_helper = new AppHelpers(&m_pars, "BINNED");
   }
   set_phi_status();
   m_helper->checkOutputFile();
   ExposureCube & ltcube = 
      const_cast<ExposureCube &>(m_helper->observation().expCube());
   ltcube.readExposureCube(ltcube_file);

   if ( !cmap_is_null ) {
      m_helper->checkTimeCuts(cmap_file, "", ltcube_file, "Exposure");
      m_projMethod = AppHelpers::checkProjectionMethod(cmap_file, "");
   } else {
     if ( cmap_file == "none" ) {
       m_projMethod = astro::ProjBase::WCS;
     } else if ( cmap_file == "hpx" ) {
       m_projMethod = astro::ProjBase::HEALPIX;
     }
   }

   switch ( m_projMethod ) {
   case astro::ProjBase::WCS:
     run_wcs(cmap_is_null ? 0 : &cmap_file);
     break;
   case astro::ProjBase::HEALPIX:
     run_hpx(cmap_is_null ? 0 : &cmap_file);
     break;
   default:
     throw std::runtime_error("Unknown projection type");
   }

   copyHeaderKeywords();
   copyGtis();
}


void ExpCube::run_wcs(const std::string* cmap_file) {

   // Conventional map geometries, with possible overrides by user.
   bool user_energies(have_user_energies());
   bool user_sky_geom(have_user_sky_geom_wcs());

   if (user_energies || user_sky_geom ) {
      // Create a local copy of the m_pars object.
      st_app::AppParGroup & pars(m_pars);
      // Create map for user-defined geometry.
      std::vector<double> energies;
      generateEnergies( energies, user_energies ? 0 : cmap_file);
      if (!m_useEbounds) {
	for (size_t k(0); k < energies.size() - 1; k++) {
	  energies[k] = std::sqrt(energies[k]*energies[k+1]);
	}
	energies.pop_back();
      } 
      if (!user_sky_geom) {
	if ( cmap_file == 0 ) {
	  throw std::runtime_error("Geometry not fully defined and no input counts map file");
	}
	// Harvest sky geometry from counts map.
	CountsMap cmap(*cmap_file);
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
      if ( cmap_file == 0 ) {
	 throw std::runtime_error("Geometry not fully defined and no input counts map file");
      }
      CountsMap cmap(*cmap_file);
      BinnedExposure bexpmap(cmap, m_helper->observation(), 
                             m_useEbounds, &m_pars);
      bexpmap.writeOutput(m_pars["outfile"]);
   }
}


void ExpCube::run_hpx(const std::string* cmap_file) {

   // Conventional map geometries, with possible overrides by user.
   bool user_energies(have_user_energies());
   bool user_sky_geom(have_user_sky_geom_hpx());

   if (user_energies || user_sky_geom ) {
      // Create a local copy of the m_pars object.
      st_app::AppParGroup & pars(m_pars);
      // Create map for user-defined geometry.
      std::vector<double> energies;
      generateEnergies( energies, user_energies ? 0 : cmap_file);
      if (!m_useEbounds) {
	for (size_t k(0); k < energies.size() - 1; k++) {
	  energies[k] = std::sqrt(energies[k]*energies[k+1]);
	}
	energies.pop_back();
      } 
      if (!user_sky_geom) {
         // Harvest sky geometry from counts map.
	 if ( cmap_file == 0 ) {
	   throw std::runtime_error("Geometry not fully defined and no input counts map file");
	 }
         CountsMapHealpix cmap(*cmap_file);
         pars["hpx_order"] = cmap.healpixProj()->healpix().Order();
         pars["hpx_ordering_scheme"] = cmap.healpixProj()->healpix().Scheme() == RING ? "RING" : "NESTED";
         if (cmap.isGalactic()) {
            pars["coordsys"] = "GAL";
         } else {
            pars["coordsys"] = "CEL";
         }
      }
      BinnedHealpixExposure bexpmap(energies, m_helper->observation(), &pars);
      bexpmap.writeOutput(m_pars["outfile"]);
   } else {
      if ( cmap_file == 0 ) {
        throw std::runtime_error("Geometry not fully defined and no input counts map file");
      }
      evtbin::HealpixMap cmap(*cmap_file);
      BinnedHealpixExposure bexpmap(cmap, m_helper->observation(), 
				    m_useEbounds, &m_pars);
      bexpmap.writeOutput(m_pars["outfile"]);
   }
}




bool ExpCube::have_user_sky_geom_wcs() const {
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

bool ExpCube::have_user_sky_geom_hpx() const {
   // Check for INDEFs in m_pars.
   try {
      std::string hpx_scheme = m_pars["hpx_ordering_scheme"];
      unsigned int hpx_order = m_pars["hpx_order"];
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

void ExpCube::generateEnergies(std::vector<double> & energies,
			       const std::string* cmap_file) const {
   energies.clear();

   if ( cmap_file != 0 ){
     // Harvest energies from the counts map.
     const tip::Table * 
       ebounds(tip::IFileSvc::instance().readTable(*cmap_file, "EBOUNDS"));
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
   } else {
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
}

void ExpCube::promptForParameters() {
   m_pars.Prompt("infile");
   m_pars.Prompt("cmap");
   m_pars.Prompt("outfile");
   m_pars.Prompt("irfs");
   std::string cmap = m_pars["cmap"];
   bool cmap_is_null = check_cmap_null(cmap);

   if ( cmap_is_null ) {
     if (cmap == "none") {
       m_pars.Prompt("nxpix");
       m_pars.Prompt("nypix");
       m_pars.Prompt("binsz");
       m_pars.Prompt("xref");
       m_pars.Prompt("yref");
       m_pars.Prompt("axisrot");
       m_pars.Prompt("proj");
      } else if (cmap == "hpx") {
       m_pars.Prompt("hpx_ordering_scheme");
       m_pars.Prompt("hpx_order");
     } 
     m_pars.Prompt("coordsys");
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
   
   if (check_cmap_null(cmap_file)) {
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

