#include <cstdlib>

#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "fitsio.h"

#include "facilities/Util.h"

#include "st_facilities/Util.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Header.h"

#include "astro/SkyDir.h"

#include "irfInterface/Irfs.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/BinnedLikelihood.h"

using namespace Likelihood;

class gtsrcmaps : public st_app::StApp {
public:
   gtsrcmaps();
   virtual ~gtsrcmaps() throw() {
      try {
         delete m_helper;
         delete m_binnedLikelihood;
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
         std::cerr << "gtsrcmaps::~gtsrcmaps: "
                   << "unknown exception encountered."
                   << std::endl;
      }
   }
   virtual void run();
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;

   BinnedLikelihood * m_binnedLikelihood;

   void getRefCoord(const std::string & countsMapFile, 
                    double & ra, double & dec) const;

};

st_app::StAppFactory<gtsrcmaps> myAppFactory;

gtsrcmaps::gtsrcmaps() 
   : st_app::StApp(), m_helper(0),
     m_pars(st_app::StApp::getParGroup("gtsrcmaps")),
     m_binnedLikelihood(0) {
   try {
      m_pars.Prompt();
      m_pars.Save();
      m_helper = new AppHelpers(m_pars);
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
      std::exit(1);
   } catch (...) {
      std::cerr << "Caught unknown exception in gtsrcmaps constructor." 
                << std::endl;
      std::exit(1);
   }
}

void gtsrcmaps::run() {
   std::string expcube_file = m_pars["exposure_cube_file"];
   if (expcube_file == "none") {
      throw std::runtime_error("Please specify an exposure cube file.");
   }
   ExposureCube::readExposureCube(expcube_file);
   std::string cntsMapFile = m_pars["counts_map_file"];

   CountsMap dataMap(cntsMapFile);
   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

   double ra, dec;
   getRefCoord(cntsMapFile, ra, dec);
// NB: energies from EBOUNDS are in keV.
   RoiCuts::instance()->setCuts(ra, dec, 20., energies.front()/1e3, 
                                energies.back()/1e3);

   m_binnedLikelihood = new BinnedLikelihood(dataMap, cntsMapFile);

   std::string srcModelFile = m_pars["source_model_file"];
   m_binnedLikelihood->readXml(srcModelFile, m_helper->funcFactory(), false);

   std::string srcMapsFile = m_pars["output_file"];

   dataMap.writeOutput("gtsrcmaps", srcMapsFile);

   m_binnedLikelihood->saveSourceMaps(srcMapsFile);
}

void gtsrcmaps::getRefCoord(const std::string & countsMapFile, 
                            double & ra, double & dec) const {
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(countsMapFile, ""));
   const tip::Header & header = image->getHeader();
   std::string ctype;
   header["CTYPE1"].get(ctype);
   if (ctype.find("RA") != std::string::npos) {
      header["CRVAL1"].get(ra);
      header["CRVAL2"].get(dec);
   } else if (ctype.find("GLON") != std::string::npos) {
      double glon, glat;
      header["CRVAL1"].get(glon);
      header["CRVAL2"].get(glat);
      astro::SkyDir dir(glon, glat, astro::SkyDir::GALACTIC);
      ra = dir.ra();
      dec = dir.dec();
   } else {
      throw std::runtime_error("gtsrcmaps: Unknown coordinate system in "
                               + countsMapFile);
   }
}
