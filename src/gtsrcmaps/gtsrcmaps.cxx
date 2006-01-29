/**
 * @file gtsrcmaps.cxx
 * @brief Compute SourceMaps for use by binned likelihood. Inputs include
 * a counts map and a source model xml file.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/gtsrcmaps/gtsrcmaps.cxx,v 1.21 2005/06/04 20:05:09 jchiang Exp $
 */

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

#include "dataSubselector/Cuts.h"

#include "irfInterface/Irfs.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/RoiCuts.h"

#include "Verbosity.h"

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
   virtual void banner() const;
private:
   AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;

   BinnedLikelihood * m_binnedLikelihood;

   void getRefCoord(const std::string & countsMapFile, 
                    double & ra, double & dec) const;

   static std::string s_cvs_id;

};

st_app::StAppFactory<gtsrcmaps> myAppFactory("gtsrcmaps");

gtsrcmaps::gtsrcmaps() 
   : st_app::StApp(), m_helper(0),
     m_pars(st_app::StApp::getParGroup("gtsrcmaps")),
     m_binnedLikelihood(0) {
   setVersion(s_cvs_id);
}

std::string gtsrcmaps::s_cvs_id("$Name$");

void gtsrcmaps::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void gtsrcmaps::run() {
   m_pars.Prompt();
   m_pars.Save();
   Likelihood::Verbosity::instance(m_pars["chatter"]);
   m_helper = new AppHelpers(&m_pars);
   m_helper->readScData();
   m_helper->checkOutputFile();
   m_helper->checkTimeCuts(m_pars["counts_map_file"], "",
                           m_pars["exposure_cube_file"], "Exposure");

   std::string expCubeFile = m_pars["exposure_cube_file"];
   ExposureCube & expCube = 
      const_cast<ExposureCube &>(m_helper->observation().expCube());
   expCube.readExposureCube(expCubeFile);

   std::string cntsMapFile = m_pars["counts_map_file"];
   dataSubselector::Cuts my_cuts(cntsMapFile, "", false);
   CountsMap dataMap(cntsMapFile);
   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

   double ra, dec;
   getRefCoord(cntsMapFile, ra, dec);

   RoiCuts &roiCuts = const_cast<RoiCuts &>(m_helper->observation().roiCuts());
   roiCuts.setCuts(ra, dec, 20., energies.front(), energies.back());

   std::string binnedMap = m_pars["binned_exposure_map"];
   if (binnedMap != "none" && binnedMap != "") {
      SourceMap::setBinnedExposure(binnedMap);
   }
   bool computePointSources =
      AppHelpers::param(m_pars, "compute_point_sources", true);
   bool psf_corrections =
      AppHelpers::param(m_pars, "apply_psf_corrections", true);
   m_binnedLikelihood = 
      new BinnedLikelihood(dataMap, m_helper->observation(),
                           cntsMapFile, computePointSources, psf_corrections);

   std::string srcModelFile = m_pars["source_model_file"];
   m_binnedLikelihood->readXml(srcModelFile, m_helper->funcFactory(), false);

   std::string srcMapsFile = m_pars["outfile"];

   dataMap.writeOutput("gtsrcmaps", srcMapsFile);

   m_binnedLikelihood->saveSourceMaps(srcMapsFile);

   std::auto_ptr<tip::Image>
      image(tip::IFileSvc::instance().editImage(srcMapsFile, ""));
   my_cuts.writeDssKeywords(image->getHeader());
   my_cuts.writeGtiExtension(srcMapsFile);
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
      throw std::runtime_error("gtsrcmaps::getRefCoord: " + 
                               std::string("Unknown coordinate system in ") +
                               countsMapFile);
   }
}
