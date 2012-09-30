/**
 * @file gtsrcmaps.cxx
 * @brief Compute SourceMaps for use by binned likelihood. Inputs include
 * a counts map and a source model xml file.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/gtsrcmaps/gtsrcmaps.cxx,v 1.41 2012/04/17 20:28:14 jchiang Exp $
 */

#include <cstdlib>

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
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

std::string gtsrcmaps::s_cvs_id("$Name:  $");

void gtsrcmaps::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void gtsrcmaps::run() {
   m_pars.Prompt();
   m_pars.Save();
   m_helper = new AppHelpers(&m_pars, "BINNED");
   m_helper->checkOutputFile();
   m_helper->checkTimeCuts(m_pars["cmap"], "",
                           m_pars["expcube"], "Exposure");

   std::string expCubeFile = m_pars["expcube"];
   ExposureCube & expCube = 
      const_cast<ExposureCube &>(m_helper->observation().expCube());
   expCube.readExposureCube(expCubeFile);

   std::string cntsMapFile = m_pars["cmap"];
   std::string irfs = m_pars["irfs"];
   dataSubselector::Cuts::checkIrfs(cntsMapFile, "", irfs);
   dataSubselector::Cuts my_cuts(cntsMapFile, "", false);
   CountsMap dataMap(cntsMapFile);
   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

   double ra, dec;
   getRefCoord(cntsMapFile, ra, dec);

   RoiCuts &roiCuts = const_cast<RoiCuts &>(m_helper->observation().roiCuts());
   roiCuts.setCuts(ra, dec, 20., energies.front(), energies.back());

   std::string binnedMap = m_pars["bexpmap"];
   AppHelpers::checkExposureMap(m_pars["cmap"], m_pars["bexpmap"]);

   if (!st_facilities::Util::fileExists(binnedMap)) {
      std::ostringstream message;
      message << "Binned exposure map file named "
              << binnedMap << " does not exist.";
      throw std::runtime_error(message.str());
   } else {
      bool enforce_boundaries = m_pars["emapbnds"];
      m_helper->observation().bexpmap().setBoundaryFlag(enforce_boundaries);
   }
   bool computePointSources = AppHelpers::param(m_pars, "ptsrc", true);
   bool psf_corrections = AppHelpers::param(m_pars, "psfcorr", true);
   bool perform_convolution = AppHelpers::param(m_pars, "convol", true);

   bool resample = m_pars["resample"];
   int resamp_factor = m_pars["rfactor"];
   double minbinsz = m_pars["minbinsz"];

   m_binnedLikelihood = 
      new BinnedLikelihood(dataMap, m_helper->observation(),
                           cntsMapFile, computePointSources, psf_corrections,
                           perform_convolution, resample, resamp_factor,
                           minbinsz);

   std::string srcModelFile = m_pars["srcmdl"];
   bool loadMaps, createAllMaps;
   try {
/// Turn off loading of maps when xml file is read in.  Instead, read maps
/// when they are first accessed by MapBase objects.
      m_binnedLikelihood->readXml(srcModelFile, m_helper->funcFactory(), false,
                                  computePointSources, loadMaps=false,
                                  createAllMaps=true);
   } catch(std::runtime_error & eObj) {
      std::string message("Request for exposure at a sky position "
                          "that is outside of the map boundaries.");
      if (st_facilities::Util::expectedException(eObj, message)) {
         std::ostringstream app_message;
         app_message << "\n" << message << "\n\n"
                     << "The contribution of the diffuse source outside of "
                     << "the exposure \nand counts map boundaries is being "
                     << "computed to account for PSF \nleakage into the "
                     << "analysis region.  To handle this, use an all-sky\n"
                     << "binned exposure map.  Alternatively, to neglect "
                     << "contributions \n"
                     << "outside of the counts map region, use the "
                     << "emapbnds=no option when \nrunning gtsrcmaps.";
         throw std::runtime_error(app_message.str());
      }
      throw;
   }

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
