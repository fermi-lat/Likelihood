/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube.
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.16 2004/12/06 20:20:30 jchiang Exp $
 */

#include <cstdlib>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

#include "map_tools/ExposureHyperCube.h"

#include "Likelihood/LikeExposure.h"
#include "Likelihood/RoiCuts.h"

#include "Verbosity.h"

/**
 * @class ExposureCube
 * @brief Class to encapsulate methods for creating an exposure
 * hypercube in (ra, dec, cos_theta) using the LikeExposure class.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.16 2004/12/06 20:20:30 jchiang Exp $
 */
class ExposureCube : public st_app::StApp {
public:
   ExposureCube() : st_app::StApp(), 
                    m_pars(st_app::StApp::getParGroup("makeExposureCube")), 
                    m_exposure(0) {}
   virtual ~ExposureCube() throw() {
      try {
         delete m_exposure;
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
    }
   virtual void run();
private:
   st_app::AppParGroup & m_pars;
   Likelihood::LikeExposure * m_exposure;
   void promptForParameters();
   void readRoiCuts() const;
   void createDataCube();
};

st_app::StAppFactory<ExposureCube> myAppFactory;

void ExposureCube::run() {
   promptForParameters();
   readRoiCuts();
   std::string output_file = m_pars["outfile"];
   if (st_facilities::Util::fileExists(output_file)) {
      if (m_pars["clobber"]) {
         std::remove(output_file.c_str());
      } else {
         std::cout << "Output file " << output_file 
                   << " already exists and you have set 'clobber' to 'no'.\n"
                   << "Please provide a different output file name." 
                   << std::endl;
         std::exit(1);
      }
   }
   Likelihood::Verbosity::instance(m_pars["chatter"]);
   createDataCube();
   map_tools::ExposureHyperCube cube(*m_exposure, output_file);
   cube.save();
   tip::Image * image
      = tip::IFileSvc::instance().editImage(output_file, "hypercube");
   Likelihood::RoiCuts::instance()->writeDssKeywords(image->getHeader());
   Likelihood::RoiCuts::instance()->writeGtiExtension(output_file);
}

void ExposureCube::promptForParameters() {
   m_pars.Prompt("evfile");
   std::string event_file = m_pars["evfile"];
   if (event_file == "none" || event_file == "") {
      m_pars.Prompt("ROI_file");
      std::string Roi_file = m_pars["ROI_file"];
      if (!st_facilities::Util::fileExists(Roi_file)) {
         throw std::runtime_error("ROI file " + Roi_file +
                                  " does not exist.  Please specify a " +
                                  std::string("valid event file or ") +
                                  "an ROI file.");
      }
   } else {
      st_facilities::Util::file_ok(m_pars["evfile"]);
   }
   m_pars.Prompt("scfile");
   m_pars.Prompt("outfile");
   m_pars.Save();
}

void ExposureCube::readRoiCuts() const {
   std::string event_file = m_pars["evfile"];
   if (event_file == "none" || event_file == "") {
      std::string roi_file = m_pars["ROI_file"];
      Likelihood::RoiCuts::setCuts(roi_file);
   } else {
      Likelihood::RoiCuts::instance()->readCuts(m_pars["evfile"]);
   }
}

void ExposureCube::createDataCube() {
   m_exposure = new Likelihood::LikeExposure(m_pars["pixel_size"], 
                                             m_pars["cos_theta_step"]);
   std::string scFile = m_pars["scfile"];
   st_facilities::Util::file_ok(scFile);
   std::vector<std::string> scFiles;
   st_facilities::Util::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      if (Likelihood::print_output()) {
         std::cerr << "Working on file " << *scIt << std::endl;
      }
      tip::Table * scData = 
         tip::IFileSvc::instance().editTable(*scIt, m_pars["sctable"]);
      m_exposure->load(scData, Likelihood::print_output());
      delete scData;
   }
}

