/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube.
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.31 2005/03/08 06:22:09 jchiang Exp $
 */

#include <cstdlib>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.31 2005/03/08 06:22:09 jchiang Exp $
 */
class ExposureCube : public st_app::StApp {
public:
   ExposureCube() : st_app::StApp(), 
                    m_pars(st_app::StApp::getParGroup("gtlivetimecube")), 
                    m_exposure(0), m_roiCuts(0) {}
   virtual ~ExposureCube() throw() {
      try {
         delete m_exposure;
         delete m_roiCuts;
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
    }
   virtual void run();
   virtual void banner() const {}
private:
   st_app::AppParGroup & m_pars;
   Likelihood::LikeExposure * m_exposure;
   Likelihood::RoiCuts * m_roiCuts;
   void promptForParameters();
   void readRoiCuts();
   void createDataCube();
};

st_app::StAppFactory<ExposureCube> myAppFactory("gtlivetimecube");

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
   m_exposure->write(output_file);
   std::auto_ptr<tip::Table> 
      table(tip::IFileSvc::instance().editTable(output_file, "Exposure"));
   m_roiCuts->writeDssKeywords(table->getHeader());
   m_roiCuts->writeGtiExtension(output_file);
}

void ExposureCube::promptForParameters() {
   m_pars.Prompt("evfile");
   m_pars.Prompt("scfile");
   m_pars.Prompt("outfile");
   m_pars.Prompt("cos_theta_step");
   m_pars.Prompt("pixel_size");
   m_pars.Save();
}

void ExposureCube::readRoiCuts() {
   std::string event_file = m_pars["evfile"];
   m_roiCuts = new Likelihood::RoiCuts();
   m_roiCuts->readCuts(m_pars["evfile"], "EVENTS", false);
}

void ExposureCube::createDataCube() {
   std::vector<std::pair<double, double> > timeCuts;
   m_roiCuts->getTimeCuts(timeCuts);
   m_exposure = new Likelihood::LikeExposure(m_pars["pixel_size"], 
                                             m_pars["cos_theta_step"],
                                             m_roiCuts->timeRangeCuts(),
                                             m_roiCuts->gtis());
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
