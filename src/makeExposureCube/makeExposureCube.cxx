/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube.
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.10 2004/09/03 14:11:03 jchiang Exp $
 */

#include <cstdlib>

#include <iostream>
#include <sstream>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

#include "map_tools/ExposureHyperCube.h"

#include "Likelihood/RoiCuts.h"
#include "Likelihood/LikeExposure.h"

/**
 * @class ExposureCube
 * @brief Class to encapsulate methods for creating an exposure
 * hypercube in (ra, dec, cos_theta) using the LikeExposure class.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.10 2004/09/03 14:11:03 jchiang Exp $
 */
class ExposureCube : public st_app::StApp {
public:
   ExposureCube() : st_app::StApp(), 
                    m_pars(st_app::StApp::getParGroup("makeExposureCube")), 
                    m_exposure(0) {
      try {
         m_pars.Prompt();
         m_pars.Save();
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
         std::exit(1);
      } catch (...) {
         std::cerr << "Caught unknown exception in ExposureCube constructor." 
                   << std::endl;
         std::exit(1);
      }
   }

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
   void createDataCube();
   void addRoiHistory(map_tools::ExposureHyperCube & cube);
};

st_app::StAppFactory<ExposureCube> myAppFactory;

void ExposureCube::run() {
   createDataCube();
   map_tools::ExposureHyperCube cube(*m_exposure, m_pars["Output_file"]);
   addRoiHistory(cube);
}

void ExposureCube::createDataCube() {
   m_exposure = new Likelihood::LikeExposure(m_pars["pixel_size"], 
                                             m_pars["cos_theta_step"], 
                                             m_pars["ROI_file"]);
   std::string scFile = m_pars["Spacecraft_file"];
   st_facilities::Util::file_ok(scFile);
   std::vector<std::string> scFiles;
   st_facilities::Util::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      std::cerr << "Working on file " << *scIt << std::endl;
      tip::Table * scData = 
         tip::IFileSvc::instance().editTable(*scIt, "Ext1");
      m_exposure->load(scData);
      delete scData;
   }
}

void ExposureCube::addRoiHistory(map_tools::ExposureHyperCube & cube) {
   Likelihood::RoiCuts * roiCuts = Likelihood::RoiCuts::instance();
   std::ostringstream roi_xml;
   roiCuts->writeXml(roi_xml);
   cube.setKey("HISTORY", roi_xml.str());
}
