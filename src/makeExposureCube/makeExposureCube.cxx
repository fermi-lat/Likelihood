/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube.
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.4 2004/04/05 22:01:57 jchiang Exp $
 */

#include <sstream>

#include "tip/IFileSvc.h"
#include "st_app/IApp.h"
#include "hoops/hoops_prompt_group.h"

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
 * $Header$
 */
class ExposureCube {
public:
   ExposureCube(hoops::IParGroup & pars) : m_pars(pars), m_exposure(0) {}
   virtual ~ExposureCube() {
      delete m_exposure;
   }
   virtual void run();
private:
   hoops::IParGroup & m_pars;
   Likelihood::LikeExposure * m_exposure;
   void createDataCube();
   void addRoiHistory(map_tools::ExposureHyperCube & cube);
};

/**
 * @class app
 * @brief Class (and object declaration) of boiler-plate code expected 
 * by st_app.
 */
class app : public st_app::IApp {
public:
   app() : st_app::IApp("makeExposureCube") {}
   virtual ~app() throw() {}
   virtual void run() {
      hoopsPrompt();
      hoopsSave();
      hoops::IParGroup & pars = hoopsGetParGroup();
      ExposureCube my_cube(pars);
      my_cube.run();
   }
} my_app;

void ExposureCube::run() {
   createDataCube();
   map_tools::ExposureHyperCube cube(*m_exposure, m_pars["Output file"]);
   addRoiHistory(cube);
}

void ExposureCube::createDataCube() {
   m_exposure = new Likelihood::LikeExposure(m_pars["pixel size"], 
                                             m_pars["cos_theta step"], 
                                             m_pars["ROI_file"]);
   tip::Table * scData = 
      tip::IFileSvc::instance().editTable(m_pars["Spacecraft file"], "Ext1");
   m_exposure->load(scData);
}

void ExposureCube::addRoiHistory(map_tools::ExposureHyperCube & cube) {
   Likelihood::RoiCuts * roiCuts = Likelihood::RoiCuts::instance();
   std::ostringstream roi_xml;
   roiCuts->writeXml(roi_xml);
   cube.setKey("HISTORY", roi_xml.str());
}
