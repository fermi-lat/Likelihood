/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube. This program is based on
 * map_tools/exposure_cube.cxx
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.2 2004/04/04 01:22:25 jchiang Exp $
 */

#include <sstream>

#include "facilities/Util.h"
#include "tip/IFileSvc.h"
#include "hoops/hoops_prompt_group.h"
#include "st_app/IApp.h"
#include "map_tools/ExposureHyperCube.h"

#include "Likelihood/RoiCuts.h"
#include "Likelihood/LikeExposure.h"

class ExposureCube : public st_app::IApp {

public:

   ExposureCube() : st_app::IApp("makeExposureCube"), m_exposure(0) {}

   virtual ~ExposureCube() throw() {}

   virtual void run();

private:

   Likelihood::LikeExposure * m_exposure;

   void setUp();
   void createDataCube(hoops::IParGroup & pars);
   void addRoiHistory(map_tools::ExposureHyperCube & cube);
   void tearDown();

};

st_app::IApp * my_application = new ExposureCube();

void ExposureCube::setUp() {
   hoopsPrompt();
   hoopsSave();
}

void ExposureCube::tearDown() {
   delete m_exposure;
}

void ExposureCube::createDataCube(hoops::IParGroup & pars) {
   m_exposure = new Likelihood::LikeExposure(pars["pixel size"], 
                                             pars["cos_theta step"], 
                                             pars["ROI_file"]);
   tip::Table * scData = 
      tip::IFileSvc::instance().editTable(pars["Spacecraft file"], "Ext1");
   m_exposure->load(scData);
}

void ExposureCube::addRoiHistory(map_tools::ExposureHyperCube & cube) {
   Likelihood::RoiCuts * roiCuts = Likelihood::RoiCuts::instance();
   std::ostringstream roi_xml;
   roiCuts->writeXml(roi_xml);
   cube.setKey("HISTORY", roi_xml.str());
}

void ExposureCube::run() {
   setUp();
   hoops::IParGroup & pars = hoopsGetParGroup();
   createDataCube(pars);
   map_tools::ExposureHyperCube cube(*m_exposure, pars["Output file"]);
   addRoiHistory(cube);
   tearDown();
}
