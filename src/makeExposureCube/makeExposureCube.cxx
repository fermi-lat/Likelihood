/** 
 * @file makeExposureCube.cxx
 * @brief Create an Exposure hypercube. This program is based on
 * map_tools/exposure_cube.cxx
 * @author J. Chiang
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/makeExposureCube/makeExposureCube.cxx,v 1.3 2004/04/05 18:31:11 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

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
   void createDataCube();
   void addRoiHistory(map_tools::ExposureHyperCube & cube);
   void tearDown();

};

st_app::IApp * my_application = new ExposureCube();

void ExposureCube::setUp() {
#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
   hoopsPrompt();
   hoopsSave();
}

void ExposureCube::tearDown() {
   delete m_exposure;
}

void ExposureCube::createDataCube() {
   hoops::IParGroup & pars = hoopsGetParGroup();
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
   hoops::IParGroup & pars = hoopsGetParGroup();
   setUp();
   createDataCube();
   map_tools::ExposureHyperCube cube(*m_exposure, pars["Output file"]);
   addRoiHistory(cube);
   tearDown();
}
