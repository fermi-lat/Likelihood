/**
 * @file AppBase.cxx
 * @brief Base class for Likelihood applications.
 * @author J. Chiang
 *
 * $Header$
 */

#include <map>
#include <stdexcept>
#include <vector>

#include "hoops/hoops_prompt_group.h"

#include "latResponse/IrfsFactory.h"

#include "Likelihood/ExposureMap.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/Util.h"

#include "Likelihood/AppBase.h"

using latResponse::irfsFactory;

namespace Likelihood {

void AppBase::setUp() {
#ifdef TRAP_FPE
   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
   hoopsPrompt();
   hoopsSave();
   prepareFunctionFactory();
   setRoi();
   readScData();
   createResponseFuncs();
}

void AppBase::prepareFunctionFactory() {
   bool makeClone(false);
   m_funcFactory.addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
   m_funcFactory.addFunc("SpatialMap", new SpatialMap(), makeClone);
}

void AppBase::tearDown() {
}

void AppBase::setRoi() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string roiCutsFile = pars["ROI_cuts_file"];
   Util::file_ok(roiCutsFile);
   RoiCuts::setCuts(roiCutsFile);
}

void AppBase::readScData() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string scFile = pars["Spacecraft_file"];
   Util::file_ok(scFile);
   long scHdu = pars["Spacecraft_file_hdu"];
   std::vector<std::string> scFiles;
   Util::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      Util::file_ok(*scIt);
      ScData::readData(*scIt, scHdu);
   }
}

void AppBase::readExposureMap() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string exposureFile = pars["Exposure_map_file"];
   if (exposureFile != "none") {
      Util::file_ok(exposureFile);
      ExposureMap::readExposureFile(exposureFile);
   }
}

void AppBase::createResponseFuncs() {
   hoops::IParGroup & pars = hoopsGetParGroup();
   std::string responseFuncs = pars["Response_functions"];
   std::map< std::string, std::vector<std::string> > responseIds;
   responseIds["FRONT"].push_back("DC1::Front");
   responseIds["BACK"].push_back("DC1::Back");
   responseIds["FRONT/BACK"].push_back("DC1::Front");
   responseIds["FRONT/BACK"].push_back("DC1::Back");
   responseIds["GLAST25"].push_back("Glast25::Front");
   responseIds["GLAST25"].push_back("Glast25::Back");
   if (responseIds.count(responseFuncs)) {
      std::vector<std::string> &resps = responseIds[responseFuncs];
      for (unsigned int i = 0; i < resps.size(); i++) {
         ResponseFunctions::addRespPtr(i, irfsFactory().create(resps[i]));
      }
   } else {
      throw std::invalid_argument("Invalid response function choice: "
                                  + responseFuncs);
   }
}

} // namespace Likelihood
