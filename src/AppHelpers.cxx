/**
 * @file AppHelpers.cxx
 * @brief Class of "helper" methods for Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AppHelpers.cxx,v 1.2 2004/05/24 23:51:31 jchiang Exp $
 */

#include <map>
#include <stdexcept>
#include <vector>

//#include "latResponse/IrfsFactory.h"
#include "dc1Response/loadIrfs.h"
#include "g25Response/loadIrfs.h"
#include "irfInterface/IrfsFactory.h"

#include "Likelihood/ExposureMap.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/Util.h"

#include "Likelihood/AppHelpers.h"

//using latResponse::irfsFactory;
using irfInterface::IrfsFactory;

namespace Likelihood {

optimizers::FunctionFactory & AppHelpers::funcFactory() {
   return *m_funcFactory;
}

void AppHelpers::prepareFunctionFactory() {
   bool makeClone(false);
   m_funcFactory = new optimizers::FunctionFactory;
   m_funcFactory->addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
   m_funcFactory->addFunc("SpatialMap", new SpatialMap(), makeClone);
}

void AppHelpers::setRoi() {
   std::string roiCutsFile = m_pars["ROI_cuts_file"];
   Util::file_ok(roiCutsFile);
   RoiCuts::setCuts(roiCutsFile);
}

void AppHelpers::readScData() {
   std::string scFile = m_pars["Spacecraft_file"];
   Util::file_ok(scFile);
   long scHdu = m_pars["Spacecraft_file_hdu"];
   std::vector<std::string> scFiles;
   Util::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      Util::file_ok(*scIt);
      ScData::readData(*scIt, scHdu);
   }
}

void AppHelpers::readExposureMap() {
   std::string exposureFile = m_pars["Exposure_map_file"];
   if (exposureFile != "none") {
      Util::file_ok(exposureFile);
      ExposureMap::readExposureFile(exposureFile);
   }
}

void AppHelpers::createResponseFuncs() {
   dc1Response::loadIrfs();
   g25Response::loadIrfs();

   IrfsFactory * myFactory = IrfsFactory::instance();

   std::string responseFuncs = m_pars["Response_functions"];
   std::map< std::string, std::vector<std::string> > responseIds;
   responseIds["FRONT"].push_back("DC1::Front");
   responseIds["BACK"].push_back("DC1::Back");
   responseIds["FRONT/BACK"].push_back("DC1::Front");
   responseIds["FRONT/BACK"].push_back("DC1::Back");
   responseIds["GLAST25"].push_back("Glast25::Front");
   responseIds["GLAST25"].push_back("Glast25::Back");
   responseIds["GLAST25_10"].push_back("Glast25::Front_10");
   responseIds["GLAST25_10"].push_back("Glast25::Back_10");
   if (responseIds.count(responseFuncs)) {
      std::vector<std::string> &resps = responseIds[responseFuncs];
      for (unsigned int i = 0; i < resps.size(); i++) {
//          ResponseFunctions::addRespPtr(i, irfsFactory().create(resps[i]));
         ResponseFunctions::addRespPtr(i, myFactory->create(resps[i]));
      }
   } else {
      throw std::invalid_argument("Invalid response function choice: "
                                  + responseFuncs);
   }
}

} // namespace Likelihood
