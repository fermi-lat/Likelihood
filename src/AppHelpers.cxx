/**
 * @file AppHelpers.cxx
 * @brief Class of "helper" methods for Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AppHelpers.cxx,v 1.7 2004/08/25 15:27:32 jchiang Exp $
 */

#include <map>
#include <stdexcept>
#include <vector>

#include "irfLoader/Loader.h"
#include "irfInterface/IrfsFactory.h"

#include "st_facilities/Util.h"

#include "Likelihood/ExposureMap.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SpatialMap.h"

#include "Likelihood/AppHelpers.h"

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
   st_facilities::Util::file_ok(roiCutsFile);
   RoiCuts::setCuts(roiCutsFile);
}

void AppHelpers::readScData() {
   std::string scFile = m_pars["Spacecraft_file"];
   st_facilities::Util::file_ok(scFile);
   long scHdu = m_pars["Spacecraft_file_hdu"];
   st_facilities::Util::resolve_fits_files(scFile, m_scFiles);
   std::vector<std::string>::const_iterator scIt = m_scFiles.begin();
   for ( ; scIt != m_scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      ScData::readData(*scIt, scHdu);
   }
}

void AppHelpers::readExposureMap() {
   std::string exposureFile = m_pars["Exposure_map_file"];
   if (exposureFile != "none") {
      st_facilities::Util::file_ok(exposureFile);
      ExposureMap::readExposureFile(exposureFile);
   }
}

void AppHelpers::createResponseFuncs() {
   irfLoader::Loader::go();
   IrfsFactory * myFactory = IrfsFactory::instance();

   std::string responseFuncs = m_pars["Response_functions"];

   typedef std::map< std::string, std::vector<std::string> > respMap;
   const respMap & responseIds = irfLoader::Loader::respIds();
   respMap::const_iterator it;
   if ( (it = responseIds.find(responseFuncs)) != responseIds.end() ) {
      const std::vector<std::string> & resps = it->second;
      for (unsigned int i = 0; i < resps.size(); i++) {
         ResponseFunctions::addRespPtr(i, myFactory->create(resps[i]));
      }
   } else {
      throw std::invalid_argument("Invalid response function choice: "
                                  + responseFuncs);
   }
}

} // namespace Likelihood
