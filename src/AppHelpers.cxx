/**
 * @file AppHelpers.cxx
 * @brief Class of "helper" methods for Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AppHelpers.cxx,v 1.13 2004/12/07 01:30:10 jchiang Exp $
 */

#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "irfLoader/Loader.h"
#include "irfInterface/IrfsFactory.h"

#include "st_facilities/Util.h"

#include "dataSubselector/Cuts.h"

#include "Likelihood/ExposureMap.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SpatialMap.h"

#include "Likelihood/AppHelpers.h"

using irfInterface::IrfsFactory;

namespace Likelihood {

AppHelpers::AppHelpers(st_app::AppParGroup & pars)
   : m_pars(pars), m_funcFactory(0) {
   prepareFunctionFactory();
   createResponseFuncs();
}

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
   std::string event_file = m_pars["evfile"];
   if (event_file == "none" || event_file == "") {
      std::string roi_file = m_pars["ROI_file"];
      st_facilities::Util::file_ok(roi_file);
      RoiCuts::setCuts(roi_file);
   } else {
      RoiCuts::instance()->readCuts(m_pars["evfile"]);
   }
}

void AppHelpers::readScData() {
   std::string scFile = m_pars["scfile"];
   st_facilities::Util::file_ok(scFile);
   st_facilities::Util::resolve_fits_files(scFile, m_scFiles);
   std::vector<std::string>::const_iterator scIt = m_scFiles.begin();
   for ( ; scIt != m_scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      ScData::readData(*scIt);
   }
}

void AppHelpers::readExposureMap() {
   std::string exposureFile = m_pars["exposure_map_file"];
   if (exposureFile != "none") {
      st_facilities::Util::file_ok(exposureFile);
      ExposureMap::readExposureFile(exposureFile);
   }
}

void AppHelpers::createResponseFuncs() {
   irfLoader::Loader::go();
   IrfsFactory * myFactory = IrfsFactory::instance();

   std::string responseFuncs = m_pars["rspfunc"];

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

void AppHelpers::checkOutputFile(bool clobber, const std::string & file) {
   if (!clobber) {
      if (st_facilities::Util::fileExists(file)) {
         std::cout << "Output file " << file 
                   << " already exists and you have set 'clobber' to 'no'.\n"
                   << "Please provide a different output file name."
                   << std::endl;
         std::exit(1);
      }
   }
}

void AppHelpers::checkCuts(const std::string & file1, const std::string ext1,
                           const std::string & file2, const std::string ext2) {
   dataSubselector::Cuts cuts1(file1, ext1, false);
   dataSubselector::Cuts cuts2(file2, ext2, false);
   if (!(cuts1 == cuts2)) {
// try comparing output streams
      std::ostringstream c1, c2;
      cuts1.writeCuts(c1);
      cuts2.writeCuts(c2);
      if (c1.str() != c2.str()) {
         std::ostringstream message;
         message << "AppHelpers::checkCuts:\n" 
                 << "DSS keywords in " << file1;
         if (ext1 != "") message << "[" << ext1 << "] ";
         message << "do not match those in " << file2;
         if (ext2 != "") message << "[" << ext2 << "] ";
         throw std::runtime_error(message.str());
      }
   }
}

} // namespace Likelihood
