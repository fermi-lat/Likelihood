/**
 * @file AppHelpers.cxx
 * @brief Class of "helper" methods for Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AppHelpers.cxx,v 1.27 2005/03/04 22:08:26 jchiang Exp $
 */

#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "irfLoader/Loader.h"
#include "irfInterface/IrfsFactory.h"

#include "st_facilities/Util.h"

#include "dataSubselector/CutBase.h"
#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BandFunction.h"
#include "Likelihood/EventContainer.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/MapCubeFunction.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SpatialMap.h"

using irfInterface::IrfsFactory;

namespace Likelihood {

AppHelpers::AppHelpers(st_app::AppParGroup & pars)
   : m_pars(pars), m_funcFactory(0) {
   prepareFunctionFactory();
   createResponseFuncs();

   m_roiCuts = new RoiCuts();
   m_scData = new ScData();
   m_expCube = new ExposureCube();
   m_expMap = new ExposureMap();
   m_eventCont = new EventContainer(*m_respFuncs, *m_roiCuts, *m_scData);
   m_observation = new Observation(m_respFuncs,
                                   m_scData,
                                   m_roiCuts,
                                   m_expCube,
                                   m_expMap,
                                   m_eventCont);
}

AppHelpers::~AppHelpers() {
   delete m_funcFactory;
   delete m_roiCuts;
   delete m_scData;
   delete m_expCube;
   delete m_expMap;
   delete m_eventCont;
   delete m_respFuncs;
   delete m_observation;
}

optimizers::FunctionFactory & AppHelpers::funcFactory() {
   return *m_funcFactory;
}

void AppHelpers::prepareFunctionFactory() {
   bool makeClone(false);
   m_funcFactory = new optimizers::FunctionFactory;
   m_funcFactory->addFunc("SkyDirFunction", new SkyDirFunction(), makeClone);
   m_funcFactory->addFunc("SpatialMap", new SpatialMap(), makeClone);
   m_funcFactory->addFunc("BandFunction", new BandFunction(), makeClone);
   m_funcFactory->addFunc("MapCubeFunction", new MapCubeFunction(), makeClone);
}

void AppHelpers::setRoi(const std::string & filename,
                        const std::string & ext, bool strict) {
   RoiCuts & roiCuts = const_cast<RoiCuts &>(m_observation->roiCuts());
   if (filename != "") {
      roiCuts.readCuts(filename, ext, strict);
      return;
   }
   std::string event_file = m_pars["evfile"];
   roiCuts.readCuts(event_file, "EVENTS", strict);
}

void AppHelpers::readScData() {
   std::string scFile = m_pars["scfile"];
   st_facilities::Util::file_ok(scFile);
   st_facilities::Util::resolve_fits_files(scFile, m_scFiles);
   std::vector<std::string>::const_iterator scIt = m_scFiles.begin();
   for ( ; scIt != m_scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      m_scData->readData(*scIt);
   }
}

void AppHelpers::readExposureMap() {
   std::string exposureFile = m_pars["exposure_map_file"];
   if (exposureFile != "none") {
      st_facilities::Util::file_ok(exposureFile);
      m_expMap->readExposureFile(exposureFile);
   }
}

void AppHelpers::createResponseFuncs() {
   m_respFuncs = new ResponseFunctions();
   std::string responseFuncs = m_pars["rspfunc"];
   m_respFuncs->load(responseFuncs);
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

void AppHelpers::checkCuts(const std::string & file1,
                           const std::string & ext1,
                           const std::string & file2,
                           const std::string & ext2) {
   dataSubselector::Cuts cuts1(file1, ext1, false);
   dataSubselector::Cuts cuts2(file2, ext2, false);
   if (!(cuts1 == cuts2)) {
// Try comparing output streams.
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

void AppHelpers::checkTimeCuts(const std::string & file1,
                               const std::string & ext1,
                               const std::string & file2,
                               const std::string & ext2) {
   dataSubselector::Cuts cuts1(file1, ext1, false);
   dataSubselector::Cuts cuts2(file2, ext2, false);
// This is a bit fragile as one must assume the ordering of the 
// cuts is the same for both Cuts objects.
   std::vector<const dataSubselector::CutBase *> time_cuts1;
   std::vector<const dataSubselector::CutBase *> time_cuts2;
   gatherTimeCuts(cuts1, time_cuts1);
   gatherTimeCuts(cuts2, time_cuts2);
   bool ok(true);
   if (time_cuts1.size() == time_cuts2.size()) {
      for (unsigned int i = 0; i < time_cuts1.size(); i++) {
         ok = ok && *(time_cuts1[i]) == *(time_cuts2[i]);
      }
   } else {
      ok = false;
   }
   if (!ok) {
      std::ostringstream message;
      message << "AppHelpers::checkTimeCuts:\n" 
              << "Time cuts in files " << file1;
      if (ext1 != "") message << "[" << ext1 << "]";
      message << " and " << file2;
      if (ext2 != "") message << "[" << ext2 << "]";
      message << " do not agree.";
      throw std::runtime_error(message.str());
   }                
}

void AppHelpers::
gatherTimeCuts(dataSubselector::Cuts & cuts,
               std::vector<const dataSubselector::CutBase *> time_cuts) {
   for (unsigned int i = 0; i < cuts.size(); i++) {
      if ( cuts[i].type() == "GTI" || 
           (cuts[i].type() == "range" &&
            dynamic_cast<dataSubselector::RangeCut &>(
               const_cast<dataSubselector::CutBase &>(cuts[i])).colname() 
            == "TIME") ) {
         time_cuts.push_back(&cuts[i]);
      }
   }
}

} // namespace Likelihood
