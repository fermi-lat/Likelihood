/**
 * @file AppHelpers.cxx
 * @brief Class of "helper" methods for Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AppHelpers.cxx,v 1.56 2006/12/04 01:38:22 jchiang Exp $
 */

#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Util.h"

#include "dataSubselector/CutBase.h"
#include "dataSubselector/Cuts.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BandFunction.h"
#include "Likelihood/BrokenPowerLaw2.h"
#include "Likelihood/BrokenPowerLawExpCutoff.h"
#include "Likelihood/EventContainer.h"
#include "Likelihood/ExpCutoff.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/FileFunction.h"
#include "Likelihood/LogParabola.h"
#include "Likelihood/MapCubeFunction.h"
#include "Likelihood/Observation.h"
#include "Likelihood/PowerLaw2.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SpatialMap.h"

namespace {
   void getRangeBounds(const std::vector<dataSubselector::RangeCut *> & cuts,
                       double & xmin, double & xmax) {
      xmin = cuts.at(0)->minVal();
      xmax = cuts.at(0)->maxVal();
      for (size_t j = 0; j < cuts.size(); j++) {
         if (xmin < cuts.at(j)->minVal()) {
            xmin = cuts.at(j)->minVal();
         }
         if (xmax > cuts.at(j)->maxVal()) {
            xmax = cuts.at(j)->maxVal();
         }
      }
   }
}

namespace Likelihood {

AppHelpers::AppHelpers(st_app::AppParGroup * pars,
                       const std::string & analysisType) 
   : m_pars(pars), m_funcFactory(0), m_respFuncs(0) {
   prepareFunctionFactory();
   createResponseFuncs(analysisType);

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
   m_funcFactory->addFunc("LogParabola", new LogParabola(), makeClone);
   m_funcFactory->addFunc("MapCubeFunction", new MapCubeFunction(), makeClone);
   m_funcFactory->addFunc("PowerLaw2", new PowerLaw2(), makeClone);
   m_funcFactory->addFunc("BrokenPowerLaw2", new BrokenPowerLaw2(), makeClone);
   m_funcFactory->addFunc("FileFunction", new FileFunction(), makeClone);
   m_funcFactory->addFunc("ExpCutoff", new ExpCutoff(), makeClone);
   m_funcFactory->addFunc("BPLExpCutoff", new BrokenPowerLawExpCutoff(),
                          makeClone);
}

void AppHelpers::setRoi(const std::string & filename,
                        const std::string & ext, bool strict) {
   RoiCuts & roiCuts = const_cast<RoiCuts &>(m_observation->roiCuts());
   if (filename != "") {
      roiCuts.readCuts(filename, ext, strict);
      return;
   }
   st_app::AppParGroup & pars(*m_pars);
   std::string event_file = pars["evfile"];
   std::string evtable = pars["evtable"];
   std::vector<std::string> eventFiles;
   st_facilities::Util::resolve_fits_files(event_file, eventFiles);
   roiCuts.readCuts(eventFiles, evtable, strict);
}

std::string AppHelpers::responseFuncs(const std::string & file,
                                      const std::string & respBase) {
   static char * respcombos[] = {"DC2", "DC2FA", "DC2BA", "DC2FB", "DC2BB", 
                                 "DC2_A"};
   if (respBase != "DC2") {
      return respBase;
   }
// If just DC2 is given, then parse the EVENT_CLASS keyword to find
// the proper combination.
   dataSubselector::Cuts cuts;
   try {
      cuts = dataSubselector::Cuts(file, "EVENTS");
   } catch (tip::TipException & eObj) { 
      if (st_facilities::Util::expectedException(eObj, "Could not open")) {
// Probably have a counts map file:
         cuts = dataSubselector::Cuts(file, "", false);
      } else {
         throw;
      }
   }
   for (unsigned int i = 0; i < cuts.size(); i++) {
      if (cuts[i].type() == "range") {
         const dataSubselector::RangeCut & myCut = 
            dynamic_cast<dataSubselector::RangeCut &>(
               const_cast<dataSubselector::CutBase &>(cuts[i]));
         if (myCut.colname() == "EVENT_CLASS") {
            if (myCut.intervalType() == dataSubselector::RangeCut::MAXONLY
                && myCut.maxVal() == 1) {
               return "DC2_A";
            } else {
               return respcombos[static_cast<int>(myCut.minVal())+1];
            }
         }
      }
   }
// Default action: use all event classes
   return respBase;
}

void AppHelpers::readScData() {
   st_app::AppParGroup & pars(*m_pars);
   std::string scFile = pars["scfile"];
   std::string sctable = pars["sctable"];
   st_facilities::Util::file_ok(scFile);
   st_facilities::Util::resolve_fits_files(scFile, m_scFiles);
   std::vector<std::string>::const_iterator scIt = m_scFiles.begin();
   for ( ; scIt != m_scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      m_scData->readData(*scIt, false, sctable);
   }
}

void AppHelpers::readExposureMap() {
   st_app::AppParGroup & pars(*m_pars);
   std::string exposureFile = pars["exposure_map_file"];
   if (exposureFile != "none") {
      st_facilities::Util::file_ok(exposureFile);
      m_expMap->readExposureFile(exposureFile);
   }
}

void AppHelpers::createResponseFuncs(const std::string & analysisType) {
   m_respFuncs = new ResponseFunctions();
   st_app::AppParGroup & pars(*m_pars);
   std::string respBase = pars["rspfunc"];
   std::string evfile;
   if (analysisType == "UNBINNED") {
      std::string myfile = pars["evfile"];
      evfile = myfile;
   } else if (analysisType == "BINNED") {
      std::string myfile = pars["counts_map_file"];
      evfile = myfile;
   } else {
      m_respFuncs->load(respBase);
      return;
   }
   std::vector<std::string> files;
   st_facilities::Util::resolve_fits_files(evfile, files);
   if (respBase == "DSS") {
      std::string respFuncs = responseFuncs(files.front(), "DC2");
      m_respFuncs->load(respFuncs, "DC2");
   } else {
      m_respFuncs->load(respBase);
   }      
}

void AppHelpers::checkOutputFile(bool clobber, const std::string & file) {
   if (!clobber) {
      if (file != "none" && st_facilities::Util::fileExists(file)) {
         st_stream::StreamFormatter formatter("AppHelpers", 
                                              "checkOutputFile", 2);
         formatter.info() << "Output file " << file << " already exists,\n" 
                          << "and you have set 'clobber' to 'no'.\n"
                          << "Please provide a different output file name."
                          << std::endl;
         std::exit(1);
      }
   }
}

void AppHelpers::checkCuts(const std::string & file1,
                           const std::string & ext1,
                           const std::string & file2,
                           const std::string & ext2,
                           bool compareGtis,
                           bool relyOnStreams,
                           bool skipEventClassCuts,
                           bool gtiWarningOnly) {
   bool checkColumns(false);
   bool skipTimeRangeCuts(false);
   dataSubselector::Cuts cuts1(file1, ext1, checkColumns,
                               skipTimeRangeCuts, skipEventClassCuts);
   dataSubselector::Cuts cuts2(file2, ext2, checkColumns,
                               skipTimeRangeCuts, skipEventClassCuts);

   dataSubselector::Cuts gtiCuts1 = gtiCuts(cuts1);
   dataSubselector::Cuts gtiCuts2 = gtiCuts(cuts2);

   if (!checkCuts(cuts1, cuts2, false, relyOnStreams)) {
      std::ostringstream message;
      message << "AppHelpers::checkCuts:\n" 
              << "DSS keywords ";
      message << "in " << file1;
      if (ext1 != "") {
         message << "[" << ext1 << "] ";
      }
      message << "do not match those in " << file2;
      if (ext2 != "") {
         message << "[" << ext2 << "] ";
      }
      throw std::runtime_error(message.str());
   }
   if (compareGtis && !checkCuts(cuts1, cuts2, true, relyOnStreams)) {
      std::ostringstream message;
      message << "AppHelpers::checkCuts:\n"
              << "GTIs in " << file1;
      if (ext1 != "") {
         message << "[" << ext1 << "] ";
      }
      message << "do not match those in " << file2;
      if (ext2 != "") {
         message << "[" << ext2 << "] ";
      }
      if (gtiWarningOnly) {
         st_stream::StreamFormatter formatter("AppHelpers", "checkCuts", 2);
         formatter.warn() << "\nWARNING: \n" << message.str() << "\n\n";
      } else {
            throw std::runtime_error(message.str());
      }
   }
}

void AppHelpers::checkCuts(const std::vector<std::string> & files1,
                           const std::string & ext1,
                           const std::string & file2,
                           const std::string & ext2,
                           bool compareGtis, 
                           bool relyOnStreams, 
                           bool skipEventClassCuts,
                           bool gtiWarningOnly) {
   bool checkColumns(false);
   bool skipTimeRangeCuts(false);
   dataSubselector::Cuts cuts1(files1, ext1, checkColumns,
                               skipTimeRangeCuts, skipEventClassCuts);
   dataSubselector::Cuts cuts2(file2, ext2, checkColumns,
                               skipTimeRangeCuts, skipEventClassCuts);

   dataSubselector::Cuts gtiCuts1 = gtiCuts(cuts1);
   dataSubselector::Cuts gtiCuts2 = gtiCuts(cuts2);

   if (!checkCuts(cuts1, cuts2, false, relyOnStreams)) {
      std::ostringstream message;
      message << "AppHelpers::checkCuts:\n" 
              << "DSS keywords ";
      message << "in \n";
      for (unsigned int i = 0; i < files1.size(); i++) {
         message << files1.at(i) << "\n";
      }
      if (ext1 != "") {
         message << "in extension " << ext1 << "\n";
      }
      message << "do not match those in " << file2;
      if (ext2 != "") {
         message << "[" << ext2 << "] ";
      }
      throw std::runtime_error(message.str());
   }
   if (compareGtis && !checkCuts(cuts1, cuts2, true, relyOnStreams)) {
      std::ostringstream message;
      message << "AppHelpers::checkCuts:\n" 
              << "GTIs in\n";
      for (unsigned int i = 0; i < files1.size(); i++) {
         message << files1.at(i) << "\n";
      }
      if (ext1 != "") {
         message << "in extension " << ext1 << "\n";
      }
      message << "do not match those in " << file2;
      if (ext2 != "") {
         message << "[" << ext2 << "] ";
      }
      if (gtiWarningOnly) {
         st_stream::StreamFormatter formatter("AppHelpers", "checkCuts", 2);
         formatter.warn() << "\nWARNING: \n" << message.str() << "\n\n";
      } else {
         throw std::runtime_error(message.str());
      }
   }
}

bool AppHelpers::checkCuts(const dataSubselector::Cuts & cuts1,
                           const dataSubselector::Cuts & cuts2,
                           bool compareGtis, bool relyOnStreams) {
   bool standardTest;
   if (relyOnStreams) {
      std::ostringstream c1, c2;
      cuts1.writeCuts(c1);
      cuts2.writeCuts(c2);
      standardTest = (c1.str() == c2.str());
   } else {
      if (compareGtis) {
         standardTest = (cuts1 == cuts2);
      } else {
         standardTest = cuts1.compareWithoutGtis(cuts2);
      }
   }
   return standardTest;
}

void AppHelpers::checkTimeCuts(const std::string & file1,
                               const std::string & ext1,
                               const std::string & file2,
                               const std::string & ext2,
                               bool compareGtis,
                               bool gtiWarningOnly) {
   dataSubselector::Cuts cuts1(file1, ext1, false, false, true);
   dataSubselector::Cuts cuts2(file2, ext2, false, false, true);
   if (!checkTimeCuts(cuts1, cuts2, compareGtis)) {
      std::ostringstream message;
      message << "AppHelpers::checkTimeCuts:\n" 
              << "Time range cuts ";
      if (compareGtis) {
         message << "and GTI extensions ";
      }
      message << "in files " << file1;
      if (ext1 != "") {
         message << "[" << ext1 << "]";
      }
      message << "and " << file2;
      if (ext2 != "") {
         message << "[" << ext2 << "]";
      }
      message << "do not agree.";
      if (gtiWarningOnly) {
         st_stream::StreamFormatter formatter("AppHelpers", 
                                              "checkTimeCuts", 2);
         formatter.warn() << "\nWARNING: \n" << message.str() << "\n\n";
      } else {
         throw std::runtime_error(message.str());
      }
   }
}

void AppHelpers::checkTimeCuts(const std::vector<std::string> & files1,
                               const std::string & ext1,
                               const std::string & file2,
                               const std::string & ext2,
                               bool compareGtis,
                               bool gtiWarningOnly) {
   dataSubselector::Cuts cuts1(files1, ext1, false, false, true);
   dataSubselector::Cuts cuts2(file2, ext2, false, false, true);
   if (!checkTimeCuts(cuts1, cuts2, compareGtis)) {
      std::ostringstream message;
      message << "AppHelpers::checkTimeCuts:\n" 
              << "Time range cuts ";
      if (compareGtis) {
         message << "and GTI extensions ";
      }
      message << "in files \n";
      for (unsigned int i = 0; i < files1.size(); i++) {
         message << files1.at(i);
         if (ext1 != "") {
            message << "[" << ext1 << "]";
         }
      }
      message << "\n";
      message << "and " << file2;
      if (ext2 != "") {
         message << "[" << ext2 << "]\n";
      }
      message << "do not agree.";
      if (gtiWarningOnly) {
         st_stream::StreamFormatter formatter("AppHelpers", 
                                              "checkTimeCuts", 2);
         formatter.warn() << "\nWARNING: \n" << message.str() << "\n"
                          << std::endl;
      } else {
         throw std::runtime_error(message.str());
      }
   }
}

bool AppHelpers::checkTimeCuts(const dataSubselector::Cuts & cuts1,
                               const dataSubselector::Cuts & cuts2,
                               bool compareGtis) {
// Assume GTIs encapsulate all of the time range information, so that
// individual time range cuts need not be checked.

   (void)(compareGtis);
   std::vector<const dataSubselector::CutBase *> time_cuts1;
   std::vector<const dataSubselector::CutBase *> time_cuts2;
   gatherGtiCuts(cuts1, time_cuts1);
   gatherGtiCuts(cuts2, time_cuts2);
   bool ok(true);
   if (time_cuts1.size() == time_cuts2.size()) {
      for (unsigned int i = 0; i < time_cuts1.size(); i++) {
         ok = ok && *(time_cuts1[i]) == *(time_cuts2[i]);
      }
   } else {
      ok = false;
   }
   return ok;
}

void AppHelpers::
checkExpMapCuts(const std::vector<std::string> & evFiles,
                const std::string & expMap,
                const std::string & evfileExt,
                const std::string & expMapExt) {
   dataSubselector::Cuts evCuts(evFiles, evfileExt, false);
   dataSubselector::Cuts expMapCuts(expMap, expMapExt, false);

   std::vector<dataSubselector::RangeCut *> evEnergyCuts;
   std::vector<dataSubselector::RangeCut *> expMapEnergyCuts;

   evCuts.removeRangeCuts("ENERGY", evEnergyCuts);
   expMapCuts.removeRangeCuts("ENERGY", expMapEnergyCuts);

   if (evCuts != expMapCuts) {
      std::ostringstream message;
      message << "AppHelpers::checkExpMapCuts:\n"
              << "Inconsistent DSS keywords in "
              << "event file(s) and unbinned exposure map file."
              << std::endl;
      throw std::runtime_error(message.str());
   }
   double ev_emin, ev_emax;
   ::getRangeBounds(evEnergyCuts, ev_emin, ev_emax);

   double expMap_emin, expMap_emax;
   ::getRangeBounds(expMapEnergyCuts, expMap_emin, expMap_emax);

   if (expMap_emin > ev_emin || expMap_emax < ev_emax) {
      std::ostringstream message;
      message << "AppHelpers::checkExpMapCuts:\n"
              << "Energy ranges for event file(s) and "
              << "unbinned exposure map are not consistent.";
      throw std::runtime_error(message.str());
   }
}

void AppHelpers::
gatherGtiCuts(const dataSubselector::Cuts & cuts,
              std::vector<const dataSubselector::CutBase *> & gti_cuts) {
   for (unsigned int i = 0; i < cuts.size(); i++) {
      if (cuts[i].type() == "GTI") {
         gti_cuts.push_back(&cuts[i]);
      }
   }
}

dataSubselector::Cuts 
AppHelpers::gtiCuts(const dataSubselector::Cuts & cuts) {
   std::vector<const dataSubselector::CutBase *> gti_cuts;
   gatherGtiCuts(cuts, gti_cuts);
   
   dataSubselector::Cuts my_gtiCuts;
   for (size_t i = 0; i < gti_cuts.size(); i++) {
      my_gtiCuts.addCut(*gti_cuts.at(i));
   }
   return my_gtiCuts;
}

} // namespace Likelihood
