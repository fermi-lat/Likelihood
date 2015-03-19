/**
 * @file AppHelpers.cxx
 * @brief Class of "helper" methods for Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/AppHelpers.cxx,v 1.121 2015/02/28 03:27:31 jchiang Exp $
 */

#include <cmath>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "st_stream/StreamFormatter.h"

#include "astro/SkyDir.h"

#include "st_facilities/Util.h"

#include "evtbin/Gti.h"

#include "dataSubselector/BitMaskCut.h"
#include "dataSubselector/CutBase.h"
#include "dataSubselector/Cuts.h"

#include "optimizers/Gaussian.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BandFunction.h"
#include "Likelihood/BinnedExposure.h"
#include "Likelihood/BrokenPowerLaw2.h"
#include "Likelihood/BrokenPowerLaw3.h"
#include "Likelihood/BrokenPowerLawExpCutoff.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/DMFitFunction.h"
#include "Likelihood/EblAtten.h"
#include "Likelihood/EnergyBand.h"
#include "Likelihood/EventContainer.h"
#include "Likelihood/ExpCutoff.h"
#include "Likelihood/ExpCutoffSEDPeak.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/FileFunction.h"
#include "Likelihood/LogGaussian.h"
#include "Likelihood/LogNormal.h"
#include "Likelihood/LogNormalLog.h"
#include "Likelihood/LogParabola.h"
#include "Likelihood/MapCubeFunction2.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/MultipleBrokenPowerLaw.h"
#include "Likelihood/Observation.h"
#include "Likelihood/PiecewisePowerLaw.h"
#include "Likelihood/PowerLawSuperExpCutoff.h"
#include "Likelihood/PowerLaw2.h"
#include "Likelihood/RadialProfile.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/ScaleFactor.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/SmoothBrokenPowerLaw.h"
#include "Likelihood/SmoothDoubleBrokenPowerLaw.h"
#include "Likelihood/SpatialMap.h"
#include "Likelihood/WcsMap2.h"

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
   void strip_at_sign(std::string & input) {
      if (input.find_first_of("@") == 0) {
         std::string output = "";
         std::string::iterator it = input.begin() + 1;
         for ( ; it != input.end(); ++it) {
            output += *it;
         }
         input = output;
      }
   }
}

namespace Likelihood {

AppHelpers::AppHelpers(st_app::AppParGroup * pars,
                       const std::string & analysisType) 
   : m_pars(pars), m_funcFactory(0), m_scData(0), m_expCube(0),
     m_expMap(0), m_respFuncs(0), m_roiCuts(0), m_eventCont(0),
     m_bexpmap(0), m_phased_expmap(0), m_meanpsf(0), m_irfsName(""),
     m_respFuncCuts(0) {
   st_app::AppParGroup & my_pars(*m_pars);
   prepareFunctionFactory();
   createResponseFuncs(analysisType);

   m_roiCuts = new RoiCuts();
   m_scData = new ScData();
   m_expCube = new ExposureCube();
   m_expCube->setEfficiencyFactor(m_respFuncs->efficiencyFactor());
   try {
      std::string expcube = my_pars["expcube"];
      if (expcube != "" && expcube != "none") {
         m_expCube->readExposureCube(expcube);
      }
   } catch (hoops::Hexception &) {
   }
   m_expMap = new ExposureMap();
   m_eventCont = new EventContainer(*m_respFuncs, *m_roiCuts, *m_scData);
   if (analysisType == "BINNED") {
      try {
         std::string bexpmap = my_pars["bexpmap"];
         m_bexpmap = new BinnedExposure(bexpmap);
      } catch (hoops::Hexception &) {
      }
      try {
         std::string phased_expmap = my_pars["phased_expmap"];
         if (phased_expmap != "none" && phased_expmap != "") {
            m_phased_expmap = new WcsMap2(phased_expmap);
         }
      } catch (hoops::Hexception &) {
      }
   }
   m_observation = new Observation(m_respFuncs,
                                   m_scData,
                                   m_roiCuts,
                                   m_expCube,
                                   m_expMap,
                                   m_eventCont,
                                   m_bexpmap,
                                   m_phased_expmap);
   if (analysisType == "BINNED") {
      try {
         // If there is no livetime cube, skip the setting of the mean psf.
         my_pars["expcube"];
         std::string cmapfile;
         // gtmodel still uses "srcmaps" as the parameter name of the
         // counts cube that is used to determine the map geometry,
         // so need to specialize here.
         try {
            std::string tmp = my_pars["cmap"];
            cmapfile = tmp;
         } catch (hoops::Hexception & ee) {
            std::string tmp = my_pars["srcmaps"];
            cmapfile = tmp;
         }
         CountsMap cmap(cmapfile);
         m_meanpsf = new MeanPsf(cmap.refDir(), cmap.energies(),
                                 *m_observation);
         m_observation->setMeanPsf(m_meanpsf);
      } catch (hoops::Hexception &) {
         // Do nothing.
      }
   }
}

AppHelpers::~AppHelpers() {
   delete m_funcFactory;
   delete m_roiCuts;
   delete m_scData;
   delete m_expCube;
   delete m_expMap;
   delete m_eventCont;
   delete m_respFuncs;
   delete m_bexpmap;
   delete m_phased_expmap;
   delete m_meanpsf;
   delete m_observation;
   delete m_respFuncCuts;
}

optimizers::FunctionFactory & AppHelpers::funcFactory() {
   return *m_funcFactory;
}

void AppHelpers::prepareFunctionFactory() {
   m_funcFactory = new optimizers::FunctionFactory;
   addFunctionPrototypes(m_funcFactory);
}

void AppHelpers::
addFunctionPrototypes(optimizers::FunctionFactory * funcFactory) {
   bool makeClone(false);
   funcFactory->addFunc(new SkyDirFunction(), makeClone);
   funcFactory->addFunc(new SpatialMap(), makeClone);
   funcFactory->addFunc(new BandFunction(), makeClone);
   funcFactory->addFunc(new LogParabola(), makeClone);
   funcFactory->addFunc(new LogGaussian(), makeClone);
   funcFactory->addFunc(new LogNormal(), makeClone);
   funcFactory->addFunc(new LogNormalLog(), makeClone);
   funcFactory->addFunc(new MapCubeFunction2(), makeClone);
   funcFactory->addFunc(new RadialProfile(), makeClone);
   funcFactory->addFunc(new PowerLaw2(), makeClone);
   funcFactory->addFunc(new BrokenPowerLaw2(), makeClone);
   funcFactory->addFunc(new BrokenPowerLaw3(), makeClone);
   funcFactory->addFunc(new MultipleBrokenPowerLaw(), makeClone);
   funcFactory->addFunc(new PiecewisePowerLaw(), makeClone);
   funcFactory->addFunc(new SmoothBrokenPowerLaw(), makeClone);
   funcFactory->addFunc(new SmoothDoubleBrokenPowerLaw(), makeClone);
   funcFactory->addFunc(new FileFunction(), makeClone);
   funcFactory->addFunc(new ExpCutoff(), makeClone);
   funcFactory->addFunc(new ExpCutoffSEDPeak(), makeClone);
   funcFactory->addFunc(new BrokenPowerLawExpCutoff(), makeClone);
   funcFactory->addFunc(new PowerLawSuperExpCutoff(), makeClone);
   funcFactory->addFunc(new DMFitFunction(), makeClone);

   funcFactory->addFunc(new EblAtten(), makeClone);
   funcFactory->addFunc(new EblAtten(BrokenPowerLaw2()), makeClone);
   funcFactory->addFunc(new EblAtten(LogParabola()), makeClone);
   funcFactory->addFunc(new EblAtten(BandFunction()), makeClone);
   funcFactory->addFunc(new EblAtten(SmoothBrokenPowerLaw()), makeClone);
   funcFactory->addFunc(new EblAtten(FileFunction()), makeClone);
   funcFactory->addFunc(new EblAtten(ExpCutoff()), makeClone);
   funcFactory->addFunc(new EblAtten(BrokenPowerLawExpCutoff()), makeClone);
   funcFactory->addFunc(new EblAtten(PowerLawSuperExpCutoff()), makeClone);

   funcFactory->addFunc(new EnergyBand(), makeClone);
   funcFactory->addFunc(new EnergyBand(BrokenPowerLaw2()), makeClone);
   funcFactory->addFunc(new EnergyBand(LogParabola()), makeClone);
   funcFactory->addFunc(new EnergyBand(BandFunction()), makeClone);
   funcFactory->addFunc(new EnergyBand(SmoothBrokenPowerLaw()), makeClone);
   funcFactory->addFunc(new EnergyBand(FileFunction()), makeClone);
   funcFactory->addFunc(new EnergyBand(ExpCutoff()), makeClone);
   funcFactory->addFunc(new EnergyBand(BrokenPowerLawExpCutoff()), makeClone);
   funcFactory->addFunc(new EnergyBand(PowerLawSuperExpCutoff()), makeClone);

   funcFactory->addFunc(new ScaleFactor(FileFunction()), makeClone);
   funcFactory->addFunc(new ScaleFactor(PowerLaw2()), makeClone);
   funcFactory->addFunc(new ScaleFactor(PowerLawSuperExpCutoff()), makeClone);
   funcFactory->addFunc(new ScaleFactor(optimizers::Gaussian()), makeClone);
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
   static const char * respcombos[] = {"DC2", "DC2FA", "DC2BA", 
                                       "DC2FB", "DC2BB", "DC2_A"};
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
   double tmin(m_observation->roiCuts().minTime());
   double tmax(m_observation->roiCuts().maxTime());
   st_app::AppParGroup & pars(*m_pars);
   std::string scFile = pars["scfile"];
   std::string sctable = pars["sctable"];
   ::strip_at_sign(scFile);
   st_facilities::Util::file_ok(scFile);
   st_facilities::Util::resolve_fits_files(scFile, m_scFiles);
   m_scData->readData(m_scFiles, tmin, tmax);
}

void AppHelpers::readExposureMap() {
   st_app::AppParGroup & pars(*m_pars);
   std::string exposureFile = pars["expmap"];
   if (exposureFile != "none") {
      st_facilities::Util::file_ok(exposureFile);
      m_expMap->readExposureFile(exposureFile);
   }
}

void AppHelpers::createResponseFuncs(const std::string & analysisType) {
   // Infer the response functions to use as specified by event_class
   // and event_type selections in the input ft1, cmap, srcmaps, or
   // bexpmap files or from the gt-tool pars.  Do this via a
   // dataSubselector::Cuts object and its CALDB_implied_irfs() member
   // function.
   st_stream::StreamFormatter formatter("AppHelpers", 
                                        "createResponseFuncs", 2);
   m_respFuncs = new ResponseFunctions();
   st_app::AppParGroup & pars(*m_pars);
   std::string respBase = pars["irfs"];
   std::string evfile;
   std::string extname;
   dataSubselector::Cuts * my_cuts(0);
   if (analysisType == "UNBINNED") {
      try {
         std::string my_file = pars["expmap"];
         evfile = my_file;
         extname = "PRIMARY";
      } catch (hoops::Hexception & eObj) {
         evfile = "none";
      }

      // Need to explicitly test for evfile=="none", since
      // that is a possible value for expmap.
      if (evfile == "none") {
         std::string my_file = pars["evfile"];
         evfile = my_file;
         extname = "EVENTS";
      }
   } else if (analysisType == "BINNED") {
      try {
         std::string my_file = pars["bexpmap"];
         evfile = my_file;
      } catch (hoops::Hexception &) {
         // Must be running gtexpcube2.
         std::string my_file = pars["cmap"];
         evfile = my_file;
      }
      extname = "";
   }
   if (respBase == "CALDB") {
      std::vector<std::string> files;
      st_facilities::Util::resolve_fits_files(evfile, files);
      my_cuts = new dataSubselector::Cuts(files.at(0), extname,
                                          false, true, true);
      dataSubselector::BitMaskCut * evclass_cut(0);
      if (!(evclass_cut = my_cuts->bitMaskCut("EVENT_CLASS"))) {
         // The binned expposure map may not have DSS keywords
         // specifying an event class cut. If so, look in the cmap
         // file or the srcmaps file.
         delete my_cuts;
         try {
            std::string my_file = pars["cmap"];
            evfile = my_file;
         } catch (hoops::Hexception &) {
            // Running gtmodel
            std::string my_file = pars["srcmaps"];
            evfile = my_file;
         }
         my_cuts = new dataSubselector::Cuts(evfile, "", false, true, true);
      }
      // If the pars["evtype"] option is available and specified (i.e., not
      // INDEF) and if my_cuts does not already have one, add an EVENT_TYPE
      // bit-mask cut with the provided value.
      dataSubselector::BitMaskCut * evtype_cut(0);
      if (!(evtype_cut = my_cuts->bitMaskCut("EVENT_TYPE"))) {
         try {
            unsigned int evtype = pars["evtype"];
            my_cuts->addBitMaskCut("EVENT_TYPE", evtype, my_cuts->pass_ver());
         } catch (hoops::Hexception &) {
            /// Do nothing.
         }
      } else {
         delete evtype_cut;
      }
   } else {
      // User has specified a value for irfs other than the default of
      // "CALDB".  In this case, generate a my_cuts object and add
      // event_class and event_type BitMaskCuts based on the values in
      // the "evtype" parameter.
      //
      // Disallow event_type qualifiers on irfs name.
      std::string::size_type pos(respBase.find_first_of(":"));
      if (pos != std::string::npos) {
         throw std::runtime_error("event_type qualifiers should not be "
                                  "included in the irfs specification.\n"
                                  "event_type information should be provided "
                                  "via DSS keywords or an evtype command line "
                                  "option, if available.");
      }
      my_cuts = new dataSubselector::Cuts();

      // The member function dataSubselector::Cuts::setIrfs adds the
      // event_class BitMaskCut and sets the pass version based on the
      // irf name.
      my_cuts->setIrfs(respBase);
      
      // Add event_type BitMaskCut.
      try {
         unsigned int evtype = pars["evtype"];
         if (evtype == 3) {
            formatter.info() << "Using evtype=3 (i.e., FRONT/BACK irfs)\n";
         }
         my_cuts->addBitMaskCut("EVENT_TYPE", evtype, my_cuts->pass_ver());
      } catch (hoops::Hexception &) {
         // evtype not a command line option, so do not add this cut.
      }
   }
   std::vector<unsigned int> selectedEvtTypes;
   if (my_cuts->pass_ver() != "NONE") {
      // We have Pass 7 or later, so reset the response function name,
      // identifying the event_type partition to use.
      respBase = my_cuts->CALDB_implied_irfs();
      getSelectedEvtTypes(*my_cuts, selectedEvtTypes);
   }

   m_respFuncs->load(respBase, "", selectedEvtTypes);

   m_irfsName = respBase;

//   delete my_cuts;
   m_respFuncCuts = my_cuts;
}

void AppHelpers::
getSelectedEvtTypes(const std::string & evfile,
                    const std::string & extname,
                    std::vector<unsigned int> & selectedEvtTypes,
                    unsigned int evtype_bit_mask) {
   dataSubselector::Cuts my_cuts(evfile, extname, false);
   getSelectedEvtTypes(my_cuts, selectedEvtTypes, evtype_bit_mask);
}

void AppHelpers::
getSelectedEvtTypes(const dataSubselector::Cuts & cuts,
                    std::vector<unsigned int> & selectedEvtTypes,
                    unsigned int evtype_bit_mask) {
   selectedEvtTypes.clear();

   // Find the event_type BitMaskCut, if it exists, and get mask.
   for (size_t i(0); i < cuts.size(); i++) {
      if (cuts[i].type() == "bit_mask") {
         const dataSubselector::BitMaskCut & event_type_cut
            = dynamic_cast<const dataSubselector::BitMaskCut &>(cuts[i]);
         if (event_type_cut.colname() == "EVENT_TYPE") {
            evtype_bit_mask = event_type_cut.mask();
         }
      }
   }
   // Test each bit in the mask and save non-zero positions, which
   // correspond to the selected event_types.
   for (size_t j(0); j < 32; j++) {
      if ((evtype_bit_mask & (1 << j)) != 0) {
         selectedEvtTypes.push_back(j);
      }
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
   double gti_start_diffs, gti_stop_diffs;
   if (!checkTimeCuts(cuts1, cuts2, compareGtis, gti_start_diffs,
                      gti_stop_diffs)) {
      std::ostringstream message;
      message << "AppHelpers::checkTimeCuts:\n" 
              << "Time range cuts ";
      if (compareGtis) {
         message << "and GTI extensions ";
      }
      message << "in files " << file1;
      if (ext1 != "") {
         message << "[" << ext1 << "] ";
      }
      message << "and " << file2;
      if (ext2 != "") {
         message << "[" << ext2 << "] ";
      }
      message << "do not agree." << std::endl;
      message << "Aggregate absolute differences in GTI start times (s): "
              << gti_start_diffs << std::endl;
      message << "Aggregate absolute differences in GTI stop times (s): "
              << gti_stop_diffs << std::endl;
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
   double gti_start_diffs, gti_stop_diffs;
   if (!checkTimeCuts(cuts1, cuts2, compareGtis, gti_start_diffs,
                      gti_stop_diffs)) {
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
      message << "do not agree.\n";
      message << "Aggregate absolute differences in GTI start times (s): "
              << gti_start_diffs << std::endl;
      message << "Aggregate absolute differences in GTI stop times (s): "
              << gti_stop_diffs << std::endl;
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
                               bool compareGtis,
                               double & gti_start_diffs,
                               double & gti_stop_diffs) {
// Assume GTIs encapsulate all of the time range information, so that
// individual time range cuts need not be checked.

   (void)(compareGtis);
   std::vector<const dataSubselector::CutBase *> time_cuts1;
   std::vector<const dataSubselector::CutBase *> time_cuts2;
   gatherGtiCuts(cuts1, time_cuts1);
   gatherGtiCuts(cuts2, time_cuts2);
   bool ok(true);
   if (time_cuts1.size() == time_cuts2.size()) {
      gti_start_diffs = 0;
      gti_stop_diffs = 0;
      for (unsigned int i = 0; i < time_cuts1.size(); i++) {
         ok = ok && *(time_cuts1[i]) == *(time_cuts2[i]);
         if (!(*(time_cuts1[i]) == *(time_cuts2[i]))) {
            std::cout << "difference in gti interval found" << std::endl;
            const dataSubselector::GtiCut * gticut1 
               = dynamic_cast<const dataSubselector::GtiCut *>(time_cuts1[i]);
            const dataSubselector::GtiCut * gticut2
               = dynamic_cast<const dataSubselector::GtiCut *>(time_cuts2[i]);
            evtbin::Gti::ConstIterator it1(gticut1->gti().begin());
            evtbin::Gti::ConstIterator it2(gticut2->gti().begin());
            for ( ; it1 != gticut1->gti().end() && it2 != gticut2->gti().end();
                  ++it1, ++it2) {
               gti_start_diffs += std::fabs(it1->first - it2->first);
               gti_stop_diffs += std::fabs(it1->second - it2->second);
            }
         }
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

   evCuts.removeVersionCut("IRF_VERSION");
   expMapCuts.removeVersionCut("IRF_VERSION");

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

void AppHelpers::
checkExposureMap(const std::string & cmapfile,
                 const std::string & emapfile) {
   CountsMap cmap(cmapfile);
   BinnedExposure emap(emapfile);
   emap.setBoundaryFlag(true);
   
   double energy(emap.energies()[0]);

   std::vector<size_t> ii, jj;

   size_t i;
   size_t j(1);
   for (i=1; i < cmap.naxis1() + 1; i++) {
      ii.push_back(i);
      jj.push_back(j);
   }

   i = cmap.naxis1();
   for (; j < cmap.naxis2() + 1; j++) {
      ii.push_back(i);
      jj.push_back(j);
   }

   j = cmap.naxis2();
   for (; i > 0; i--) {
      ii.push_back(i);
      jj.push_back(j);
   }
   
   i = 1;
   for (; j > 0; j--) {
      ii.push_back(i);
      jj.push_back(j);
   }
   for (size_t indx(0); indx < ii.size(); indx++) {
      astro::SkyDir my_dir;
      st_facilities::Util::pixel2SkyDir(cmap.projection(), 
                                        ii[indx], jj[indx], my_dir);
      try {
         emap(energy, my_dir.ra(), my_dir.dec());
      } catch (std::runtime_error & eObj) {
         if (st_facilities::Util::
             expectedException(eObj, "outside of the map boundaries")) {
            throw std::runtime_error("Counts map not covered by exposure map.");
         } else {
            throw;
         }
      }
   }
}

void AppHelpers::setBitMaskCuts(dataSubselector::Cuts & other_cuts) const {
   other_cuts.setBitMaskCut(m_respFuncCuts->bitMaskCut("EVENT_CLASS"));
   other_cuts.setBitMaskCut(m_respFuncCuts->bitMaskCut("EVENT_TYPE"));
}

} // namespace Likelihood
