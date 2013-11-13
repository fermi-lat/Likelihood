/**
 * @file AppHelpers.h
 * @brief Class of "helper" methods for the Likelihood applications.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/AppHelpers.h,v 1.38 2013/08/26 22:55:31 jchiang Exp $
 */

#ifndef Likelihood_AppHelpers
#define Likelihood_AppHelpers

#include <iostream>
#include <string>

#include "st_app/AppParGroup.h"

#include "optimizers/FunctionFactory.h"

#include "dataSubselector/Cuts.h"

namespace dataSubselector {
   class CutBase;
}

namespace Likelihood {

   class BinnedExposure;
   class EventContainer;
   class ExposureCube;
   class ExposureMap;
   class MeanPsf;
   class Observation;
   class ResponseFunctions;
   class RoiCuts;
   class ScData;
   class WcsMap2;

/**
 * @class AppHelpers
 * @brief The methods in this class call various static methods for
 * reading in spacecraft data, defining the region-of-interest, and
 * preparing the response functions --- all standard tasks which must
 * be performed as part of any Likelihood analysis.
 *
 */

class AppHelpers {

public:

   AppHelpers() : m_pars(0), m_funcFactory(0), m_observation(0),
                  m_scData(0), m_expCube(0), m_expMap(0), m_respFuncs(0),
                  m_roiCuts(0), m_eventCont(0), m_bexpmap(0), 
                  m_phased_expmap(0), m_meanpsf(0), m_irfsName("") {
      prepareFunctionFactory();
   }
#ifndef SWIG
   AppHelpers(st_app::AppParGroup * pars, 
              const std::string & analysisType="UNBINNED");

   ~AppHelpers();

   void readScData();
   void readExposureMap();
   void setRoi(const std::string & filename="",
               const std::string & ext="EVENTS",
               bool strict=true);

   void checkOutputFile() {
      st_app::AppParGroup & pars(*m_pars);
      checkOutputFile(pars["clobber"], pars["outfile"]);
   }

   const std::vector<std::string> & scFiles() const {
      return m_scFiles;
   }

   const Observation & observation() const {
      return *m_observation;
   }

   Observation & observation() {
      return *m_observation;
   }

   const std::string & irfsName() const {
      return m_irfsName;
   }
                  
#endif // SWIG

   optimizers::FunctionFactory & funcFactory();

   template<typename T>
   static T param(st_app::AppParGroup & pars, const std::string & parname,
                  const T & def_value);

   static void checkOutputFile(bool clobber, const std::string & filename);

   static void checkCuts(const std::string & file1, const std::string & ext1,
                         const std::string & file2, const std::string & ext2,
                         bool compareGtis=true,
                         bool relyOnStreams=false,
                         bool skipEventClassCuts=false,
                         bool gtiWarningOnly=true);

   static void checkCuts(const std::vector<std::string> & files1,
                         const std::string & ext1,
                         const std::string & file2,
                         const std::string & ext2,
                         bool compareGtis=true,
                         bool relyOnStreams=false,
                         bool skipEventClassCuts=false,
                         bool gtiWarningOnly=true);

   static void checkTimeCuts(const std::string & file1, 
                             const std::string & ext1,
                             const std::string & file2,
                             const std::string & ext2,
                             bool compareGtis=true,
                             bool gtiWarningOnly=true);

   static void checkTimeCuts(const std::vector<std::string> & files1,
                             const std::string & ext1,
                             const std::string & file2,
                             const std::string & ext2,
                             bool compareGtis=true,
                             bool gtiWarningOnly=true);

   static void checkExpMapCuts(const std::vector<std::string> & evfiles,
                               const std::string & expMap,
                               const std::string & evfileExt="EVENTS",
                               const std::string & expMapExt="");

   static std::string responseFuncs(const std::string & file,
                                    const std::string & respBase);

   static void getSelectedEvtTypes(const std::string & evfile,
                                   const std::string & extname,
                                   std::vector<size_t> & selectedEvtTypes);

   static void addFunctionPrototypes(optimizers::FunctionFactory * funcFactory);

   /// Compare the geometry of a binned exposure map to the counts
   /// (our source) map it is intended to serve.  Raise an exception
   /// if the exposure map does not cover the counts map.
   static void checkExposureMap(const std::string & cmapfile,
                                const std::string & emapfile);

protected:

   st_app::AppParGroup * m_pars;
   optimizers::FunctionFactory * m_funcFactory;
   std::vector<std::string> m_scFiles;

   Observation * m_observation;

   ScData * m_scData;
   ExposureCube * m_expCube;
   ExposureMap * m_expMap;
   ResponseFunctions * m_respFuncs;
   RoiCuts * m_roiCuts;
   EventContainer * m_eventCont;

   BinnedExposure * m_bexpmap;
   WcsMap2 * m_phased_expmap;
   MeanPsf * m_meanpsf;

   std::string m_irfsName;

   void prepareFunctionFactory();
   void createResponseFuncs(const std::string & analysisType);

   static bool checkCuts(const dataSubselector::Cuts & cuts1,
                         const dataSubselector::Cuts & cuts2,
                         bool compareGtis, bool relyOnStreams);

   static bool checkTimeCuts(const dataSubselector::Cuts & cuts1,
                             const dataSubselector::Cuts & cuts2,
                             bool compareGtis,
                             double & gti_start_diffs,
                             double & gti_stop_diffs);

   static dataSubselector::Cuts gtiCuts(const dataSubselector::Cuts &);

   static void 
   gatherGtiCuts(const dataSubselector::Cuts & cuts,
                 std::vector<const dataSubselector::CutBase *> & time_cuts);
};

template<typename T>
T AppHelpers::param(st_app::AppParGroup & pars, const std::string & parname,
                    const T & def_value) {
   try {
      T value = pars[parname];
      return value;
   } catch (std::exception & eObj) {
      std::cerr << eObj.what() << "\n"
                << "Using default value of " << def_value
                << " for parameter " <<  parname << std::endl;
      return def_value;
   }
}

} // namespace Likelihood

#endif // Likelihood_AppHelpers
