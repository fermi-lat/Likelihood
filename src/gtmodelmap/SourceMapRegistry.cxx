/**
 * @file SourceMapRegistry.cxx
 * @brief Implementation for a class that manages SourceMaps for
 * gtmodelmap.  This class attempts to be lightweight by using lazy
 * evaluation of the SourceMaps and not storing the maps internally
 * in BinnedLikelhood.
 * @author J. Chiang
 *
 * $Header$
 */

#include "st_facilities/Util.h"

#include "optimizers/FunctionFactory.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/EventContainer.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/SourceMap.h"

#include "SourceMapRegistry.h"

SourceMapRegistry::SourceMapRegistry(const std::string & countsMap,
                                     const std::string & xmlFile,
                                     const std::string & irfs,
                                     const std::string & expCubeFile,
                                     const std::string & binnedExpMap,
                                     optimizers::FunctionFactory & funcFactory)
   : m_observation(0), m_countsMap(0), m_logLike(0), m_sourceMap(0) {

   Likelihood::ResponseFunctions * respFuncs = 
      new Likelihood::ResponseFunctions();
   respFuncs->load(irfs);

   Likelihood::ExposureMap * expMap = new Likelihood::ExposureMap();

   Likelihood::ScData * scData = new Likelihood::ScData();

   Likelihood::RoiCuts * roiCuts = new Likelihood::RoiCuts();
   roiCuts->readCuts(countsMap, "", false);

   Likelihood::ExposureCube * expCube = new Likelihood::ExposureCube();
   expCube->readExposureCube(expCubeFile);

   Likelihood::EventContainer * eventCont = 
      new Likelihood::EventContainer(*respFuncs, *roiCuts, *scData);

   m_observation = new Likelihood::Observation(respFuncs, scData, roiCuts,
                                               expCube, expMap, eventCont);

   if (binnedExpMap != "" && binnedExpMap != "none") {
      Likelihood::SourceMap::setBinnedExposure(binnedExpMap);
   }

   m_countsMap = new Likelihood::CountsMap(countsMap);

   bool computePointSources(false);
   bool applyPsfCorrections(false);
   m_logLike = new Likelihood::BinnedLikelihood(*m_countsMap, *m_observation,
                                                countsMap, computePointSources,
                                                applyPsfCorrections);

   m_logLike->readXml(xmlFile, funcFactory, false);
}

SourceMapRegistry::~SourceMapRegistry() {
   delete &(m_observation->respFuncs());
   delete &(m_observation->scData());
   delete &(m_observation->roiCuts());
   delete &(m_observation->expCube());
   delete &(m_observation->expMap());
   delete &(m_observation->eventCont());

   delete m_observation;
   delete m_countsMap;
   delete m_logLike;
}

const std::vector<float> & 
SourceMapRegistry::sourceMap(const std::string & srcName) {
   try {
      return m_logLike->sourceMap(srcName).model();
   } catch (std::runtime_error & eObj) {
      if (!st_facilities::Util::
          expectedException(eObj, "Cannot find source map named")) {
         throw;
      }
   }
// Re-evaluate the source map with each call to save memory.  This is
// ok since gtmodelmap will call this method only once per source.
   delete m_sourceMap;
   m_sourceMap = m_logLike->createSourceMap(srcName);

   return m_sourceMap->model();
}
