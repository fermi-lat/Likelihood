/**
 * @file Observation.h
 * @brief Class to contain all of the information associated with a given
 * observation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Observation.h,v 1.2 2005/03/03 00:46:51 jchiang Exp $
 */

#ifndef Likelihood_Observation_h
#define Likelihood_Observation_h

#include "Likelihood/ExposureCube.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"

namespace Likelihood {

/**
 * @class Observation
 * @brief A container class composed of all of the classes associated
 * with a particular observation, with the exception of the event data,
 * which are contained in the classes in the SourceModel hierarchy.
 *
 * The contained classes were(are) all originally implemented as
 * Singletons.  This class encapsulates the observation information
 * comprising those classes, providing a single access point whilst
 * allowing multiple observations to be considered within the same
 * program instance.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Observation.h,v 1.2 2005/03/03 00:46:51 jchiang Exp $
 */

class Observation {

public:

   Observation() : m_respFuncs(0), m_scData(0), m_roiCuts(0), m_expCube(0),
      m_expMap(0) {}

   Observation(ResponseFunctions * respFuncs, 
               ScData * scData,
               RoiCuts * roiCuts,
               ExposureCube * expCube,
               ExposureMap * expMap) :
      m_respFuncs(respFuncs), m_scData(scData), m_roiCuts(roiCuts),
      m_expCube(expCube), m_expMap(expMap) {}

   const ResponseFunctions & respFuncs() const {
      return *m_respFuncs;
   }

   ResponseFunctions & respFuncs() {
      return *m_respFuncs;
   }

   const ScData & scData() const {
      return *m_scData;
   }

   ScData & scData() {
      return *m_scData;
   }

   const RoiCuts & roiCuts() const {
      return *m_roiCuts;
   }

   RoiCuts & roiCuts() {
      return *m_roiCuts;
   }

   const ExposureCube & expCube() const {
      return *m_expCube;
   }

   ExposureCube & expCube() {
      return *m_expCube;
   }

   const ExposureMap & expMap() const {
      return *m_expMap;
   }

   ExposureMap & expMap() {
      return *m_expMap;
   }

private:

   /// @todo Replace all these pointers with objects and convert all
   /// these classes from Singletons to proper classes.
   ResponseFunctions * m_respFuncs;
   ScData * m_scData;
   RoiCuts * m_roiCuts;
   ExposureCube * m_expCube;
   ExposureMap * m_expMap;

};

} //namespace Likelihood

#endif // Likelihood_Observation_h
