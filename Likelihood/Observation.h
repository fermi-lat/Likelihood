/**
 * @file Observation.h
 * @brief Class to contain all of the information associated with a given
 * observation.
 * @author J. Chiang
 *
 * $Header$
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
 * $Header$
 */

class Observation {

public:

//   Observation(st_app::AppParGroup & pars);

   Observation() : m_respFuncs(0), m_scData(0), m_roiCuts(0), m_expCube(0),
      m_expMap(0) {}

   Observation(const ResponseFunctions * respFuncs, 
               const ScData * scData,
               const RoiCuts * roiCuts,
               const ExposureCube * expCube,
               const ExposureMap * expMap) :
      m_respFuncs(respFuncs), m_scData(scData), m_roiCuts(roiCuts),
      m_expCube(expCube), m_expMap(expMap) {}

   const ResponseFunctions & respFuncs() const {
      return *m_respFuncs;
   }

   const ScData & scData() const {
      return *m_scData;
   }

   const RoiCuts & roiCuts() const {
      return *m_roiCuts;
   }

   const ExposureCube & expCube() const {
      return *m_expCube;
   }

   const ExposureMap & expMap() const {
      return *m_expMap;
   }

private:

   /// @todo Replace all these pointers with objects and convert all
   /// these classes from Singletons to proper classes.
   const ResponseFunctions * m_respFuncs;
   const ScData * m_scData;
   const RoiCuts * m_roiCuts;
   const ExposureCube * m_expCube;
   const ExposureMap * m_expMap;

};

} //namespace Likelihood

#endif // Likelihood_Observation_h
