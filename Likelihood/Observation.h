/**
 * @file Observation.h
 * @brief Class to contain all of the information associated with a given
 * observation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Observation.h,v 1.5 2005/10/10 21:42:00 jchiang Exp $
 */

#ifndef Likelihood_Observation_h
#define Likelihood_Observation_h

#include "Likelihood/EventContainer.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"

namespace Likelihood {

   class BinnedExposure;
   class WcsMap2;
   class MeanPsf;

/**
 * @class Observation
 * @brief A container class composed of all of the data-related
 * classes associated with a particular observation.
 *
 * The contained classes were all originally implemented as
 * Singletons.  This class encapsulates the observation information
 * comprising those classes, providing a single access point whilst
 * allowing multiple observations to be considered within the same
 * program instance.
 *
 */

class Observation {

public:

   Observation(ResponseFunctions * respFuncs=0, 
               ScData * scData=0,
               RoiCuts * roiCuts=0,
               ExposureCube * expCube=0,
               ExposureMap * expMap=0,
               EventContainer * eventCont=0,
               BinnedExposure * bexpmap=0,
               WcsMap2 * phased_expmap=0) :
      m_respFuncs(respFuncs), m_scData(scData), m_roiCuts(roiCuts),
      m_expCube(expCube), m_expMap(expMap), m_eventCont(eventCont),
      m_bexpmap(bexpmap), m_phased_expmap(phased_expmap),
      m_meanpsf(0) {
   }

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

   const EventContainer & eventCont() const {
      return *m_eventCont;
   }

   EventContainer & eventCont() {
      return *m_eventCont;
   }

   const BinnedExposure & bexpmap() const {
      return *m_bexpmap;
   }

   BinnedExposure & bexpmap() {
      return *m_bexpmap;
   }
   
   const WcsMap2 & phased_expmap() const {
      return *m_phased_expmap;
   }

   WcsMap2 & phased_expmap() {
      return *m_phased_expmap;
   }

   void setMeanPsf(MeanPsf * meanpsf) {
      m_meanpsf = meanpsf;
   }

   const MeanPsf & meanpsf() const {
      return *m_meanpsf;
   }

   MeanPsf & meanpsf() {
      return *m_meanpsf;
   }

private:

   ResponseFunctions * m_respFuncs;
   ScData * m_scData;
   RoiCuts * m_roiCuts;
   ExposureCube * m_expCube;
   ExposureMap * m_expMap;
   EventContainer * m_eventCont;
   BinnedExposure * m_bexpmap;
   WcsMap2 * m_phased_expmap;
   MeanPsf * m_meanpsf;

};

} //namespace Likelihood

#endif // Likelihood_Observation_h
