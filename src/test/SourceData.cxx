/**
 * @file SourceData.cxx
 * @brief Implementation for Likelihood test class data.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/test/SourceData.cxx,v 1.1 2004/02/23 17:52:57 jchiang Exp $
 */

#include "SourceData.h"

SourceData::SourceData() {
   char * srcNames[] = {"Extragalactic Diffuse", "Galactic Diffuse",
                              "PKS 0528+134", "Crab Pulsar", "Geminga"};
   const char * srcTypes[] = {"Diffuse", "Diffuse", "Point", "Point", "Point"};
   const char * spatialModels[] = {"ConstantValue", "SpatialMap", 
                                   "SkyDirFunction", "SkyDirFunction", 
                                   "SkyDirFunction"};
   for (unsigned int i = 0; i < 5; i++) {
      m_srcId[srcNames[i]] = i;
      m_srcTypes.push_back(srcTypes[i]);
      m_spatialModels.push_back(spatialModels[i]);
   }
   setParameters();
}

void SourceData::setParameters() {
// We need to do this by hand to ensure an unbiased test of the source
// data. Here we accumulate parameter data only for the spectral model
// components.

   m_parameters.resize(5);
   unsigned int id = m_srcId["Extragalactic Diffuse"];

   m_parameters[id].push_back(Parameter("Prefactor", 1.32, 1e-5, 1e2, true));
   m_parameters[id].push_back(Parameter("Index", -2.1, -1., -3.5, false));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 50., 2e2, false));

   id = m_srcId["Galactic Diffuse"];
   m_parameters[id].push_back(Parameter("Prefactor", 11., 1e-3, 1e3, true));
   m_parameters[id].push_back(Parameter("Index", -2.1, -1., -3.5, false));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 50., 2e2, false));

   id = m_srcId["PKS 0528+134"];
   m_parameters[id].push_back(Parameter("Prefactor", 13.65, 1e-3, 1e3, true));
   m_parameters[id].push_back(Parameter("Index", -2.46, -1., -3.5, true));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 30., 2e3, false));

   id = m_srcId["Crab Pulsar"];
   m_parameters[id].push_back(Parameter("Prefactor", 27., 1e-3, 1e3, true));
   m_parameters[id].push_back(Parameter("Index", -2.19, -1., -3.5, true));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 30., 2e3, false));

   id = m_srcId["Geminga"];
   m_parameters[id].push_back(Parameter("Prefactor", 23.29, 1e-3, 1e3, true));
   m_parameters[id].push_back(Parameter("Index", -1.66, -1., -3.5, true));
   m_parameters[id].push_back(Parameter("Scale", 1e2, 30., 2e3, false));

   m_paramNames.push_back("Prefactor");
   m_paramNames.push_back("Index");
   m_paramNames.push_back("Scale");
}
