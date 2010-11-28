/**
 * @file ExposureCube.cxx
 * @brief Implementation for ExposureCube wrapper class of map_tools::Exposure
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/ExposureCube.cxx,v 1.9 2010/11/28 03:52:08 jchiang Exp $
 */

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "Likelihood/ExposureCube.h"
#include "Likelihood/Observation.h"

namespace Likelihood {

ExposureCube::ExposureCube(const ExposureCube & other) 
   : m_exposure(new map_tools::Exposure(*(other.m_exposure))),
     m_weightedExposure(0), m_efficiencyFactor(0),
     m_haveFile(other.m_haveFile),
     m_fileName(other.m_fileName), 
     m_hasPhiDependence(other.m_hasPhiDependence) {
   if (other.m_weightedExposure) {
      m_weightedExposure = new map_tools::Exposure(*(other.m_weightedExposure));
   }
   if (other.m_efficiencyFactor) {
      m_efficiencyFactor = other.m_efficiencyFactor->clone();
   }
}

void ExposureCube::readExposureCube(std::string filename) {
   facilities::Util::expandEnvVar(&filename);
   m_fileName = filename;
   m_exposure = new map_tools::Exposure(filename);
   try {
      m_weightedExposure = new map_tools::Exposure(filename,
                                                   "WEIGHTED_EXPOSURE");
   } catch(tip::TipException &) {
      m_weightedExposure = 0;
   }
   m_haveFile = true;
   m_hasPhiDependence = phiDependence(filename);
}

bool ExposureCube::phiDependence(const std::string & filename) const {
   const tip::Table * table 
      = tip::IFileSvc::instance().readTable(filename, "EXPOSURE");
   const tip::Header & header(table->getHeader());
   long nphibins;
   try {
      header["PHIBINS"].get(nphibins);
   } catch (tip::TipException &) {
      nphibins = 0;
   }
   return nphibins > 0;
}

ExposureCube::Aeff::Aeff(double energy, int evtType, 
                         const Observation & observation) 
   : m_energy(energy), m_evtType(evtType), m_observation(observation) {
// Turn off phi-dependence if omitted from livetime cube.
   bool phi_dependence(m_observation.expCube().hasPhiDependence());
   if (!phi_dependence) {
      std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
         = m_observation.respFuncs().begin();
      for ( ; respIt != m_observation.respFuncs().end(); ++respIt) {
         respIt->second->aeff()->setPhiDependence(false);
      }
   }
}

double ExposureCube::Aeff::operator()(double cosTheta, double phi) const {
   double inclination = acos(cosTheta)*180./M_PI;
// Check if we need to truncate the inclination past 70 deg as in 
// MeanPsf::Aeff::operator()(...)
   // if (inclination > 70.) {
   //    return 0;
   // }
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = m_observation.respFuncs().begin();
   for ( ; respIt != m_observation.respFuncs().end(); ++respIt) {
      if (respIt->second->irfID() == m_evtType) {
         irfInterface::IAeff * aeff = respIt->second->aeff();
         double aeff_val = aeff->value(m_energy, inclination, phi);
         return aeff_val;
      }
   }
   return 0;
}

} // namespace Likelihood
