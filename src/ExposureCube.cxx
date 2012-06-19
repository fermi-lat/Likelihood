/**
 * @file ExposureCube.cxx
 * @brief Implementation for ExposureCube wrapper class of map_tools::Exposure
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/ExposureCube.cxx,v 1.14 2012/06/19 04:03:50 jchiang Exp $
 */

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "healpix/CosineBinner.h"

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

double ExposureCube::livetime(const astro::SkyDir & dir,
                              double costheta, double phi) const {
   size_t ci(healpix::CosineBinner::cosine_index(costheta));
   const healpix::CosineBinner & binner(m_exposure->data()[dir]);
   if (phi < 0) {
      return binner.at(ci);
   }
   size_t index((healpix::CosineBinner::phi_index(phi*M_PI/180.) + 1)*
                healpix::CosineBinner::nbins() + ci);
   return binner.at(index);
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

double ExposureCube::AeffBase::operator()(double cosTheta, double phi) const {
   std::pair<double, double> key(cosTheta, phi);
   AeffCacheMap_t::const_iterator it(m_cache.find(key));
   if (it != m_cache.end()) {
      return it->second;
   }
   double my_value(value(cosTheta, phi));
   m_cache.insert(std::make_pair(key, my_value));
   return my_value;
}

double ExposureCube::Aeff::value(double cosTheta, double phi) const {
   double inclination = acos(cosTheta)*180./M_PI;
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
