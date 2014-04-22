/** 
 * @file DiffuseSource.cxx
 * @brief DiffuseSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/DiffuseSource.cxx,v 1.62 2013/09/04 05:30:45 jchiang Exp $
 */

#include <algorithm>
#include <string>
#include <valarray>
#include <vector>

#include "optimizers/Function.h"
#include "optimizers/Exception.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/RadialProfile.h"
#include "Likelihood/Observation.h"

namespace Likelihood {

DiffuseSource::DiffuseSource(optimizers::Function * spatialDist,
                             const Observation & observation,
                             bool requireExposure,
                             bool mapBasedIntegral)
   : Source(&observation), 
     m_spatialDist(spatialDist->clone()),
     m_mapBasedIntegral(mapBasedIntegral) {
   m_functions["SpatialDist"] = m_spatialDist;
   m_useEdisp = observation.respFuncs().useEdisp();

// In order to compute exposure, RoiCuts and spacecraft data must be
// available; furthermore, the ExposureMap object must have been
// created.
   if (requireExposure) {
      integrateSpatialDist();
   }
   m_srcType = "Diffuse";
}

void DiffuseSource::integrateSpatialDist() {
   const Observation & obs(*observation());
   m_energies = obs.roiCuts().energies();
   computeExposure(m_energies);
}

void DiffuseSource::computeExposure(const std::vector<double> & energies,
                                    bool verbose) {
   m_energies = energies;
   const Observation & obs(*observation());
   if (obs.expMap().haveMap()) {
      if (m_mapBasedIntegral || ::getenv("MAP_BASED_NPRED")) {
         try {
            // Integrate using the map pixels for the quadrature
            mapBaseObject()->integrateSpatialDist(energies, obs.expMap(),
                                                  m_exposure);
            // Delete internal representation of the map to save memory.
            mapBaseObject()->deleteMap();
         } catch (MapBaseException & eObj) {
            // We do not have a map representation of the sources so
            // integrate using ExposureMap implementation.
            obs.expMap().integrateSpatialDist(energies, m_spatialDist, 
                                              m_exposure);
         }
      } else {
         obs.expMap().integrateSpatialDist(energies, m_spatialDist, 
                                           m_exposure);
         try {
            // Delete internal representation of the map to save memory.
            mapBaseObject()->deleteMap();
         } catch (MapBaseException & eObj) {
            // We do not have a map representation of the sources so
            // do nothing.
         }
      }
   } else {
      std::string what("DiffuseSource: An exposure map must be defined "
                       "if diffuse sources are in the model.");
      throw std::runtime_error(what);
   }
}

DiffuseSource::DiffuseSource(const DiffuseSource & rhs) 
   : Source(rhs),
     m_spatialDist(rhs.m_spatialDist->clone()),
     m_mapBasedIntegral(rhs.m_mapBasedIntegral) {
   m_functions["SpatialDist"] = m_spatialDist;
   m_functions["Spectrum"] = m_spectrum;
}

double DiffuseSource::fluxDensity(const Event & evt,
				  CachedResponse* cResp) const {
   (void)(cResp);
   double my_fluxDensity;
   if (m_useEdisp) {
      const std::vector<double> & trueEnergies = evt.trueEnergies();
      const std::vector<double> & diffuseResponses 
         = evt.diffuseResponse(getName());
      unsigned int npts(trueEnergies.size());
      std::vector<double> my_integrand(npts);
      for (unsigned int k = 0; k < npts; k++) {
         optimizers::dArg energy_arg(trueEnergies[k]);
         my_integrand[k] = (*m_spectrum)(energy_arg)*diffuseResponses[k];
      }
      bool useLog;
      TrapQuad trapQuad(trueEnergies, my_integrand, useLog=true);
      my_fluxDensity = trapQuad.integral();
   } else {
      double trueEnergy = evt.getEnergy(); 
      optimizers::dArg energy_arg(trueEnergy);
      my_fluxDensity = (*m_spectrum)(energy_arg)
         *evt.diffuseResponse(trueEnergy, getName());
   }
   return my_fluxDensity;
}

double DiffuseSource::fluxDensityDeriv(const Event &evt, 
                                       const std::string &paramName,
				       CachedResponse* cResp) const {
   (void)(cResp);
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude".
//
// This method needs to be generalized for spectra that are
// CompositeFunctions.
   double my_fluxDensityDeriv;
   double prefactor;
   if (paramName == "Prefactor" && 
       (prefactor = m_spectrum->getParamValue("Prefactor")) != 0) {
      my_fluxDensityDeriv = fluxDensity(evt)/prefactor;
   } else {
      if (m_useEdisp) {
         const std::vector<double> & trueEnergies = evt.trueEnergies();
         const std::vector<double> & diffuseResponses 
            = evt.diffuseResponse(getName());
         unsigned int npts(trueEnergies.size());
         std::vector<double> my_integrand(npts);
         for (unsigned int k = 0; k < npts; k++) {
            optimizers::dArg energy_arg(trueEnergies[k]);
            my_integrand[k] = m_spectrum->derivByParam(energy_arg, paramName)
               *diffuseResponses[k];
         }
         bool useLog;
         TrapQuad trapQuad(trueEnergies, my_integrand, useLog=true);
         my_fluxDensityDeriv = trapQuad.integral();
      } else {
         double trueEnergy = evt.getEnergy(); 
         optimizers::dArg energy_arg(trueEnergy);
         my_fluxDensityDeriv = m_spectrum->derivByParam(energy_arg, paramName)
            *evt.diffuseResponse(trueEnergy, getName());
      }
   }
   return my_fluxDensityDeriv;
}

const MapBase * DiffuseSource::mapBaseObject() const {
   optimizers::Function * foo = 
      const_cast<optimizers::Function *>(this->spatialDist());
   const MapBase * mapBaseObject = dynamic_cast<MapBase *>(foo);
   if (!mapBaseObject) {
      throw MapBaseException("MapBase object not found for this "
                             + ("diffuse source: " + getName()));
   }
   return mapBaseObject;
}

MapBase * DiffuseSource::mapBaseObject() {
   optimizers::Function * foo = 
      const_cast<optimizers::Function *>(this->spatialDist());
   MapBase * mapBaseObject = dynamic_cast<MapBase *>(foo);
   if (!mapBaseObject) {
      throw MapBaseException("MapBase object not found for this "
                             + ("diffuse source: " + getName()));
   }
   return mapBaseObject;
}

double DiffuseSource::angularIntegral(double energy) const {
   if (spatialDist()->genericName() == "ConstantValue") { 
// Here we have an isotropic source
      return 4*M_PI;
   } else if (spatialDist()->genericName() == "RadialProfile") {
      optimizers::Function * foo = 
         const_cast<optimizers::Function *>(this->spatialDist());
      const RadialProfile * profile = dynamic_cast<RadialProfile *>(foo);
      return profile->angularIntegral();
   }
   return mapBaseObject()->mapIntegral(energy);
}

double DiffuseSource::flux() const {
   return computeEnergyIntegral(*m_spectrum, m_energies);
}

double DiffuseSource::fluxDeriv(const std::string & parName) const {
   FluxDeriv my_functor(*m_spectrum, parName);
   return computeEnergyIntegral(my_functor, m_energies);
}

double DiffuseSource::flux(double emin, double emax, size_t npts) const {
   return computeEnergyIntegral(*m_spectrum, emin, emax, npts);
}

double DiffuseSource::fluxDeriv(const std::string & parName,
                                double emin, double emax, size_t npts) const {
   FluxDeriv my_functor(*m_spectrum, parName);
   return computeEnergyIntegral(my_functor, emin, emax, npts);
}

double DiffuseSource::energyFlux() const {
   EnergyFlux my_functor(*m_spectrum);
   return computeEnergyIntegral(my_functor, m_energies);
}

double DiffuseSource::energyFluxDeriv(const std::string & parName) const {
   EnergyFluxDeriv my_functor(*m_spectrum, parName);
   return computeEnergyIntegral(my_functor, m_energies);
}

double DiffuseSource::energyFlux(double emin, double emax, size_t npts) const {
   EnergyFlux my_functor(*m_spectrum);
   return computeEnergyIntegral(my_functor, emin, emax, npts);
}

double DiffuseSource::energyFluxDeriv(const std::string & parName,
                                      double emin, double emax, 
                                      size_t npts) const {
   EnergyFluxDeriv my_functor(*m_spectrum, parName);
   return computeEnergyIntegral(my_functor, emin, emax, npts);
}

double DiffuseSource::diffuseResponse(const Event & evt) const {
   double trueEnergy(evt.getEnergy());
   const ResponseFunctions & respFuncs(m_observation->respFuncs());
   const WcsMap2 & wcsmap(mapBaseObject()->wcsmap());
   const std::vector< std::vector<float> > & solidAngles(wcsmap.solidAngles());
   double my_value(0);
   for (size_t i(0); i < solidAngles.size(); i++) {
      for (size_t j(0); j < solidAngles.at(i).size(); j++) {
         double psf_range(psfRange(evt.getEnergy()));
         if (evt.getDir().difference(wcsmap.skyDir(i+1, j+1)) < psf_range) {
            // WcsMap::skyDir uses wcslib pixel numbering, i.e., starting with 1
            astro::SkyDir srcDir(wcsmap.skyDir(i+1, j+1));
            double mapValue(spatialDist(SkyDirArg(srcDir, trueEnergy)));
            my_value += (respFuncs.totalResponse(trueEnergy, evt.getEnergy(), 
                                                 evt.zAxis(), evt.xAxis(),
                                                 srcDir, evt.getDir(),
                                                 evt.getType(),
                                                 evt.getArrTime())
                         *mapValue*solidAngles.at(i).at(j));
         }
      }
   }
   return my_value;
}

bool DiffuseSource::mapBasedIntegral() const {
   return m_mapBasedIntegral;
}

double DiffuseSource::psfRange(double energy) const {
   static double min_angle(5.*M_PI/180.);
   static double ebreak(1e4);
   if (energy > ebreak) {
      /// PSF containment flattens above ebreak.
      return min_angle;
   }
   /// Empirical power-law scaling.
   return min_angle*std::pow(energy/ebreak, -0.8);
}

} // namespace Likelihood
