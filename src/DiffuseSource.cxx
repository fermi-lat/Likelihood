/** 
 * @file DiffuseSource.cxx
 * @brief DiffuseSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/DiffuseSource.cxx,v 1.49 2009/03/26 01:32:43 jchiang Exp $
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
#include "Likelihood/Observation.h"

namespace Likelihood {

DiffuseSource::DiffuseSource(optimizers::Function * spatialDist,
                             const Observation & observation,
                             bool requireExposure) 
   : Source(&observation) {
   m_spatialDist = spatialDist->clone();
   m_functions["SpatialDist"] = m_spatialDist;
   m_useEdisp = observation.respFuncs().useEdisp();

// In order to compute exposure, RoiCuts and spacecraft data must be
// available; furthermore, the ExposureMap object must have been
// created.
   if (requireExposure) {
      const RoiCuts & roiCuts = observation.roiCuts();
      const std::vector<double> & energies = roiCuts.energies();
      if (observation.expMap().haveMap()) {
         observation.expMap().integrateSpatialDist(energies, spatialDist,
                                                   m_exposure);
      } else {
         std::string what("DiffuseSource: An exposure map must be defined ");
         what += "if diffuse sources are in the model.";
         throw std::runtime_error(what);
      }
   }
   m_srcType = "Diffuse";
}

DiffuseSource::DiffuseSource(const DiffuseSource &rhs) : Source(rhs) {
   m_spatialDist = rhs.m_spatialDist->clone();
   m_functions["SpatialDist"] = m_spatialDist;

   m_spectrum = rhs.m_spectrum->clone();
   m_functions["Spectrum"] = m_spectrum;

   m_exposure = rhs.m_exposure;

   m_observation = rhs.m_observation;
}

double DiffuseSource::fluxDensity(const Event &evt,
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
      TrapQuad trapQuad(trueEnergies, my_integrand);
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
         TrapQuad trapQuad(trueEnergies, my_integrand);
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
      throw MapBaseException("Flux calculations are not available for this "
                             + ("diffuse source: " + getName()));
   }
   return mapBaseObject;
}

double DiffuseSource::angularIntegral(double energy) const {
   if (spatialDist()->genericName() == "ConstantValue") { 
// Here we have an isotropic source
      return 4*M_PI;
   }
   return mapBaseObject()->mapIntegral(energy);
}

double DiffuseSource::flux() const {
   return computeEnergyIntegral(*m_spectrum, 
                                m_observation->roiCuts().energies());
}

double DiffuseSource::fluxDeriv(const std::string & parName) const {
   FluxDeriv my_functor(*m_spectrum, parName);
   return computeEnergyIntegral(my_functor,
                                m_observation->roiCuts().energies());
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
   return computeEnergyIntegral(my_functor, 
                                m_observation->roiCuts().energies());
}

double DiffuseSource::energyFluxDeriv(const std::string & parName) const {
   EnergyFluxDeriv my_functor(*m_spectrum, parName);
   return computeEnergyIntegral(my_functor, 
                                m_observation->roiCuts().energies());
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
   const WcsMap & wcsmap(mapBaseObject()->wcsmap());
   const std::vector< std::vector<double> > & solidAngles(wcsmap.solidAngles());
   double my_value(0);
   for (size_t i(0); i < solidAngles.size(); i++) {
      for (size_t j(0); j < solidAngles.at(i).size(); j++) {
         // WcsMap::skyDir uses wcslib pixel numbering, i.e., starting with 1
         astro::SkyDir srcDir(wcsmap.skyDir(i+1, j+1));
         double mapValue(spatialDist(SkyDirArg(srcDir, trueEnergy)));
         my_value += (respFuncs.totalResponse(trueEnergy, evt.getEnergy(), 
                                              evt.zAxis(), evt.xAxis(),
                                              srcDir, evt.getDir(),
                                              evt.getType())
                      *mapValue*solidAngles.at(i).at(j));
      }
   }
   return my_value;
}

} // namespace Likelihood
