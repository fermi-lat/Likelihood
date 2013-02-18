/** 
 * @file DiffuseSource.cxx
 * @brief DiffuseSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/DiffuseSource.cxx,v 1.61 2013/01/09 00:44:41 jchiang Exp $
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

void PointSource::computeResponse(Response& resp, const Event & evt) const {
   const double& energy         = evt.getEnergy();

   double minusone(-1);
   double one(1);
   double mumin(minusone);
   double mumax(one);

   EquinoxRotation eqRot(dir.ra(), dir.dec());
   try {
      srcs.at(i)->mapBaseObject()->getDiffRespLimits(getDir(), 
						     mumin, mumax,
						     phimin, phimax);
   } catch (MapBaseException &) {
      // do nothing
   }

   if(m_edisp() == ED_NONE) {
      resp.resize(1);
      resp[0] = calculateRespValue(event, energy, /* useEdisp = */ false,
				   eqRot, mumin, mumax, phimin, phimax);
   } else {
      const std::vector<double> & true_energies = evt.trueEnergies();
      resp.resize(true_energies.size());
      for(unsigned ienergy = 0; ienergy<resp.size(); ienergy++) {
	resp[ienergy] = 
	  calculateRespValue(event, energy, /* useEdisp = */ false,
			     eqRot, mumin, mumax, phimin, phimax);
      }
      if(m_edisp() == ED_GAUSSIAN) {
	 computeGaussianParams(resp, resp, true_energies);
      }
   }
}

double DiffuseSource::
calculateRespValue(const Event& event, double trueEnergy, bool useEdisp,
		   const EquinoxRotation& eqRot, double mumin, double mumax,
		   double phimin, double phimax) const
{
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
   double respValue(0);
   if (mapBasedIntegral() || 
       (::getenv("MAP_BASED_DIFFRSP") 
	&& (mumin != minusone || mumax != one))) {
      respValue = calculateMapBasedRespValue(event, trueEnergy, useEdisp);
   } else if (::getenv("USE_OLD_DIFFRSP")) {
      /// Old integration scheme with the phi integral
      /// evaluated inside the theta integral.
      respValue = DiffRespIntegrand::
	do2DIntegration(*this, respFuncs, *srcs.at(i), eqRot, 
			trueEnergy, useEdisp,
			mumin, mumax, phimin, phimax, 0.01, 0.1);
   } else {
      /// Steve's integration scheme with the theta integral
      /// evaluated inside the phi integral.  The produces
      /// much more accurate results.
      respValue = DiffRespIntegrand2::
	do2DIntegration(*this, respFuncs, *srcs.at(i), eqRot,
			trueEnergy, useEdisp,
			mumin, mumax, phimin, phimax, 0.001, 0.01);
   }
   return respValue;
}

double DiffuseSource::
calculateMapBasedRespValue(const Event & evt,
			   double trueEnergy, bool useEdisp) const {
   const ResponseFunctions & respFuncs(m_observation->respFuncs());
   const WcsMap2 & wcsmap(mapBaseObject()->wcsmap());
   const std::vector< std::vector<double> > & solidAngles(wcsmap.solidAngles());
   double my_value(0);
   for (size_t i(0); i < solidAngles.size(); i++) {
      for (size_t j(0); j < solidAngles.at(i).size(); j++) {
         // WcsMap::skyDir uses wcslib pixel numbering, i.e., starting with 1
         astro::SkyDir srcDir(wcsmap.skyDir(i+1, j+1));
         double mapValue(spatialDist(SkyDirArg(srcDir, trueEnergy)));
	 if(mapValue != 0) {
	    // This test potentially saves multiple calls to calculate
	    // the response functions
	    my_value += 
	      respFuncs.totalResponse(trueEnergy, evt.getEnergy(), 
				      evt.zAxis(), evt.xAxis(), 
				      srcDir, evt.getDir(),
				      evt.getType(), evt.getArrTime(),
				      useEdisp)
	      *mapValue*solidAngles.at(i).at(j));
	 }
      }
   }
   return my_value;
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

bool DiffuseSource::mapBasedIntegral() const {
   return m_mapBasedIntegral;
}

} // namespace Likelihood
