/** 
 * @file DiffuseSource.cxx
 * @brief DiffuseSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/DiffuseSource.cxx,v 1.34 2006/12/19 20:59:47 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <string>
#include <valarray>
#include <vector>

#include "optimizers/Function.h"
#include "optimizers/dArg.h"
#include "optimizers/Exception.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Event.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/Observation.h"
#include "Likelihood/TrapQuad.h"

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

double DiffuseSource::fluxDensity(const Event &evt) const {
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
                                       const std::string &paramName) const {
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
   
double DiffuseSource::pixelCounts(double emin, double emax,
                                  double wtMin, double wtMax) const {
   optimizers::Function & spectrum = *m_spectrum;
   optimizers::dArg eminArg(emin);
   optimizers::dArg emaxArg(emax);
   return (spectrum(emaxArg)*wtMax + spectrum(eminArg)*wtMin)*(emax - emin)/2.;
}

double DiffuseSource::pixelCountsDeriv(double emin, double emax,
                                       double wtMin, double wtMax,
                                       const std::string & paramName) const {
   optimizers::Function & spectrum = *m_spectrum;
   optimizers::dArg eminArg(emin);
   optimizers::dArg emaxArg(emax);
   return (spectrum.derivByParam(emaxArg, paramName)*wtMax +
           spectrum.derivByParam(eminArg, paramName)*wtMin)*(emax - emin)/2.;
}

} // namespace Likelihood
