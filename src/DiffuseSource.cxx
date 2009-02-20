/** 
 * @file DiffuseSource.cxx
 * @brief DiffuseSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/DiffuseSource.cxx,v 1.44 2009/02/18 18:13:38 jchiang Exp $
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
   
double DiffuseSource::pixelCounts(double emin, double emax,
                                  double wtMin, double wtMax) const {
   optimizers::Function & spectrum = *m_spectrum;
   optimizers::dArg eminArg(emin);
   optimizers::dArg emaxArg(emax);

   double f1(spectrum(eminArg));
   double f2(spectrum(emaxArg));

   double y1(f1*wtMin);
   double y2(f2*wtMax);
   if (::getenv("USE_OLD_PIX_EST") || y1 == 0 || y2 == 0) {
      return (y1 + y2)*(emax - emin)/2.;
   }

   double gam(std::log(y2/y1)/std::log(emax/emin));
   double y0(y2/std::pow(emax, gam));
   if (gam == -1) {
      return y0*std::log(emax/emin);
   }
   return y0/(gam + 1.)*(std::pow(emax, gam + 1.) - std::pow(emin, gam + 1.));
}

double DiffuseSource::pixelCountsDeriv(double emin, double emax,
                                       double wtMin, double wtMax,
                                       const std::string & paramName) const {
   optimizers::Function & spectrum = *m_spectrum;
   optimizers::dArg eminArg(emin);
   optimizers::dArg emaxArg(emax);

   double y1(spectrum(eminArg)*wtMin);
   double y2(spectrum(emaxArg)*wtMax);

   double f1(spectrum.derivByParam(eminArg, paramName));
   double f2(spectrum.derivByParam(emaxArg, paramName));

   double dy1dp(f1*wtMin);
   double dy2dp(f2*wtMax);
   if (::getenv("USE_OLD_PIX_EST") || y1 == 0 || y2 == 0) {
      return (dy1dp + dy2dp)*(emax - emin)/2.;
   }

   double gam(std::log(y2/y1)/std::log(emax/emin));
   double y0(y2/std::pow(emax, gam));
   double dgamdp((dy2dp/y2 - dy1dp/y1)/std::log(emax/emin));
   double dy0dp((dy2dp - y2*dgamdp*std::log(emax))/std::pow(emax, gam));
   if (gam == -1) {
      return dy0dp*std::log(emax/emin);
   }
   return (dy0dp*(std::pow(emax, gam+1.) - std::pow(emin, gam+1.))/(gam+1.) +
           y0*dgamdp/(gam+1.)*((std::pow(emax, gam+1.)*std::log(emax)
                                - std::pow(emin, gam+1.)*std::log(emin))
                               - (std::pow(emax,gam+1.)-std::pow(emin,gam+1.))
                               /(gam+1.)));
}

double DiffuseSource::flux() const {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   std::vector<double> integrand;
   for (size_t k(0); k < energies.size(); k++) {
      
   }
   
   TrapQuad fluxIntegral(m_spectrum);
   return fluxIntegral.integral(energies);
   
}

} // namespace Likelihood
