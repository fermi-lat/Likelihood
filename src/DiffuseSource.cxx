/** 
 * @file DiffuseSource.cxx
 * @brief DiffuseSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/DiffuseSource.cxx,v 1.25 2005/03/02 01:10:53 jchiang Exp $
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
   : m_spectrum(0), m_observation(&observation) {
   m_spatialDist = spatialDist->clone();
   m_functions["SpatialDist"] = m_spatialDist;
   m_useEdisp = observation.respFuncs().useEdisp();

// In order to compute exposure, RoiCuts and spacecraft data must be
// available; furthermore, the ExposureMap object must have been
// created.
   if (requireExposure) {
      RoiCuts & roiCuts = const_cast<RoiCuts &>(observation.roiCuts());
      const std::vector<double> & energies = roiCuts.energies();
      observation.expMap().integrateSpatialDist(energies, spatialDist,
                                                m_exposure);
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
   if (paramName == "Prefactor") {
      my_fluxDensityDeriv 
         = fluxDensity(evt)/m_spectrum->getParamValue("Prefactor");
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

double DiffuseSource::Npred() {
   optimizers::Function *specFunc = m_functions["Spectrum"];

// Evaluate the Npred integrand at the abscissa points contained in
// RoiCuts::energies().

   RoiCuts & roiCuts = const_cast<RoiCuts &>(m_observation->roiCuts());
   const std::vector<double> & energies = roiCuts.energies();

   std::vector<double> NpredIntegrand(energies.size());
   for (unsigned int k = 0; k < energies.size(); k++) {
      optimizers::dArg eArg(energies[k]);
      NpredIntegrand[k] = (*specFunc)(eArg)*m_exposure[k];
   }
   TrapQuad trapQuad(energies, NpredIntegrand);
   return trapQuad.integral();
}

double DiffuseSource::NpredDeriv(const std::string &paramName) {
   optimizers::Function * specFunc = m_functions["Spectrum"];

   const std::vector<double> & energies = 
      const_cast<RoiCuts &>(m_observation->roiCuts()).energies();

   if (paramName == std::string("Prefactor")) {
      return Npred()/specFunc->getParamValue("Prefactor");
   } else {  // loop over energies and fill integrand vector
      std::vector<double> myIntegrand(energies.size());
      for (unsigned int k = 0; k < energies.size(); k++) {
         optimizers::dArg eArg(energies[k]);
         myIntegrand[k] = specFunc->derivByParam(eArg, paramName)
            *m_exposure[k];
      }
      TrapQuad trapQuad(energies, myIntegrand);
      return trapQuad.integral();
   }
}

double DiffuseSource::Npred(double emin, double emax) {
   const std::vector<double> & energies = 
      const_cast<RoiCuts &>(m_observation->roiCuts()).energies();

   if (emin < energies.front() || emax > energies.back()) {
      throw std::out_of_range("DiffuseSource::Npred(emin, emax)");
   }
   std::vector<double>::const_iterator first 
      = std::upper_bound(energies.begin(), energies.end(), emin);
   std::vector<double>::const_iterator last 
      = std::upper_bound(energies.begin(), energies.end(), emax);
   std::vector<double> my_energies(last - first);
   std::copy(first, last, my_energies.begin());
   int begin_offset = first - energies.begin();
   int end_offset = last - energies.begin();
   my_energies.insert(my_energies.begin(), emin);
   my_energies.push_back(emax);
   std::vector<double> exposure(last - first);
   std::copy(m_exposure.begin() + begin_offset,
             m_exposure.begin() + end_offset,
             exposure.begin());
   double begin_exposure = (emin - energies[begin_offset - 1])
      /(energies[begin_offset] - energies[begin_offset - 1])
      *(m_exposure[begin_offset] - m_exposure[begin_offset - 1])
      + m_exposure[begin_offset - 1];
   double end_exposure = (emin - energies[end_offset - 1])
      /(energies[end_offset] - energies[end_offset - 1])
      *(m_exposure[end_offset] - m_exposure[end_offset - 1])
      + m_exposure[end_offset - 1];
   exposure.insert(exposure.begin(), begin_exposure);
   exposure.push_back(end_exposure);
   optimizers::Function & specFunc = *m_functions["Spectrum"];
   std::vector<double> integrand(my_energies.size());
   for (unsigned int k = 0; k < my_energies.size(); k++) {
      optimizers::dArg eArg(my_energies[k]);
      integrand[k] = specFunc(eArg)*exposure[k];
   }
   TrapQuad trapQuad(my_energies, integrand);
   return trapQuad.integral();
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
