/** 
 * @file DiffuseSource.cxx
 * @brief DiffuseSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/DiffuseSource.cxx,v 1.22 2005/02/27 06:42:25 jchiang Exp $
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
#include "Likelihood/RoiCuts.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/ResponseFunctions.h"

namespace Likelihood {

bool DiffuseSource::s_haveStaticMembers = false;
std::vector<double> DiffuseSource::s_energies;

DiffuseSource::DiffuseSource(optimizers::Function * spatialDist,
                             const Observation & observation,
                             bool requireExposure) : m_spectrum(0) {
   m_spatialDist = spatialDist->clone();
   m_functions["SpatialDist"] = m_spatialDist;
   m_useEdisp = observation.respFuncs().useEdisp();

   double emin = observation.roiCuts().getEnergyCuts().first;
   double emax = observation.roiCuts().getEnergyCuts().second;
   if (!s_haveStaticMembers || emin != s_energies.front() ||
       emax != s_energies.back()) {
      makeEnergyVector(emin, emax);
      s_haveStaticMembers = true;
   }

// In order to compute exposure, RoiCuts and spacecraft data must be
// available; furthermore, the ExposureMap object must have been
// created.
   if (requireExposure) {
      observation.expMap().integrateSpatialDist(s_energies, spatialDist,
                                                m_exposure);
   }
   m_srcType = "Diffuse";
}

DiffuseSource::DiffuseSource(const DiffuseSource &rhs) : Source(rhs) {
// make a deep copy
   m_spatialDist = rhs.m_spatialDist->clone();
   m_functions["SpatialDist"] = m_spatialDist;

   m_spectrum = rhs.m_spectrum->clone();
   m_functions["Spectrum"] = m_spectrum;

   m_exposure = rhs.m_exposure;
   m_srcType = rhs.m_srcType;
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
// s_energies.
   
   std::vector<double> NpredIntegrand(s_energies.size());
   for (unsigned int k = 0; k < s_energies.size(); k++) {
      optimizers::dArg eArg(s_energies[k]);
      NpredIntegrand[k] = (*specFunc)(eArg)*m_exposure[k];
   }
   TrapQuad trapQuad(s_energies, NpredIntegrand);
   return trapQuad.integral();
}

double DiffuseSource::NpredDeriv(const std::string &paramName) {
   optimizers::Function *specFunc = m_functions["Spectrum"];

   if (paramName == std::string("Prefactor")) {
      return Npred()/specFunc->getParamValue("Prefactor");
   } else {  // loop over energies and fill integrand vector
      std::vector<double> myIntegrand(s_energies.size());
      for (unsigned int k = 0; k < s_energies.size(); k++) {
         optimizers::dArg eArg(s_energies[k]);
         myIntegrand[k] = specFunc->derivByParam(eArg, paramName)
            *m_exposure[k];
      }
      TrapQuad trapQuad(s_energies, myIntegrand);
      return trapQuad.integral();
   }
}

double DiffuseSource::Npred(double emin, double emax) {
   if (emin < s_energies.front() || emax > s_energies.back()) {
      throw std::out_of_range("DiffuseSource::Npred(emin, emax)");
   }
   std::vector<double>::iterator first 
      = std::upper_bound(s_energies.begin(), s_energies.end(), emin);
   std::vector<double>::iterator last 
      = std::upper_bound(s_energies.begin(), s_energies.end(), emax);
   std::vector<double> energies(last - first);
   std::copy(first, last, energies.begin());
   int begin_offset = first - s_energies.begin();
   int end_offset = last - s_energies.begin();
   energies.insert(energies.begin(), emin);
   energies.push_back(emax);
   std::vector<double> exposure(last - first);
   std::copy(m_exposure.begin() + begin_offset,
             m_exposure.begin() + end_offset,
             exposure.begin());
   double begin_exposure = (emin - s_energies[begin_offset - 1])
      /(s_energies[begin_offset] - s_energies[begin_offset - 1])
      *(m_exposure[begin_offset] - m_exposure[begin_offset - 1])
      + m_exposure[begin_offset - 1];
   double end_exposure = (emin - s_energies[end_offset - 1])
      /(s_energies[end_offset] - s_energies[end_offset - 1])
      *(m_exposure[end_offset] - m_exposure[end_offset - 1])
      + m_exposure[end_offset - 1];
   exposure.insert(exposure.begin(), begin_exposure);
   exposure.push_back(end_exposure);
   optimizers::Function & specFunc = *m_functions["Spectrum"];
   std::vector<double> integrand(energies.size());
   for (unsigned int k = 0; k < energies.size(); k++) {
      optimizers::dArg eArg(energies[k]);
      integrand[k] = specFunc(eArg)*exposure[k];
   }
   TrapQuad trapQuad(energies, integrand);
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

void DiffuseSource::makeEnergyVector(double emin, double emax, int nee) {
   double estep = log(emax/emin)/(nee-1);
   s_energies.clear();
   s_energies.reserve(nee);
   for (int i = 0; i < nee; i++) {
      s_energies.push_back(emin*exp(i*estep));
   }
}

} // namespace Likelihood
