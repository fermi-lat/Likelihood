/** @file DiffuseSource.cxx
 * @brief DiffuseSource class implementation
 *
 * $Header$
 */

#include <vector>
#include <string>
#include <cmath>

#include "astro/SkyDir.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/Psf.h"
#include "Likelihood/Aeff.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/dArg.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/ExposureMap.h"

namespace Likelihood {

bool DiffuseSource::s_haveStaticMembers = false;
std::vector<double> DiffuseSource::s_energies;

DiffuseSource::DiffuseSource(Function* spatialDist) : m_spectrum(0) {
// The spatial distribution of emission is required for instantiation.
   m_spatialDist = spatialDist->clone();
   m_functions["SpatialDist"] = m_spatialDist;

// Exposure is computed automatically, so RoiCuts and spacecraft data
// must be available; furthermore, the ExposureMap object must have
// been instantiated.
   computeExposure();
}

DiffuseSource::DiffuseSource(const DiffuseSource &rhs) : Source(rhs) {
// make a deep copy
   m_spatialDist = rhs.m_spatialDist->clone();
   m_functions["SpatialDist"] = m_spatialDist;

   m_spectrum = rhs.m_spectrum->clone();
   m_functions["Spectrum"] = m_spectrum;
}

double DiffuseSource::fluxDensity(const Event &evt) const {

// This trueEnergy is just a place-holder until we have a finite
// energy dispersion.  In principle, we want to return the convolution
// of the evt.diffuseResponse with m_spectrum.

   double trueEnergy = evt.getEnergy(); 

// Note that the Event-specific diffuseResponses are assumed to
// be identified by the DiffuseSource name.
   dArg energy_arg(evt.getEnergy());
   return (*m_spectrum)(energy_arg)*evt.diffuseResponse(trueEnergy, getName());
}

double DiffuseSource::fluxDensityDeriv(const Event &evt, 
                                       std::string &paramName) const {
                                   
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude"

   if (paramName == "Prefactor") {
      return fluxDensity(evt)/m_spectrum->getParamValue("Prefactor");
   } else {
      double trueEnergy = evt.getEnergy(); 
      dArg energy_arg(evt.getEnergy());
      return m_spectrum->derivByParam(energy_arg, paramName)
         *evt.diffuseResponse(trueEnergy, getName());
   }
}

double DiffuseSource::Npred() {
   Function *specFunc = m_functions["Spectrum"];

// evaluate the Npred integrand at the abscissa points contained in
// s_energies
   
   std::vector<double> NpredIntegrand(s_energies.size());
   for (unsigned int k = 0; k < s_energies.size(); k++) {
      dArg eArg(s_energies[k]);
      NpredIntegrand[k] = (*specFunc)(eArg)*m_exposure[k];
   }
   TrapQuad trapQuad(s_energies, NpredIntegrand);
   return trapQuad.integral();
}

double DiffuseSource::NpredDeriv(const std::string &paramName) {
   Function *specFunc = m_functions["Spectrum"];

   if (paramName == std::string("Prefactor")) {
      return Npred()/specFunc->getParamValue("Prefactor");
   } else {  // loop over energies and fill integrand vector
      std::vector<double> myIntegrand(s_energies.size());
      for (unsigned int k = 0; k < s_energies.size(); k++) {
         dArg eArg(s_energies[k]);
         myIntegrand[k] = specFunc->derivByParam(eArg, paramName)
            *m_exposure[k];
      }
      TrapQuad trapQuad(s_energies, myIntegrand);
      return trapQuad.integral();
   }
}

void DiffuseSource::computeExposure() {
   if (!s_haveStaticMembers) {
      makeEnergyVector();
      s_haveStaticMembers = true;
   }         

   ExposureMap *exposureMap = ExposureMap::instance();

// Fetch the pixel coordinates.
   std::vector<double> ra;
   exposureMap->fetchRA(ra);
   std::vector<double> dec;
   exposureMap->fetchDec(dec);

// Fetch the exposure multiplied by the solid angle of the associated
// pixel.
   std::vector<double> exposure;
   exposureMap->fetchExposure(exposure);

   m_exposure.clear();
   m_exposure.reserve(s_energies.size());
   for (unsigned int k = 0; k < s_energies.size(); k++) {
      double srcExposure = 0;
      for (unsigned int j = 0; j < ra.size(); j++) {
         astro::SkyDir skyDir(ra[j], dec[j]);
         SkyDirArg dir(skyDir);
         srcExposure += exposure[j]*(*m_spatialDist)(dir);
      }
      m_exposure.push_back(srcExposure);
   }
}

void DiffuseSource::makeEnergyVector(int nee) {
   RoiCuts *roiCuts = RoiCuts::instance();
   
// set up a logrithmic grid of energies for doing the integral over 
// the spectrum
   double emin = (roiCuts->getEnergyCuts()).first;
   double emax = (roiCuts->getEnergyCuts()).second;
   double estep = log(emax/emin)/(nee-1);
   
   s_energies.reserve(nee);
   for (int i = 0; i < nee; i++)
      s_energies.push_back(emin*exp(i*estep));
}

} // namespace Likelihood
