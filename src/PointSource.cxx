/** 
 * @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PointSource.cxx,v 1.29 2003/10/25 00:22:52 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>

#include "astro/SkyDir.h"

#include "optimizers/dArg.h"

#include "latResponse/Irfs.h"
#include "latResponse/../src/Glast25.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/TrapQuad.h"

namespace {

   double totalResponse(double energy, double time, 
                        const astro::SkyDir &srcDir,
                        const astro::SkyDir &appDir) {
// This implementation neglects energy dispersion.
      Likelihood::ResponseFunctions * respFuncs 
         = Likelihood::ResponseFunctions::instance();
      Likelihood::ScData * scData = Likelihood::ScData::instance();
   
      astro::SkyDir zAxis = scData->zAxis(time);
      astro::SkyDir xAxis = scData->xAxis(time);

      double myResponse = 0;
      std::map<unsigned int, latResponse::Irfs *>::iterator respIt
         = respFuncs->begin();
      for ( ; respIt != respFuncs->end(); respIt++) {
         
         latResponse::IPsf *psf = respIt->second->psf();
         latResponse::IAeff *aeff = respIt->second->aeff();

         double psf_val = psf->value(appDir, energy, srcDir, zAxis, xAxis);
         double aeff_val = aeff->value(energy, srcDir, zAxis, xAxis);
         myResponse += psf_val*aeff_val;
      }
      return myResponse;
   }

   double sourceEffArea(double energy, double time,
                        const astro::SkyDir &srcDir) {
      Likelihood::ResponseFunctions * respFuncs 
         = Likelihood::ResponseFunctions::instance();
      Likelihood::ScData * scData = Likelihood::ScData::instance();
   
      astro::SkyDir zAxis = scData->zAxis(time);
      astro::SkyDir xAxis = scData->xAxis(time);

      Likelihood::RoiCuts *roiCuts = Likelihood::RoiCuts::instance();
      std::vector<latResponse::AcceptanceCone *> cones;
      cones.push_back(const_cast<latResponse::AcceptanceCone *>
                      (&(roiCuts->extractionRegion())));

      double myEffArea = 0;
      std::map<unsigned int, latResponse::Irfs *>::iterator respIt
         = respFuncs->begin();
      for ( ; respIt != respFuncs->end(); respIt++) {
         
         latResponse::IPsf *psf = respIt->second->psf();
         latResponse::IAeff *aeff = respIt->second->aeff();

         double psf_val = psf->angularIntegral(energy, srcDir, 
                                               zAxis, xAxis, cones);
         double aeff_val = aeff->value(energy, srcDir, zAxis, xAxis);
         myEffArea += psf_val*aeff_val;
      }
      return myEffArea;
   }

}

namespace Likelihood {

bool PointSource::s_haveStaticMembers = false;
std::vector<double> PointSource::s_energies;

PointSource::PointSource(const PointSource &rhs) : Source(rhs) {
// make a deep copy
   m_dir = rhs.m_dir;
   m_functions["Position"] = &m_dir;

   m_spectrum = rhs.m_spectrum->clone();
   m_functions["Spectrum"] = m_spectrum;

   m_exposure = rhs.m_exposure;
   m_srcType = rhs.m_srcType;
}

double PointSource::fluxDensity(double energy, double time,
                                const astro::SkyDir &dir) const {

// Scale the energy spectrum by the psf value and the effective area
// and convolve with the energy dispersion (now a delta-function in
// energy), all of which are functions of time and spacecraft attitude
// and orbital position.

   optimizers::dArg energy_arg(energy);
   double spectrum = (*m_spectrum)(energy_arg);

   return ::totalResponse(energy, time, m_dir.getDir(), dir)*spectrum;
}

double PointSource::fluxDensityDeriv(double energy, double time,
                                     const astro::SkyDir &dir,
                                     const std::string &paramName) const {
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude"

   if (paramName == "Prefactor") {
      return fluxDensity(energy, time, dir)
         /m_spectrum->getParamValue("Prefactor");
   } else {
      optimizers::dArg energy_arg(energy);
      return ::totalResponse(energy, time, m_dir.getDir(), dir)
         *m_spectrum->derivByParam(energy_arg, paramName);
   }
}

double PointSource::Npred() {
   optimizers::Function *specFunc = m_functions["Spectrum"];

// Evaluate the Npred integrand at the abscissa points contained in
// s_energies
   
   std::vector<double> NpredIntegrand(s_energies.size());
   for (unsigned int k = 0; k < s_energies.size(); k++) {
      optimizers::dArg eArg(s_energies[k]);
      NpredIntegrand[k] = (*specFunc)(eArg)*m_exposure[k];
   }
   TrapQuad trapQuad(s_energies, NpredIntegrand);
   return trapQuad.integral();
}

double PointSource::NpredDeriv(const std::string &paramName) {
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

void PointSource::computeExposure(int verbose) {
   computeExposure(s_energies, m_exposure, verbose);
}

void PointSource::computeExposure(std::vector<double> &energies,
                                  std::vector<double> &exposure,
                                  int verbose) {
   if (!s_haveStaticMembers) {
      makeEnergyVector();
      s_haveStaticMembers = true;
   }

   ScData *scData = ScData::instance();
   RoiCuts *roiCuts = RoiCuts::instance();

// Initialize the exposure vector with zeros
   exposure = std::vector<double>(energies.size(), 0);

   if (verbose) {
      std::cerr << "Computing exposure at (" 
                << getDir().ra() << ", " 
                << getDir().dec() << ")";
   }
   unsigned int npts = scData->vec.size()-1;
   for (unsigned int it = 0; it < npts; it++) {
      if (npts/20 > 0 && ((it % (npts/20)) == 0) && verbose) std::cerr << ".";
      bool includeInterval = true;

// Check if this interval passes the time cuts
      std::vector< std::pair<double, double> > timeCuts;
      roiCuts->getTimeCuts(timeCuts);
      for (unsigned int itcut = 0; itcut < timeCuts.size(); itcut++) {
         if (scData->vec[it].time < timeCuts[itcut].first ||
             scData->vec[it].time > timeCuts[itcut].second) {
            includeInterval = false;
            break;
         }
      }

// Check for SAA passage
      if (scData->vec[it].inSaa) includeInterval = false;

// Compute the inclination and check if it's within response matrix
// cut-off angle
      double inc = getSeparation(scData->vec[it].zAxis)*180/M_PI;
      if (inc > latResponse::Glast25::incMax()) includeInterval = false;

// Having checked for relevant constraints, add the exposure
// contribution for each energy
      if (includeInterval) {
         for (unsigned int k = 0; k < energies.size(); k++) {
            double time = (scData->vec[it+1].time + scData->vec[it].time)/2.;
            exposure[k] += 
               ::sourceEffArea(energies[k], time, m_dir.getDir())
               *(scData->vec[it+1].time - scData->vec[it].time);
         }
      }
   }
   if (verbose) std::cerr << "!" << std::endl;
}

void PointSource::makeEnergyVector(int nee) {
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
