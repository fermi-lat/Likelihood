/** @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header:
 */

#include <vector>
#include <string>
#include <cmath>

#include "astro/SkyDir.h"
#include "../Likelihood/PointSource.h"
#include "../Likelihood/Psf.h"
#include "../Likelihood/Aeff.h"
#include "../Likelihood/ScData.h"
#include "../Likelihood/RoiCuts.h"
#include "../Likelihood/dArg.h"
#include "../Likelihood/TrapQuad.h"

namespace Likelihood {

bool PointSource::m_haveStaticMembers = false;
std::vector<double> PointSource::m_energies;
std::vector<double> PointSource::m_sigGauss;

void PointSource::m_makeEnergyVector(int nee) {
   RoiCuts *roiCuts = RoiCuts::instance();
   
// set up a logrithmic grid of energies for doing the integral over 
// the spectrum
   double emin = (roiCuts->getEnergyCuts()).first;
   double emax = (roiCuts->getEnergyCuts()).second;
   double estep = log(emax/emin)/(nee-1);
   
   m_energies.reserve(nee);
   for (int i = 0; i < nee; i++)
      m_energies.push_back(emin*exp(i*estep));
}

void PointSource::m_makeSigmaVector(int nsig) {
// set up a grid of sigma values in radians
   double sigmin = 0.;
   double sigmax = 50.*M_PI/180;
   double sigstep = (sigmax - sigmin)/(nsig - 1);
   
   m_sigGauss.reserve(nsig);
   for (int i = 0; i < nsig; i++)
      m_sigGauss.push_back(sigstep*i + sigmin);
}

PointSource::PointSource(const PointSource &rhs) : Source(rhs) {
   m_dir = rhs.m_dir;
   m_spectrum = rhs.m_spectrum;
   m_energies = rhs.m_energies;
}

double PointSource::fluxDensity(double energy, double time,
                                const astro::SkyDir &dir) const {

// Scale the energy spectrum by the psf value and the effective area
// and convolve with the energy dispersion (now a delta-function in
// energy), all of which are functions of time and spacecraft attitude
// and orbital position.

   Psf *psf = Psf::instance();
   Aeff *aeff = Aeff::instance();

   dArg energy_arg(energy);
   double spectrum = (*m_spectrum)(energy_arg);
   double psf_val = (*psf)(dir, energy, m_dir.getDir(), time);
   double aeff_val = (*aeff)(energy, dir, time);
   return spectrum*psf_val*aeff_val;
}

double PointSource::Npred() {
// presently implemented without psf cuts....

   Function *specFunc = m_functions["Spectrum"];

// evaluate the Npred integrand at the abscissa points contained in
// m_energies
   
   std::vector<double> NpredIntegrand(m_energies.size());
   for (unsigned int k = 0; k < m_energies.size(); k++) {
      dArg eArg(m_energies[k]);
      NpredIntegrand[k] = (*specFunc)(eArg)*m_exposure[k];
   }
   TrapQuad trapQuad(m_energies, NpredIntegrand);
   return trapQuad.integral();
}

double PointSource::NpredDeriv(const std::string &paramName) {
   Function *specFunc = m_functions["Spectrum"];

   if (paramName == std::string("Prefactor")) {
      return Npred()/specFunc->getParamValue("Prefactor");
   } else {  // loop over energies and fill integrand vector
      std::vector<double> myIntegrand(m_energies.size());
      for (unsigned int k = 0; k < m_energies.size(); k++) {
         dArg eArg(m_energies[k]);
         myIntegrand[k] = specFunc->derivByParam(eArg, paramName);
      }
      TrapQuad trapQuad(m_energies, myIntegrand);
      return trapQuad.integral();
   }
}

void PointSource::computeExposure() {
   Aeff *aeff = Aeff::instance();
   ScData *scData = ScData::instance();
   RoiCuts *roiCuts = RoiCuts::instance();

// initialize the exposure vector with zeros
   m_exposure = std::vector<double>(m_energies.size(), 0);

   for (unsigned int it = 0; it < scData->vec.size()-1; it++) {
      bool includeInterval = true;

// check if this interval passes the time cuts
      std::vector< std::pair<double, double> > timeCuts;
      roiCuts->getTimeCuts(timeCuts);
      for (unsigned int itcut = 0; itcut < timeCuts.size(); itcut++) {
	 if (scData->vec[it].time < timeCuts[itcut].first ||
             scData->vec[it].time > timeCuts[itcut].second) {
            includeInterval = false;
            break;
         }
      }

// check for SAA passage
      if (scData->vec[it].inSaa) includeInterval = false;

// compute inclination and check if within response matrix cut-off angle
      double inc = getSeparation(scData->vec[it].zAxis)*180/M_PI;
      if (inc > Response::incMax) includeInterval = false;

// having checked for relevant constraints, add up the exposure
// for each energy (still need to implement the psf fraction stuff....)
      if (includeInterval) {
         for (unsigned int k = 0; k < m_energies.size(); k++) {
            m_exposure[k] += (*aeff)(m_energies[k], inc)
               *(scData->vec[it+1].time - scData->vec[it].time);
         }
      }
   }
}

void PointSource::computeGaussFractions() {
// compute the fraction of a Gaussian centered on the PointSource position
// contained within the ROI as a function of Gaussian width sigma.

   RoiCuts *roiCuts = RoiCuts::instance();

   double psi = m_dir.SkyDir::difference(roiCuts->getExtractionRegion().first);
   double roi_radius = roiCuts->getExtractionRegion().second*M_PI/180;
   double mup = cos(roi_radius + psi);
   double mum = cos(roi_radius - psi);

   double cr = cos(roi_radius);
   double cp = cos(psi);
   double sp = sin(psi);

   m_gaussFraction.reserve(m_sigGauss.size());
   for (unsigned int i = 0; i < m_sigGauss.size(); i++) {
      double sig = m_sigGauss[i];
      if (m_sigGauss[i] == 0) {
         m_gaussFraction[i] = 0;
         break;
      }
      double denom = 1. - exp(-2/sig/sig);
      double gauss_int;
      if (psi == 0) {
         gauss_int = 0;
      } else {
// instantiate the integrand object...
	 Gint gfunc(m_sigGauss[i], cr, cp, sp);
// and the integrator object
         TrapQuad trapQuad(&gfunc);
         gauss_int = trapQuad.integral(mum, mup, 10000);
      }
      if (psi <= roi_radius) {
	 m_gaussFraction.push_back( ((1 - exp((mum-1)/sig/sig))
				     + gauss_int/M_PI/sig/sig)/denom );
      } else {
	 m_gaussFraction.push_back( gauss_int/M_PI/sig/sig/denom );
      }
   }
}

double PointSource::Gint::value(Arg &muarg) const {
   double mu = dynamic_cast<dArg &>(muarg).getValue();

   double phi;
   if (mu == 1) {
      phi = M_PI;
   } else {
      phi = acos((m_cr - mu*m_cp)/sqrt(1. - mu*mu)/m_sp);
   }
   return phi*exp((mu - 1)/m_sig/m_sig);
}

} // namespace Likelihood
