/** 
 * @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PointSource.cxx,v 1.24 2003/10/22 04:30:33 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>

#include "astro/SkyDir.h"

#include "optimizers/dArg.h"

#include "latResponse/IPsf.h"
#include "latResponse/IAeff.h"
#include "latResponse/Irfs.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/Psf.h"
#include "Likelihood/Aeff.h"
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
}

namespace Likelihood {

typedef double (*D_fp)(double*);    // "from" f2c.h

extern "C" {
   int dgaus8_(D_fp fun, double *a, double *b, 
               double *err, double *ans, long *ierr);
}

PointSource::Gint PointSource::s_gfunc;
bool PointSource::s_haveStaticMembers = false;
std::vector<double> PointSource::s_energies;
std::vector<double> PointSource::s_sigGauss;

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
      makeSigmaVector();
      s_haveStaticMembers = true;
   }
   computeGaussFractions();

   Aeff *aeff = Aeff::instance();
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
      if (inc > Response::incMax()) includeInterval = false;

// Having checked for relevant constraints, add the exposure
// contribution for each energy
      if (includeInterval) {
         for (unsigned int k = 0; k < energies.size(); k++) {
            exposure[k] += (*aeff)(energies[k], inc)
               *psfFrac(energies[k], inc)
               *(scData->vec[it+1].time - scData->vec[it].time);
         }
      }
   }
   if (verbose) std::cerr << "!" << std::endl;
}

void PointSource::computeGaussFractions() {
// compute the fraction of a Gaussian centered on the PointSource position
// contained within the ROI as a function of Gaussian width sigma.

   RoiCuts *roiCuts = RoiCuts::instance();

   double psi = getSeparation(roiCuts->getExtractionRegion().first);
   double roi_radius = roiCuts->getExtractionRegion().second*M_PI/180;
   double mup = cos(roi_radius + psi);
   double mum = cos(roi_radius - psi);

   double cr = cos(roi_radius);
   double cp = cos(psi);
   double sp = sin(psi);

   m_gaussFraction.reserve(s_sigGauss.size());
   for (unsigned int i = 0; i < s_sigGauss.size(); i++) {
      double sig = s_sigGauss[i];
      if (s_sigGauss[i] == 0) {
         m_gaussFraction.push_back(1);
      } else {
         double denom = 1. - exp(-2/sig/sig);
         double gauss_int;
         if (psi == 0) {
            gauss_int = 0;
         } else {
// set the object providing the integrand
            s_gfunc = Gint(sig, cr, cp, sp);
// use DGAUS8 to perform the integral
            double err = 1e-5;
            long ierr;
            dgaus8_(&PointSource::gfuncIntegrand, &mup, &mum, 
                    &err, &gauss_int, &ierr);
// // instantiate the integrand object...
//             Gint gfunc(sig, cr, cp, sp);
// // and the integrator object
//             TrapQuad trapQuad(&gfunc);
// //            gauss_int = trapQuad.integral(mup, mum, 10000);
//             gauss_int = trapQuad.integral(mup, mum, 1000);
         }
         if (psi <= roi_radius) {
            double value = ((1. - exp((mum-1.)/sig/sig))
                            + gauss_int/M_PI/sig/sig)/denom;
            m_gaussFraction.push_back(value);
         } else {
            double value =  gauss_int/M_PI/sig/sig/denom;
            m_gaussFraction.push_back(value);
         }
      }
//       std::cout << sig << "  " << m_gaussFraction[i] << std::endl;
   }
}

double PointSource::psfFrac(double energy, double inc) {
   Psf *psf = Psf::instance();
// Compute the fraction of the psf enclosed for this source at this
// energy and inclination
   std::vector<double> psf_params;
   (*psf).fillPsfParams(energy, inc, psf_params);
   double sig1 = psf_params[0]*M_PI/180.;
   double sig2 = psf_params[1]*M_PI/180.;
   double wt = psf_params[2];

// Interpolate the fractions of the Gaussian components contained
// within the region-of-interest
// Compute the index assuming uniform step size in s_sigGauss:
   unsigned int indx = 
      static_cast<int>((sig1 - s_sigGauss[0])
                       /(s_sigGauss[1] - s_sigGauss[0]));
   double frac1 = (sig1 - s_sigGauss[indx])
      /(s_sigGauss[indx+1] - s_sigGauss[indx])
      *(m_gaussFraction[indx+1] - m_gaussFraction[indx]) 
      + m_gaussFraction[indx];

   indx = static_cast<int>((sig2 - s_sigGauss[0])
                           /(s_sigGauss[1] - s_sigGauss[0]));
   double frac2 = (sig2 - s_sigGauss[indx])
      /(s_sigGauss[indx+1] - s_sigGauss[indx])
      *(m_gaussFraction[indx+1] - m_gaussFraction[indx]) 
      + m_gaussFraction[indx];

   return wt*frac1 + (1. - wt)*frac2;
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

void PointSource::makeSigmaVector(int nsig) {
// set up a grid of sigma values in radians
   double sigmin = 0.;
   double sigmax = 50.*M_PI/180;
   double sigstep = (sigmax - sigmin)/(nsig - 1);
   
   s_sigGauss.reserve(nsig);
   for (int i = 0; i < nsig; i++)
      s_sigGauss.push_back(sigstep*i + sigmin);
}

double PointSource::Gint::value(optimizers::Arg &muarg) const {
   double mu = dynamic_cast<optimizers::dArg &>(muarg).getValue();

   double phi;
   if (mu == 1) {
      phi = M_PI;
   } else {
      double arg = (m_cr - mu*m_cp)/sqrt(1. - mu*mu)/m_sp;
      if (arg >= 1.) {
         phi = 0;
      } else if (arg <= -1.) {
         phi = M_PI;
      } else {
         phi = acos(arg);
      }
   }
   double value = phi*exp((mu - 1.)/m_sig/m_sig);
   return value;
}

} // namespace Likelihood
