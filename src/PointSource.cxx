/** 
 * @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PointSource.cxx,v 1.66 2005/03/18 01:05:17 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "optimizers/dArg.h"

#include "irfInterface/Irfs.h"

#include "Likelihood/LikeExposure.h"
#include "Likelihood/Observation.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/TrapQuad.h"

#include "Verbosity.h"

namespace Likelihood {

std::vector<double> PointSource::s_trueEnergies(0);

PointSource::PointSource() : m_spectrum(0), m_observation(0) {
   setDir(0., 0., false);
   m_srcType = "Point";
   if (s_trueEnergies.empty()) {
      makeEnergyVector();
   }
}

PointSource::
PointSource(double ra, double dec, const Observation & observation,
            bool verbose) : m_spectrum(0), m_observation(&observation) {
   setDir(ra, dec, true, verbose);
   m_srcType = "Point";
   if (s_trueEnergies.empty()) {
      makeEnergyVector();
   }
}

PointSource::PointSource(const PointSource &rhs) : Source(rhs) {
   m_dir = rhs.m_dir;
   m_functions["Position"] = &m_dir;

   m_spectrum = rhs.m_spectrum->clone();
   m_functions["Spectrum"] = m_spectrum;

   m_exposure = rhs.m_exposure;
   m_observation = rhs.m_observation;
}

double PointSource::fluxDensity(double energy, double time,
                                const astro::SkyDir &dir,
                                int eventType) const {
   const ScData & scData = m_observation->scData();
   astro::SkyDir zAxis = scData.zAxis(time);
   astro::SkyDir xAxis = scData.xAxis(time);
   return fluxDensity(energy, zAxis, xAxis, dir, eventType);
}

double PointSource::fluxDensity(double energy, const astro::SkyDir &zAxis,
                                const astro::SkyDir &xAxis,
                                const astro::SkyDir &dir,
                                int eventType) const {
// Scale the energy spectrum by the psf value and the effective area
// and convolve with the energy dispersion, if appropriate, all of
// which are functions of time and spacecraft attitude and orbital
// position.
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
   if (respFuncs.useEdisp()) {
      unsigned int npts(s_trueEnergies.size());
      std::vector<double> my_integrand(npts);
      for (unsigned int k = 0; k < npts; k++) {
         optimizers::dArg energy_arg(s_trueEnergies[k]);
         double spectrum = (*m_spectrum)(energy_arg);
         my_integrand[k] = spectrum
            *respFuncs.totalResponse(s_trueEnergies[k], energy,
                                     zAxis, xAxis, m_dir.getDir(), 
                                     dir, eventType);
      }
      TrapQuad trapQuad(s_trueEnergies, my_integrand);
      return trapQuad.integral();
   } else {
      optimizers::dArg energy_arg(energy);
      double spectrum = (*m_spectrum)(energy_arg);
      return spectrum*respFuncs.totalResponse(energy, energy, zAxis, xAxis,
                                              m_dir.getDir(), dir, eventType);
   }
}

double PointSource::fluxDensity(double inclination, double phi, double energy, 
                                const astro::SkyDir & appDir, 
                                int evtType) const {
/// @todo add implementation for energy dispersion.
   optimizers::dArg energy_arg(energy);
   double spectrum = (*m_spectrum)(energy_arg);
   double separation = appDir.difference(getDir())*180./M_PI;
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
   return spectrum*respFuncs.totalResponse(inclination, phi, energy, energy, 
                                           separation, evtType);
}

double PointSource::fluxDensityDeriv(double energy, double time,
                                     const astro::SkyDir &dir,
                                     int eventType,
                                     const std::string &paramName) const {
   const ScData & scData = m_observation->scData();
   astro::SkyDir zAxis = scData.zAxis(time);
   astro::SkyDir xAxis = scData.xAxis(time);
   return fluxDensityDeriv(energy, zAxis, xAxis, dir, eventType, paramName);
}

double PointSource::fluxDensityDeriv(double energy, 
                                     const astro::SkyDir & zAxis,
                                     const astro::SkyDir & xAxis,
                                     const astro::SkyDir & dir,
                                     int eventType,
                                     const std::string &paramName) const {
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude"
   if (paramName == "Prefactor") {
      return fluxDensity(energy, zAxis, xAxis, dir, eventType)
         /m_spectrum->getParamValue("Prefactor");
   } else {
      if (respFuncs.useEdisp()) {
         unsigned int npts(s_trueEnergies.size());
         std::vector<double> my_integrand(npts);
         for (unsigned int k = 0; k < npts; k++) {
            optimizers::dArg energy_arg(s_trueEnergies[k]);
            my_integrand[k] = m_spectrum->derivByParam(energy_arg, paramName)
               *respFuncs.totalResponse(s_trueEnergies[k], energy,
                                        zAxis, xAxis, m_dir.getDir(), 
                                        dir, eventType);
         }
         TrapQuad trapQuad(s_trueEnergies, my_integrand);
         return trapQuad.integral();
      } else {
         optimizers::dArg energy_arg(energy);
         return respFuncs.totalResponse(energy, energy, zAxis, xAxis,
                                        m_dir.getDir(), dir, 
                                        eventType)
            *m_spectrum->derivByParam(energy_arg, paramName);
      }
   }
}

double PointSource::
fluxDensityDeriv(double inclination, double phi, double energy,
                 const astro::SkyDir & appDir, int evtType, 
                 const std::string & paramName) const {
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
   double separation = appDir.difference(getDir())*180./M_PI;
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude"
   if (paramName == "Prefactor") {
      return fluxDensity(inclination, phi, energy, separation, evtType)
         /m_spectrum->getParamValue("Prefactor");
   } else {
/// @todo Implement for finite energy resolution case.
      optimizers::dArg energy_arg(energy);
      return respFuncs.totalResponse(inclination, phi, energy, 
                                     energy, separation, evtType)
         *m_spectrum->derivByParam(energy_arg, paramName);
   }
}

double PointSource::Npred() {
   optimizers::Function *specFunc = m_functions["Spectrum"];

// Evaluate the Npred integrand at the abscissa points contained in
// RoiCuts::energies().
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   std::vector<double> NpredIntegrand(energies.size());
   for (unsigned int k = 0; k < energies.size(); k++) {
      optimizers::dArg eArg(energies[k]);
      NpredIntegrand[k] = (*specFunc)(eArg)*m_exposure[k];
   }
   TrapQuad trapQuad(energies, NpredIntegrand);
   return trapQuad.integral();
}

double PointSource::Npred(double emin, double emax) {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   if (emin < energies.front() || emax > energies.back()) {
      throw std::out_of_range("PointSource::Npred(emin, emax)");
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

double PointSource::NpredDeriv(const std::string &paramName) {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   optimizers::Function *specFunc = m_functions["Spectrum"];

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

double PointSource::pixelCounts(double emin, double emax,
                                double wtMin, double wtMax) const {
   optimizers::Function & spectrum = *m_spectrum;
   optimizers::dArg eminArg(emin);
   optimizers::dArg emaxArg(emax);
   return (spectrum(emaxArg)*wtMax + spectrum(eminArg)*wtMin)*(emax - emin)/2.;
}

double PointSource::pixelCountsDeriv(double emin, double emax,
                                     double wtMin, double wtMax,
                                     const std::string & paramName) const {
   optimizers::Function & spectrum = *m_spectrum;
   optimizers::dArg eminArg(emin);
   optimizers::dArg emaxArg(emax);
   return (spectrum.derivByParam(emaxArg, paramName)*wtMax +
           spectrum.derivByParam(eminArg, paramName)*wtMin)*(emax - emin)/2.;
}

void PointSource::computeExposure(bool verbose) {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   astro::SkyDir srcDir = getDir();
   if (m_observation->expCube().haveFile()) {
      computeExposureWithHyperCube(srcDir, energies, *m_observation,
                                   m_exposure, verbose);
   } else {
// Exposure time hypercube is not available, so perform sums using
// ScData.
      computeExposure(srcDir, energies, *m_observation,
                      m_exposure, verbose);
   }
   if (print_output() && verbose) {
      for (unsigned int i = 0; i < energies.size(); i++) {
         std::cout << energies.at(i) << "  " << m_exposure.at(i) << std::endl;
      }
   }
}

void PointSource::
computeExposureWithHyperCube(const astro::SkyDir & srcDir,
                             const std::vector<double> & energies, 
                             const Observation & observation,
                             std::vector<double> & exposure, 
                             bool verbose) {
   exposure.clear();

   if (print_output() && verbose) {
      std::cerr << "Computing exposure at (" 
                << srcDir.ra() << ", " 
                << srcDir.dec() << ")";
   }
   for (std::vector<double>::const_iterator it = energies.begin();
        it != energies.end(); it++) {
      if (print_output() && verbose) std::cerr << ".";
      PointSource::Aeff aeff(*it, srcDir, observation.roiCuts(),
                             observation.respFuncs());
      double exposure_value = observation.expCube().value(srcDir, aeff);
      exposure.push_back(exposure_value);
   }
   if (print_output() && verbose) std::cerr << "!" << std::endl;
}

void PointSource::computeExposure(const astro::SkyDir & srcDir,
                                  const std::vector<double> &energies,
                                  const Observation & observation,
                                  std::vector<double> &exposure,
                                  bool verbose) {
   const ScData & scData = observation.scData();
   const RoiCuts & roiCuts = observation.roiCuts();
   const ResponseFunctions & respFuncs = observation.respFuncs();

// Don't compute anything if there is no ScData.
   if (scData.vec.size() == 0) {
      return;
   }

// Initialize the exposure vector with zeros
   exposure.clear();
   exposure.resize(energies.size());

   if (print_output() && verbose) {
      std::cerr << "Computing exposure at (" 
                << srcDir.ra() << ", " 
                << srcDir.dec() << ")";
   }
//   unsigned int npts = scData.vec.size() - 1;
   unsigned int npts = scData.time_index(roiCuts.maxTime()) + 1;
   for (unsigned int it = 0; it < npts && it < scData.vec.size(); it++) {
      if (print_output() && 
          npts/20 > 0 && ((it % (npts/20)) == 0) && verbose) std::cerr << ".";
      double start(scData.vec[it].time);
      double stop(scData.vec[it+1].time);
      double fraction(0);

      bool includeInterval = 
         LikeExposure::acceptInterval(start, stop, roiCuts.timeRangeCuts(),
                                      roiCuts.gtis(), fraction);

// Compute the inclination and check if it's within response matrix
// cut-off angle
      double inc = srcDir.difference(scData.vec[it].zAxis)*180/M_PI;
      if (inc > 90.) {
         includeInterval = false;
      }

// Having checked for relevant constraints, add the exposure
// contribution for each energy
      if (includeInterval) {
         for (unsigned int k = 0; k < energies.size(); k++) {
            double time = (start + stop)/2.;
            exposure[k] += sourceEffArea(srcDir, energies[k], time, 
                                         scData, roiCuts, respFuncs)
               *(stop - start)*fraction;
         }
      }
   }
   if (print_output() && verbose) std::cerr << "!" << std::endl;
}

void PointSource::makeEnergyVector(int nee) {
// A logrithmic grid of true energies for convolving with energy
// dispersion.  Use hard-wired upper and lower energies.
   double trueEmin(18.);
   double trueEmax(3.17e5);
   double trueEstep = log(trueEmax/trueEmin)/(nee-1.);
   s_trueEnergies.clear();
   s_trueEnergies.reserve(nee);
   for (int i = 0; i < nee; i++) {
      s_trueEnergies.push_back(trueEmin*exp(i*trueEstep));
   }
}

double PointSource::sourceEffArea(const astro::SkyDir & srcDir, 
                                  double energy, double time,
                                  const ScData & scData,
                                  const RoiCuts & roiCuts,
                                  const ResponseFunctions & respFuncs) {

   astro::SkyDir zAxis = scData.zAxis(time);
//   astro::SkyDir xAxis = scData.xAxis(time);

   PointSource::Aeff aeff(energy, srcDir, roiCuts, respFuncs);

   double cos_theta = zAxis()*const_cast<astro::SkyDir&>(srcDir)();

   return aeff(cos_theta);
}
   
std::vector<irfInterface::AcceptanceCone *> PointSource::Aeff::s_cones;
double PointSource::Aeff::s_emin;
double PointSource::Aeff::s_emax;

PointSource::Aeff::Aeff(double energy, const astro::SkyDir &srcDir,
                        const RoiCuts & roiCuts,
                        const ResponseFunctions & respFuncs)
   : m_energy(energy), m_srcDir(srcDir), m_respFuncs(respFuncs) {
   
   if (s_cones.size() == 0) {
      s_cones.push_back(const_cast<irfInterface::AcceptanceCone *>
                        (&(roiCuts.extractionRegion())));
      s_emin = roiCuts.getEnergyCuts().first;
      s_emax = roiCuts.getEnergyCuts().second;
   }
}

double PointSource::Aeff::operator()(double cos_theta) const {
   double theta = acos(cos_theta)*180./M_PI;
   static double phi;

   double myEffArea = 0;
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt
      = m_respFuncs.begin();

   for ( ; respIt != m_respFuncs.end(); respIt++) {
      irfInterface::IPsf *psf = respIt->second->psf();
      irfInterface::IAeff *aeff = respIt->second->aeff();

      double psf_val = psf->angularIntegral(m_energy, m_srcDir,
                                            theta, phi, s_cones);
      double aeff_val = aeff->value(m_energy, theta, phi);

      if (m_respFuncs.useEdisp()) {
         irfInterface::IEdisp *edisp = respIt->second->edisp();
         double edisp_val = edisp->integral(s_emin, s_emax, m_energy, 
                                            theta, phi);
         myEffArea += psf_val*aeff_val*edisp_val;
      } else {
         myEffArea += psf_val*aeff_val;
      }
   }
   return myEffArea;
}

} // namespace Likelihood
