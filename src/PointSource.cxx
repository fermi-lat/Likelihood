/** 
 * @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PointSource.cxx,v 1.87 2006/09/19 23:07:32 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

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

namespace Likelihood {

std::vector<double> PointSource::s_trueEnergies(0);

PointSource::PointSource(const Observation * observation) 
   : Source(observation) {
   setDir(0., 0., false);
   m_srcType = "Point";
   if (s_trueEnergies.empty()) {
      makeEnergyVector();
   }
}

PointSource::PointSource(double ra, double dec,
                         const Observation & observation, bool verbose) 
   : Source(&observation) {
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

PointSource::~PointSource() {
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
      double resp(respFuncs.totalResponse(energy, energy, zAxis, xAxis,
                                          m_dir.getDir(), dir, eventType));
      return spectrum*resp;
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
   double prefactor;
   if (paramName == "Prefactor" && 
       (prefactor = m_spectrum->getParamValue("Prefactor")) !=0) {
      return fluxDensity(energy, zAxis, xAxis, dir, eventType)/prefactor;
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
   double prefactor;
   if (paramName == "Prefactor" && 
       (prefactor = m_spectrum->getParamValue("Prefactor")) != 0) {
      return fluxDensity(inclination, phi, energy, separation, evtType)
         /prefactor;
   } else {
/// @todo Implement for finite energy resolution case.
      optimizers::dArg energy_arg(energy);
      return respFuncs.totalResponse(inclination, phi, energy, 
                                     energy, separation, evtType)
         *m_spectrum->derivByParam(energy_arg, paramName);
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
   if (verbose) {
      st_stream::StreamFormatter formatter("PointSource",
                                           "computeExposure", 2);
      for (unsigned int i = 0; i < energies.size(); i++) {
         formatter.info() << energies.at(i) << "  " 
                          << m_exposure.at(i) << std::endl;
      }
   }
}

double PointSource::flux() const {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   TrapQuad fluxIntegral(m_spectrum);
   return fluxIntegral.integral(energies);
}

double PointSource::flux(double emin, double emax, size_t npts) const {
   std::vector<double> energies;
   energies.reserve(npts);
   double estep(std::log(emax/emin)/float(npts-1));
   for (size_t k=0; k < npts; k++) {
      energies.push_back(emin*std::exp(estep*k));
   }
   TrapQuad fluxIntegral(m_spectrum);
   return fluxIntegral.integral(energies);
}

void PointSource::
computeExposureWithHyperCube(const astro::SkyDir & srcDir,
                             const std::vector<double> & energies, 
                             const Observation & observation,
                             std::vector<double> & exposure, 
                             bool verbose) {
   (void)(verbose);
   exposure.clear();

   st_stream::StreamFormatter formatter("PointSource",
                                        "computeExposureWithHyperCube", 3);
   formatter.info() << "Computing exposure at (" 
                    << srcDir.ra() << ", " 
                    << srcDir.dec() << ")";
   for (std::vector<double>::const_iterator it = energies.begin();
        it != energies.end(); it++) {
      if (verbose) {
         formatter.info() << ".";
      }
      PointSource::Aeff aeff(*it, srcDir, observation.roiCuts(),
                             observation.respFuncs());
      double exposure_value = observation.expCube().value(srcDir, aeff);
      exposure.push_back(exposure_value);
   }
   formatter.info() << "!" << std::endl;
}

void PointSource::computeExposure(const astro::SkyDir & srcDir,
                                  const std::vector<double> &energies,
                                  const Observation & observation,
                                  std::vector<double> &exposure,
                                  bool verbose) {
   (void)(verbose);
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

   st_stream::StreamFormatter formatter("PointSource",
                                        "computeExposure", 3);
   formatter.info() << "Computing exposure at (" 
                    << srcDir.ra() << ", " 
                    << srcDir.dec() << ")";
   size_t npts;
   if (roiCuts.maxTime() > scData.vec.back().time) {
      npts = scData.vec.size() - 1;
   } else {
      npts = scData.time_index(roiCuts.maxTime()) + 1;
   }
   for (unsigned int it = 0; it < npts && it < scData.vec.size()-1; it++) {
      if (npts/20 > 0 && ((it % (npts/20)) == 0)) {
         formatter.info() << ".";
      }
      double start(scData.vec.at(it).time);
      double stop(scData.vec.at(it+1).time);
      double livetime(scData.vec.at(it).livetime);
      double fraction(0);

      std::vector< std::pair<double, double> > timeRanges;
      std::vector< std::pair<double, double> > gtis;
      roiCuts.getTimeCuts(timeRanges);
      roiCuts.getGtis(gtis);
      bool includeInterval = 
         LikeExposure::acceptInterval(start, stop, timeRanges, gtis, fraction);

// Compute the inclination and check if it's within response matrix
// cut-off angle
      double inc = srcDir.difference(scData.vec.at(it).zAxis)*180/M_PI;
      if (inc > 90.) {
         includeInterval = false;
      }

// Having checked for relevant constraints, add the exposure
// contribution for each energy
      if (includeInterval) {
         for (unsigned int k = 0; k < energies.size(); k++) {
            double time = (start + stop)/2.;
            double effArea = sourceEffArea(srcDir, energies[k], time, 
                                           scData, roiCuts, respFuncs);
            if (effArea < 0 || fraction < 0 || (stop-start) < 0) {
               formatter.info() << effArea << std::endl;
            }
            exposure[k] += effArea*livetime*fraction;
         }
      }
   }
   formatter.info() << "!" << std::endl;
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
   
   st_stream::StreamFormatter formatter("PointSource", "sourceEffArea", 3);

   double effArea(0);
   try {
      effArea = aeff(cos_theta);
   } catch (std::exception & eObj) {
      formatter.info() << eObj.what() << "\n"
                       << "cos_theta = " << cos_theta
                       << std::endl;
   } catch (...) {
      formatter.info() << "caught unknown exception for "
                       << "cos_theta = " << cos_theta
                       << std::endl;
   }
   return effArea;
}

PointSource::Aeff::Aeff(double energy, const astro::SkyDir &srcDir,
                        const RoiCuts & roiCuts,
                        const ResponseFunctions & respFuncs)
   : m_energy(energy), m_srcDir(srcDir), m_respFuncs(respFuncs) {
   
   m_cones.push_back(const_cast<irfInterface::AcceptanceCone *>
                     (&(roiCuts.extractionRegion())));
   m_emin = roiCuts.getEnergyCuts().first;
   m_emax = roiCuts.getEnergyCuts().second;
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

      double aeff_val = aeff->value(m_energy, theta, phi);
      if (aeff_val < 0.1) { // kluge.  Psf is likely not well defined out here.
         return 0;
      }
// Sledgehammer approach to handle IRF misbehavior, probably not needed given
// the aeff_val < 0.1 test.
      double psf_val(0);
      try {
         psf_val = psf->angularIntegral(m_energy, m_srcDir,
                                        theta, phi, m_cones);
      } catch (std::exception & eObj) { 
         st_stream::StreamFormatter formatter("PointSource::Aeff", 
                                              "operator()", 2);
         formatter.info() << eObj.what() << std::endl;
      }
      if (m_respFuncs.useEdisp()) {
         irfInterface::IEdisp *edisp = respIt->second->edisp();
         double edisp_val = edisp->integral(m_emin, m_emax, m_energy, 
                                            theta, phi);
         myEffArea += psf_val*aeff_val*edisp_val;
      } else {
         myEffArea += psf_val*aeff_val;
      }
   }
   if (myEffArea < 0) {
      return 0;
   }
   return myEffArea;
}

} // namespace Likelihood
