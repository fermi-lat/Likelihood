/** 
 * @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PointSource.cxx,v 1.107 2009/02/22 20:20:57 jchiang Exp $
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
#include "optimizers/Function.h"

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
                                int eventType,
				CachedResponse* cResp) const {
   const ScData & scData = m_observation->scData();
   astro::SkyDir zAxis = scData.zAxis(time);
   astro::SkyDir xAxis = scData.xAxis(time);
   return fluxDensity(energy, zAxis, xAxis, dir, eventType, cResp);
}

double PointSource::fluxDensity(double energy, const astro::SkyDir &zAxis,
                                const astro::SkyDir &xAxis,
                                const astro::SkyDir &dir,
                                int eventType,
				CachedResponse* cResp) const {
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
      double resp(0);
      if (cResp && cResp->first) {
	 resp = cResp->second;
      } else {
	 resp = respFuncs.totalResponse(energy, energy, zAxis, xAxis,
					m_dir.getDir(), dir, eventType);
	 if (cResp) {
            cResp->first = true;
            cResp->second = resp;
         }
      }
      return spectrum*resp;
   }
}

double PointSource::fluxDensity(double inclination, double phi, double energy, 
                                const astro::SkyDir & appDir, 
                                int evtType,
				CachedResponse* cResp) const {
   optimizers::dArg energy_arg(energy);
   double spectrum = (*m_spectrum)(energy_arg);
   double separation = appDir.difference(getDir())*180./M_PI;
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
   double resp(0);
   if (cResp && cResp->first) {
      resp = cResp->second;
   } else {
      resp = respFuncs.totalResponse(inclination, phi, energy, energy, 
				     separation, evtType);
      if (cResp) {
         cResp->first = true;
         cResp->second = resp;
      }
   }
   return spectrum*resp;
}

double PointSource::fluxDensityDeriv(double energy, double time,
                                     const astro::SkyDir &dir,
                                     int eventType,
                                     const std::string &paramName,
				     CachedResponse* cResp) const {
   const ScData & scData = m_observation->scData();
   astro::SkyDir zAxis = scData.zAxis(time);
   astro::SkyDir xAxis = scData.xAxis(time);
   return fluxDensityDeriv(energy, zAxis, xAxis, dir, eventType, paramName,
			   cResp);
}

double PointSource::fluxDensityDeriv(double energy, 
                                     const astro::SkyDir & zAxis,
                                     const astro::SkyDir & xAxis,
                                     const astro::SkyDir & dir,
                                     int eventType,
                                     const std::string &paramName,
				     CachedResponse* cResp) const {
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude"
   double prefactor;
   if (paramName == "Prefactor" && 
       (prefactor = m_spectrum->getParamValue("Prefactor")) !=0) {
     return fluxDensity(energy, zAxis, xAxis, dir, eventType, cResp)/prefactor;
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

	 double resp(0);
	 if(cResp && cResp->first) {
	   resp = cResp->second;
	 } else {
	    resp = respFuncs.totalResponse(energy, energy, zAxis, xAxis,
					   m_dir.getDir(), dir, eventType);
	    if(cResp)cResp->first=true, cResp->second=resp;
	 }

         return resp*m_spectrum->derivByParam(energy_arg, paramName);
      }
   }
}

double PointSource::
fluxDensityDeriv(double inclination, double phi, double energy,
                 const astro::SkyDir & appDir, int evtType, 
                 const std::string & paramName,
		 CachedResponse* cResp) const {
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
   double separation = appDir.difference(getDir())*180./M_PI;
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude"
   double prefactor;
   if (paramName == "Prefactor" && 
       (prefactor = m_spectrum->getParamValue("Prefactor")) != 0) {
     return fluxDensity(inclination, phi, energy, separation, evtType, cResp)
         /prefactor;
   } else {
/// @todo Implement for finite energy resolution case.
      optimizers::dArg energy_arg(energy);

      double resp(0);
      if(cResp && cResp->first)
	 resp = cResp->second;
      else {
	 resp = respFuncs.totalResponse(inclination, phi, energy, 
					energy, separation, evtType);
	 if(cResp)cResp->first=true, cResp->second=resp;
      }

      return resp*m_spectrum->derivByParam(energy_arg, paramName);
   }
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

double PointSource::fluxDeriv(const std::string & parName) const {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   FluxDeriv my_functor(*m_spectrum, parName);
   TrapQuad fluxIntegral(&my_functor);
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

double PointSource::fluxDeriv(const std::string & parName,
                              double emin, double emax, size_t npts) const {
   std::vector<double> energies;
   energies.reserve(npts);
   double estep(std::log(emax/emin)/float(npts-1));
   for (size_t k=0; k < npts; k++) {
      energies.push_back(emin*std::exp(estep*k));
   }
   FluxDeriv my_functor(*m_spectrum, parName);
   TrapQuad fluxIntegral(&my_functor);
   return fluxIntegral.integral(energies);
}

double PointSource::energyFlux() const {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   EnergyFlux my_functor(*m_spectrum);
   TrapQuad fluxIntegral(&my_functor);
   return fluxIntegral.integral(energies);
}

double PointSource::energyFluxDeriv(const std::string & parName) const {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   EnergyFluxDeriv my_functor(*m_spectrum, parName);
   TrapQuad fluxIntegral(&my_functor);
   return fluxIntegral.integral(energies);
}

double PointSource::energyFlux(double emin, double emax, size_t npts) const {
   std::vector<double> energies;
   energies.reserve(npts);
   double estep(std::log(emax/emin)/float(npts-1));
   for (size_t k=0; k < npts; k++) {
      energies.push_back(emin*std::exp(estep*k));
   }
   EnergyFlux my_functor(*m_spectrum);
   TrapQuad fluxIntegral(&my_functor);
   return fluxIntegral.integral(energies);
}

double PointSource::energyFluxDeriv(const std::string & parName,
                                    double emin, double emax, 
                                    size_t npts) const {
   std::vector<double> energies;
   energies.reserve(npts);
   double estep(std::log(emax/emin)/float(npts-1));
   for (size_t k=0; k < npts; k++) {
      energies.push_back(emin*std::exp(estep*k));
   }
   EnergyFluxDeriv my_functor(*m_spectrum, parName);
   TrapQuad fluxIntegral(&my_functor);
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
                                        "computeExposureWithHyperCube", 4);
   formatter.warn() << "Computing exposure at (" 
                    << srcDir.ra() << ", " 
                    << srcDir.dec() << ")";
   for (std::vector<double>::const_iterator it = energies.begin();
        it != energies.end(); it++) {
      if (verbose) {
         formatter.warn() << ".";
      }
      PointSource::Aeff aeff(*it, srcDir, observation.roiCuts(),
                             observation.respFuncs());
      double exposure_value = observation.expCube().value(srcDir, aeff);
      exposure.push_back(exposure_value);
   }
   formatter.warn() << "!" << std::endl;
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
                                        "computeExposure", 4);
   formatter.warn() << "Computing exposure at (" 
                    << srcDir.ra() << ", " 
                    << srcDir.dec() << ")";
   size_t npts;
   if (roiCuts.maxTime() > scData.vec.back().stoptime) {
      npts = scData.vec.size() - 1;
   } else {
      npts = scData.time_index(roiCuts.maxTime()) + 1;
   }
   for (size_t it = 0; it < npts && it < scData.vec.size(); it++) {
      if (npts/20 > 0 && ((it % (npts/20)) == 0)) {
         formatter.warn() << ".";
      }
      double start(scData.vec.at(it).time);
      double stop(scData.vec.at(it).stoptime);
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
               formatter.warn() << effArea << std::endl;
            }
            exposure[k] += effArea*livetime*fraction;
         }
      }
   }
   formatter.warn() << "!" << std::endl;
}

void PointSource::makeEnergyVector(int nee) {
// A logrithmic grid of true energies for convolving with energy
// dispersion.  Use hard-wired upper and lower energies.
   double trueEmin(10.);
   double trueEmax(6e5);
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
   astro::SkyDir xAxis = scData.xAxis(time);

   PointSource::Aeff aeff(energy, srcDir, roiCuts, respFuncs, time);

   double cos_theta = zAxis()*const_cast<astro::SkyDir&>(srcDir)();

   CLHEP::Hep3Vector yhat(zAxis().cross(xAxis()));
   double phi = 180./M_PI*std::atan2(yhat.dot(srcDir()), 
                                     xAxis().dot(srcDir()));

   double effArea(0);
   effArea = aeff(cos_theta, phi);
   return effArea;
}

PointSource::Aeff::Aeff(double energy, const astro::SkyDir &srcDir,
                        const RoiCuts & roiCuts,
                        const ResponseFunctions & respFuncs, 
                        double time)
   : m_energy(energy), m_srcDir(srcDir), m_respFuncs(respFuncs), m_time(time) {
   
   m_cones.push_back(const_cast<irfInterface::AcceptanceCone *>
                     (&(roiCuts.extractionRegion())));
   m_emin = roiCuts.getEnergyCuts().first;
   m_emax = roiCuts.getEnergyCuts().second;
}

double PointSource::Aeff::operator()(double cos_theta, double phi) const {
   double theta = acos(cos_theta)*180./M_PI;

   double myEffArea = 0;
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt
      = m_respFuncs.begin();

   std::map<int, double> psf_vals;
   for ( ; respIt != m_respFuncs.end(); ++respIt) {
      irfInterface::IPsf *psf = respIt->second->psf();
      irfInterface::IAeff *aeff = respIt->second->aeff();

      double aeff_val = aeff->value(m_energy, theta, phi, m_time);
      if (aeff_val < 0.1) { // kluge.  Psf is likely not well defined out here.
         return 0;
      }

      //
      // Use irfID % 2 to determine if the psf is front(0) or back(1).
      // Assume the IRFs are in class order and use the psf associated
      // with the lowest class number.  This should be backwards
      // compatible with Pass 5 style IRFs where there is only one
      // event class.
      //
      int id(respIt->second->irfID() % 2);

      double psf_val(0);
      std::map<int, double>::const_iterator psf_it(psf_vals.find(id));
      if (psf_it == psf_vals.end()) {
         psf_val = psf->angularIntegral(m_energy, m_srcDir,
                                        theta, phi, m_cones, m_time);
         psf_vals[id] = psf_val;
      } else {
         psf_val = psf_it->second;
      }

      if (m_respFuncs.useEdisp()) {
         irfInterface::IEdisp *edisp = respIt->second->edisp();
         double edisp_val = edisp->integral(m_emin, m_emax, m_energy, 
                                            theta, phi, m_time);
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
