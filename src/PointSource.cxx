/** 
 * @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/PointSource.cxx,v 1.119 2013/01/09 00:44:41 jchiang Exp $
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

#define DEBUGASSERT(x) assert(x)

PointSource::PointSource(const Observation * observation) 
   : Source(observation) {
   setDir(0., 0., false);
   m_srcType = "Point";
}

PointSource::PointSource(double ra, double dec,
                         const Observation & observation, bool verbose) 
   : Source(&observation) {
   setDir(ra, dec, true, verbose);
   m_srcType = "Point";
}

PointSource::PointSource(const PointSource & rhs) 
   : Source(rhs), m_dir(rhs.m_dir) {
   m_functions["Position"] = &m_dir;
   m_functions["Spectrum"] = m_spectrum;
}

PointSource::~PointSource() {
}

void PointSource::computeResponse(Response& resp, const Event & evt) const {
   const double& energy         = evt.getEnergy();
   const astro::SkyDir & zAxis  = evt.zAxis();
   const astro::SkyDir & xAxis  = evt.xAxis();
   const astro::SkyDir & dir    = evt.getDir();
   int eventType                = evt.getType();
   double time                  = evt.getArrTime();
   const ResponseFunctions & respFuncs = m_observation->respFuncs();
   if(m_edisp() == ED_NONE) {
      resp.resize(1);
      resp[0] = respFuncs.totalResponse(energy, energy, zAxis, xAxis,
					m_dir.getDir(), dir, eventType, time,
					/* useEdisp = */ false);
   } else {
      const std::vector<double> & true_energies = evt.trueEnergies();
      resp.resize(true_energies.size());
      for(unsigned ienergy = 0; ienergy<resp.size(); ienergy++) {
	 resp[ienergy] = respFuncs.totalResponse(true_energies[ienergy],
						 energy,
						 zAxis, xAxis, m_dir.getDir(),
						 dir, eventType, time,
						 /* useEdisp = */ true);
      }
      if(m_edisp() == ED_GAUSSIAN) {
	 computeGaussianParams(resp, resp, true_energies);
      }
   }
}

void PointSource::computeExposure(const std::vector<double> & energies,
                                  bool verbose) {
   m_energies = energies;
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
   const std::vector<double> & energies(m_energies);
   bool useLog;
   TrapQuad fluxIntegral(m_spectrum, useLog=true);
   return fluxIntegral.integral(energies);
}

double PointSource::fluxDeriv(const std::string & parName) const {
   const std::vector<double> & energies(m_energies);
   FluxDeriv my_functor(*m_spectrum, parName);
   bool useLog;
   TrapQuad fluxIntegral(&my_functor, useLog=true);
   return fluxIntegral.integral(energies);
}

double PointSource::flux(double emin, double emax, size_t npts) const {
   std::vector<double> energies;
   energies.reserve(npts);
   double estep(std::log(emax/emin)/float(npts-1));
   for (size_t k=0; k < npts; k++) {
      energies.push_back(emin*std::exp(estep*k));
   }
   bool useLog;
   TrapQuad fluxIntegral(m_spectrum, useLog=true);
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
   bool useLog;
   TrapQuad fluxIntegral(&my_functor, useLog=true);
   return fluxIntegral.integral(energies);
}

double PointSource::energyFlux() const {
   const std::vector<double> & energies(m_energies);
   EnergyFlux my_functor(*m_spectrum);
   bool useLog;
   TrapQuad fluxIntegral(&my_functor, useLog=true);
   return fluxIntegral.integral(energies);
}

double PointSource::energyFluxDeriv(const std::string & parName) const {
   const std::vector<double> & energies(m_energies);
   EnergyFluxDeriv my_functor(*m_spectrum, parName);
   bool useLog;
   TrapQuad fluxIntegral(&my_functor, useLog=true);
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
   bool useLog;
   TrapQuad fluxIntegral(&my_functor, useLog=true);
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
   bool useLog;
   TrapQuad fluxIntegral(&my_functor, useLog=true);
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
      double time((observation.expCube().tstart() 
                   + observation.expCube().tstop())/2.);
      PointSource::Aeff aeff(*it, srcDir, observation.roiCuts(),
                             observation.respFuncs(), time,
                             observation.expCube().hasPhiDependence());
      double exposure_value = observation.expCube().value(srcDir, aeff, *it);
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
   if (scData.numIntervals() == 0) {
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
   size_t npts(scData.numIntervals() - 1);
   if (roiCuts.maxTime() <= scData.stop(npts)) {
      npts = scData.time_index(roiCuts.maxTime()) + 1;
   }

   for (size_t it = 0; it < npts && it < scData.numIntervals(); it++) {
      if (npts/20 > 0 && ((it % (npts/20)) == 0)) {
         formatter.warn() << ".";
      }
      double start(scData.start(it));
      double stop(scData.stop(it));
      double livetime(scData.livetime(it));
      double fraction(0);
      double ltfrac(livetime/(stop - start));

      std::vector< std::pair<double, double> > timeRanges;
      std::vector< std::pair<double, double> > gtis;
      roiCuts.getTimeCuts(timeRanges);
      roiCuts.getGtis(gtis);
      bool includeInterval = 
         LikeExposure::acceptInterval(start, stop, timeRanges, gtis, fraction);

// Compute the inclination and check if it's within response matrix
// cut-off angle
      double inc = srcDir.difference(scData.zAxis(it))*180/M_PI;
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
            const irfInterface::IEfficiencyFactor * efficiency_factor
               = respFuncs.efficiencyFactor();
            double efficiency(1);
            if (efficiency_factor) {
               efficiency = efficiency_factor->value(energies[k], ltfrac);
            }
            exposure[k] += effArea*livetime*fraction*efficiency;
         }
      }
   }
   formatter.warn() << "!" << std::endl;
}

#if 0
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
#endif

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
                        double time, bool usePhiDependence)
   : ExposureCube::AeffBase(), m_energy(energy), m_srcDir(srcDir),
     m_respFuncs(respFuncs), m_time(time),
     m_usePhiDependence(usePhiDependence) {
   m_cones.push_back(const_cast<irfInterface::AcceptanceCone *>
                     (&(roiCuts.extractionRegion())));
//    m_emin = roiCuts.getEnergyCuts().first;
//    m_emax = roiCuts.getEnergyCuts().second;
}

double PointSource::Aeff::value(double cos_theta, double phi) const {
   double theta = acos(cos_theta)*180./M_PI;

   double myEffArea = 0;
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt
      = m_respFuncs.begin();

   std::map<int, double> psf_vals;
   for ( ; respIt != m_respFuncs.end(); ++respIt) {
      irfInterface::IPsf *psf = respIt->second->psf();
      irfInterface::IAeff *aeff = respIt->second->aeff();

      bool savedPhiDepState(aeff->usePhiDependence());
      aeff->setPhiDependence(m_usePhiDependence);
      double aeff_val = aeff->value(m_energy, theta, phi, m_time);
      aeff->setPhiDependence(savedPhiDepState);
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
//          irfInterface::IEdisp *edisp = respIt->second->edisp();
//          double edisp_val = edisp->integral(m_emin, m_emax, m_energy, 
//                                             theta, phi, m_time);
//          myEffArea += psf_val*aeff_val*edisp_val;
         throw std::runtime_error("energy dispersion handling disabled");
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
