/** 
 * @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PointSource.cxx,v 1.46 2004/08/03 21:28:46 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <string>
#include <vector>

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "optimizers/dArg.h"

#include "irfInterface/Irfs.h"

#include "map_tools/Exposure.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/ScData.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/TrapQuad.h"

namespace Likelihood {

map_tools::Exposure * PointSource::s_exposure = 0;

bool PointSource::s_haveStaticMembers = false;
std::vector<double> PointSource::s_energies;
std::vector<double> PointSource::s_trueEnergies;

PointSource::PointSource(const PointSource &rhs) : Source(rhs) {
// Make a deep copy.
   m_dir = rhs.m_dir;
   m_functions["Position"] = &m_dir;

   m_spectrum = rhs.m_spectrum->clone();
   m_functions["Spectrum"] = m_spectrum;

   m_exposure = rhs.m_exposure;
   m_srcType = rhs.m_srcType;
}

void PointSource::readExposureCube(std::string expCubeFile) {
   facilities::Util::expandEnvVar(&expCubeFile);
   s_exposure = new map_tools::Exposure(expCubeFile);
}

double PointSource::fluxDensity(double energy, double time,
                                const astro::SkyDir &dir,
                                int eventType) const {
// Scale the energy spectrum by the psf value and the effective area
// and convolve with the energy dispersion, if appropriate, all of
// which are functions of time and spacecraft attitude and orbital
// position.
   if (ResponseFunctions::useEdisp()) {
      unsigned int npts(s_trueEnergies.size());
      std::vector<double> my_integrand(npts);
      for (unsigned int k = 0; k < npts; k++) {
         optimizers::dArg energy_arg(s_trueEnergies[k]);
         double spectrum = (*m_spectrum)(energy_arg);
         my_integrand[k] = spectrum
            *ResponseFunctions::totalResponse(time, s_trueEnergies[k], energy,
                                              m_dir.getDir(), dir, eventType);
      }
      TrapQuad trapQuad(s_trueEnergies, my_integrand);
      return trapQuad.integral();
   } else {
      optimizers::dArg energy_arg(energy);
      double spectrum = (*m_spectrum)(energy_arg);
      return ResponseFunctions::totalResponse(time, energy, energy, 
                                              m_dir.getDir(), dir, eventType)
         *spectrum;
   }
}

double PointSource::fluxDensityDeriv(double energy, double time,
                                     const astro::SkyDir &dir,
                                     int eventType,
                                     const std::string &paramName) const {
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude"
   if (paramName == "Prefactor") {
      return fluxDensity(energy, time, dir, eventType)
         /m_spectrum->getParamValue("Prefactor");
   } else {
      if (ResponseFunctions::useEdisp()) {
         unsigned int npts(s_trueEnergies.size());
         std::vector<double> my_integrand(npts);
         for (unsigned int k = 0; k < npts; k++) {
            optimizers::dArg energy_arg(s_trueEnergies[k]);
            my_integrand[k] = m_spectrum->derivByParam(energy_arg, paramName)
               *ResponseFunctions::totalResponse(time, s_trueEnergies[k], 
                                                 energy, m_dir.getDir(), 
                                                 dir, eventType);
         }
         TrapQuad trapQuad(s_trueEnergies, my_integrand);
         return trapQuad.integral();
      } else {
         optimizers::dArg energy_arg(energy);
         return ResponseFunctions::totalResponse(time, energy, energy,
                                                 m_dir.getDir(), dir, 
                                                 eventType)
            *m_spectrum->derivByParam(energy_arg, paramName);
      }
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

void PointSource::computeExposure(bool verbose) {
   if (s_exposure == 0) {
// Exposure time hypercube is not available, so perform sums using
// ScData.
      computeExposure(s_energies, m_exposure, verbose);
   } else {
      computeExposureWithHyperCube(s_energies, m_exposure, verbose);
   }
   if (verbose) {
      for (unsigned int i = 0; i < s_energies.size(); i++) {
         std::cout << s_energies[i] << "  " << m_exposure[i] << std::endl;
      }
   }
}

void PointSource::computeExposureWithHyperCube(std::vector<double> &energies, 
                                               std::vector<double> &exposure, 
                                               bool verbose) {
   if (!s_haveStaticMembers 
       || RoiCuts::instance()->getEnergyCuts().first != s_energies.front()
       || RoiCuts::instance()->getEnergyCuts().second != s_energies.back()) {
      makeEnergyVector();
      s_haveStaticMembers = true;
   }
   exposure.clear();

   astro::SkyDir srcDir = getDir();
   if (verbose) {
      std::cerr << "Computing exposure at (" 
                << srcDir.ra() << ", " 
                << srcDir.dec() << ")";
   }
   for (std::vector<double>::const_iterator it = energies.begin();
        it != energies.end(); it++) {
      if (verbose) std::cerr << ".";
      PointSource::Aeff aeff(*it, srcDir);
      double exposure_value = (*s_exposure)(srcDir, aeff);
      exposure.push_back(exposure_value);
   }
   if (verbose) std::cerr << "!" << std::endl;
}

void PointSource::computeExposure(std::vector<double> &energies,
                                  std::vector<double> &exposure,
                                  bool verbose) {
   ScData *scData = ScData::instance();
   RoiCuts *roiCuts = RoiCuts::instance();

// Don't compute anything if there is no ScData.
   if (scData->vec.size() == 0) return;

   if (!s_haveStaticMembers 
       || RoiCuts::instance()->getEnergyCuts().first != s_energies.front()
       || RoiCuts::instance()->getEnergyCuts().second != s_energies.back()) {
      makeEnergyVector();
      s_haveStaticMembers = true;
   }

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
      std::pair<double, double> thisInterval;
      thisInterval.first = scData->vec[it].time;
      thisInterval.second = scData->vec[it+1].time;

// Check if this interval passes the time cuts
      std::vector< std::pair<double, double> > timeCuts;
      roiCuts->getTimeCuts(timeCuts);

      for (unsigned int itcut = 0; itcut < timeCuts.size(); itcut++) {
         if ( !(includeInterval 
                = overlapInterval(timeCuts[itcut], thisInterval)) )
            break;
      }

// Check for SAA passage
      if (scData->vec[it].inSaa) includeInterval = false;

// Compute the inclination and check if it's within response matrix
// cut-off angle
      double inc = getSeparation(scData->vec[it].zAxis)*180/M_PI;
      if (inc > 90.) includeInterval = false;

// Having checked for relevant constraints, add the exposure
// contribution for each energy
      if (includeInterval) {
         for (unsigned int k = 0; k < energies.size(); k++) {
            double time = (thisInterval.second + thisInterval.first)/2.;
            exposure[k] += sourceEffArea(energies[k], time)
               *(thisInterval.second - thisInterval.first);
//             exposure[k] += 
//                ::sourceEffArea(energies[k], time, m_dir.getDir())
//                *(thisInterval.second - thisInterval.first);
         }
      }
   }
   if (verbose) std::cerr << "!" << std::endl;
}

void PointSource::makeEnergyVector(int nee) {
// A logrithmic grid of true energies for convolving with energy
// dispersion.  Use hard-wired upper and lower energies.
   int npts(200);
   double trueEmin(18.);
   double trueEmax(3.17e5);
   double trueEstep = log(trueEmax/trueEmin)/(npts-1.);
   s_trueEnergies.clear();
   s_trueEnergies.reserve(npts);
   for (int i = 0; i < npts; i++) {
      s_trueEnergies.push_back(trueEmin*exp(i*trueEstep));
   }

   if (ResponseFunctions::useEdisp()) {
      s_energies = s_trueEnergies;
   } else {
      RoiCuts *roiCuts = RoiCuts::instance();
      
// set up a logrithmic grid of energies for doing the integral over 
// the spectrum
      double emin = (roiCuts->getEnergyCuts()).first;
      double emax = (roiCuts->getEnergyCuts()).second;
      double estep = log(emax/emin)/(nee-1);
   
      s_energies.clear();
      s_energies.reserve(nee);
      for (int i = 0; i < nee; i++) {
         s_energies.push_back(emin*exp(i*estep));
      }
   }
}

double PointSource::sourceEffArea(double energy, double time) const {
   Likelihood::ScData * scData = Likelihood::ScData::instance();

   astro::SkyDir zAxis = scData->zAxis(time);
//   astro::SkyDir xAxis = scData->xAxis(time);

   const astro::SkyDir & srcDir = m_dir.getDir();

   Likelihood::PointSource::Aeff aeff(energy, srcDir);

   double cos_theta = zAxis()*const_cast<astro::SkyDir&>(srcDir)();

   return aeff(cos_theta);
}
   
std::vector<irfInterface::AcceptanceCone *> PointSource::Aeff::s_cones;
double PointSource::Aeff::s_emin;
double PointSource::Aeff::s_emax;

PointSource::Aeff::Aeff(double energy, const astro::SkyDir &srcDir)
   : m_energy(energy), m_srcDir(srcDir) {
   
   if (s_cones.size() == 0) {
      RoiCuts * roiCuts = RoiCuts::instance();
      s_cones.push_back(const_cast<irfInterface::AcceptanceCone *>
                        (&(roiCuts->extractionRegion())));
      s_emin = (roiCuts->getEnergyCuts()).first;
      s_emax = (roiCuts->getEnergyCuts()).second;
   }
}

double PointSource::Aeff::operator()(double cos_theta) const {
   double theta = acos(cos_theta)*180./M_PI;
   static double phi;

   ResponseFunctions * respFuncs = ResponseFunctions::instance();

   double myEffArea = 0;
   std::map<unsigned int, irfInterface::Irfs *>::iterator respIt
      = respFuncs->begin();

   for ( ; respIt != respFuncs->end(); respIt++) {
      irfInterface::IPsf *psf = respIt->second->psf();
      irfInterface::IAeff *aeff = respIt->second->aeff();

      double psf_val = psf->angularIntegral(m_energy, m_srcDir,
                                            theta, phi, s_cones);
      double aeff_val = aeff->value(m_energy, theta, phi);

      if (ResponseFunctions::useEdisp()) {
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

bool PointSource::overlapInterval(const std::pair<double, double> & interval1,
                                  std::pair<double, double> & interval2) {
   double start = std::max(interval1.first, interval2.first);
   double stop = std::min(interval1.second, interval2.second);
   if (start < stop) {
      interval2.first = start;
      interval2.second = stop;
      return true;
   }
   return false;
}

} // namespace Likelihood
