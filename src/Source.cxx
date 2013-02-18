/** 
 * @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Source.cxx,v 1.24 2013/01/09 00:44:41 jchiang Exp $
 */

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "optimizers/dArg.h"
#include "optimizers/FunctionFactory.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/Observation.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/Source.h"
#include "Likelihood/TrapQuad.h"

namespace Likelihood {

Source::Source(const Observation * observation) 
   : m_name(""), m_srcType(""), m_useEdisp(false), m_spectrum(0), 
     m_observation(observation) {}

Source::Source(const Source & rhs)
   : m_name(rhs.m_name),
     m_srcType(rhs.m_srcType),
     m_useEdisp(rhs.m_useEdisp),
     m_spectrum(rhs.m_spectrum->clone()),
     m_observation(rhs.m_observation),
     m_exposure(rhs.m_exposure),
     m_energies(rhs.m_energies) {
// The deep copy of m_functions must be handled by the subclasses.
}

double PointSource::fluxDensity(const Event & evt, const Response& resp) {
// Scale the energy spectrum by the psf value and the effective area
// and convolve with the energy dispersion, if appropriate, all of
// which are functions of time and spacecraft attitude and orbital
// position.

   const double& energy         = evt.getEnergy();

   Response my_resp;
   const Response* the_resp = &resp;
   if(resp.empty()) {
      computeResponse(my_resp, evt);
      the_resp = &my_resp;
   }

   if(m_edisp() == ED_NONE) {
      DEBUG_ASSERT(the_resp->size() == 1);
      optimizers::dArg energy_arg(energy);
      double spectrum = (*m_spectrum)(energy_arg);
      return spectrum * (*the_resp)[0];
   } else {
      const std::vector<double> & true_energies = evt.trueEnergies();
      unsigned int nenergy(true_energies.size());
      if(m_edisp() == ED_GAUSSIAN) {
	 DEBUG_ASSERT(the_resp->size() == 3);
	 fillGaussianResponse(my_resp, *the_resp, true_energies);
	 the_resp = &my_resp;
      } else {
 	 DEBUG_ASSERT(the_resp->size() == nenergy);
      }
 
      std::vector<double> my_integrand(nenergy);
      for(unsigned ienergy=0;ienergy<nenergy;ienergy++)
	{
	   optimizers::dArg energy_arg(true_energies[ienergy]);
	   double spectrum = (*m_spectrum)(energy_arg);
 	   my_integrand[ienergy]+= spectrum * (*the_resp)[ienergy];
	}
      TrapQuad trapQuad(true_energies, my_integrand, evt.trueEnergiesUseLog());
      return trapQuad.integral();
   }
}
      
double PointSource::fluxDensityDeriv(const Event & evt, 
				     const std::string & paramName,
				     const Response & resp) const
{
// For now, just implement for spectral Parameters and neglect
// the spatial ones, "longitude" and "latitude"

   // Special shortcut for "prefactor"
   double prefactor;
   if (paramName == "Prefactor" && 
       (prefactor = m_spectrum->getParamValue("Prefactor")) !=0) {
      return fluxDensity(evt, resp)/prefactor;
   }

   const double& energy         = evt.getEnergy();

   Response my_resp;
   const Response* the_resp = &resp;
   if(resp.empty()) {
      computeResponse(my_resp, evt);
      the_resp = &my_resp;
   }

   if(m_edisp() == ED_NONE) {
      DEBUG_ASSERT(the_resp->size() == 1);
      optimizers::dArg energy_arg(energy);
      double spectrum_deriv = m_spectrum->derivByParam(energy_arg, paramName);
      return spectrum_deriv * (*the_resp)[0];
   } else {
      const std::vector<double> & true_energies = evt.trueEnergies();
      unsigned int nenergy(true_energies.size());
      if(m_edisp() == ED_GAUSSIAN) {
	 DEBUG_ASSERT(the_resp->size() == 3);
	 fillGaussianResponse(my_resp, *the_resp, true_energies);
	 the_resp = &my_resp;
      } else {
 	 DEBUG_ASSERT(the_resp->size() == nenergy);
      }
 
      std::vector<double> my_integrand(nenergy);
      for(unsigned ienergy=0;ienergy<nenergy;ienergy++)
	{
	   optimizers::dArg energy_arg(true_energies[ienergy]);
	   double spectrum_deriv =
	     m_spectrum->derivByParam(energy_arg, paramName);
 	   my_integrand[ienergy]+= spectrum_derive * (*the_resp)[ienergy];
	}
      TrapQuad trapQuad(true_energies, my_integrand, evt.trueEnergiesUseLog());
      return trapQuad.integral();
   }   
}

void Source::computeExposure(bool verbose) {
   m_energies = m_observation->roiCuts().energies();
   computeExposure(m_energies, verbose);
}

double Source::Npred() {
   optimizers::Function * specFunc = m_functions["Spectrum"];
   if (specFunc->xvalues().size() == 0) {
      const std::vector<double> & energies(m_energies);

      std::vector<double> NpredIntegrand(energies.size());
      for (unsigned int k = 0; k < energies.size(); k++) {
         optimizers::dArg eArg(energies[k]);
         NpredIntegrand[k] = (*specFunc)(eArg)*m_exposure[k];
      }
      bool useLog;
      TrapQuad trapQuad(energies, NpredIntegrand, useLog=true);
      double value(trapQuad.integral());
      return value;
   }
   const std::vector<double> & energies(specFunc->xvalues());
   std::vector<double> exposure;
   getExposureValues(energies, exposure);
   std::vector<double> NpredIntegrand(energies.size());
   for (unsigned int k = 0; k < energies.size(); k++) {
      optimizers::dArg eArg(energies[k]);
      NpredIntegrand[k] = (*specFunc)(eArg)*m_exposure[k];
   }
   bool useLog;
   TrapQuad trapQuad(energies, NpredIntegrand, useLog=true);
   double value(trapQuad.integral());
   return value;
}

void Source::getExposureValues(const std::vector<double> & energies,
                               std::vector<double> & exposures) const {
   exposures.resize(energies.size(), 0);
   double estep(std::log(m_energies[1]/m_energies[0]));
   for (size_t k(0); k < energies.size(); k++) {
      if (energies[k] < m_energies[0] || energies[k] > m_energies.back()) {
         exposures[k] = 0;
      } else {
         size_t indx = static_cast<size_t>(std::log(energies[k]/m_energies[0])
                                           /estep);
         exposures[k] = (m_exposure[indx]
                         *std::exp(std::log(energies[k]/m_energies[indx])/estep
                                   *std::log(m_exposure[indx+1]/m_exposure[indx])));
      }
   }
}

double Source::Npred(double emin, double emax) const {
   std::vector<double> energies;
   std::vector<double> exposure;

//   getExposureSubArrays(emin, emax, energies, exposure);
   getExposureArrays(emin, emax, energies, exposure);
//    std::cout << m_name << ", Npred(emin, emax): "
//              << energies.size() << "  "
//              << energies.front() << "  "
//              << energies.back() << std::endl;

   FuncMap::const_iterator my_func = m_functions.find("Spectrum");
   const optimizers::Function & specFunc = 
      const_cast<optimizers::Function &>(*my_func->second);

   std::vector<double> integrand(energies.size());
   for (unsigned int k = 0; k < energies.size(); k++) {
      optimizers::dArg eArg(energies.at(k));
      integrand.at(k) = specFunc(eArg)*exposure.at(k);
   }
   bool useLog;
   TrapQuad trapQuad(energies, integrand, useLog=true);
   return trapQuad.integral();
}

double Source::NpredDeriv(const std::string &paramName) {
   const std::vector<double> & energies(m_energies);
   optimizers::Function *specFunc = m_functions["Spectrum"];

   double prefactor;
   if (paramName == std::string("Prefactor") && 
       (prefactor = specFunc->getParamValue("Prefactor")) != 0) {
      return Npred()/prefactor;
   } else {  // loop over energies and fill integrand vector
      std::vector<double> myIntegrand(energies.size());
      for (unsigned int k = 0; k < energies.size(); k++) {
         optimizers::dArg eArg(energies[k]);
         myIntegrand[k] = specFunc->derivByParam(eArg, paramName)
            *m_exposure[k];
      }
      bool useLog;
      TrapQuad trapQuad(energies, myIntegrand, useLog=true);
      return trapQuad.integral();
   }
}

double Source::
NpredDeriv(const std::string & paramName, double emin, double emax) const {
   std::vector<double> energies;
   std::vector<double> exposures;
//   getExposureSubArrays(emin, emax, energies, exposures);
   getExposureArrays(emin, emax, energies, exposures);

   const optimizers::Function * specFunc 
      = m_functions.find("Spectrum")->second;

   double prefactor;
   if (paramName == std::string("Prefactor") && 
       (prefactor = specFunc->getParamValue("Prefactor")) != 0) {
      return Npred(emin, emax)/prefactor;
   } else {  // loop over energies and fill integrand vector
      std::vector<double> myIntegrand(energies.size());
      for (unsigned int k = 0; k < energies.size(); k++) {
         optimizers::dArg eArg(energies[k]);
         myIntegrand[k] = specFunc->derivByParam(eArg, paramName)
            *exposures[k];
      }
      bool useLog;
      TrapQuad trapQuad(energies, myIntegrand, useLog=true);
      return trapQuad.integral();
   }
}

void Source::setSpectrum(const std::string & functionName) {
   optimizers::FunctionFactory funcFactory;
   AppHelpers::addFunctionPrototypes(&funcFactory);
   try {
      m_spectrum = funcFactory.create(functionName);
      m_functions["Spectrum"] = m_spectrum;
   } catch(optimizers::Exception & eObj) {
      std::ostringstream message;
      message << eObj.what() << "\n"
              << "Available function names:\n";
      std::vector<std::string> names;
      funcFactory.getFunctionNames(names);
      for (size_t i(0); i < names.size(); i++) {
         message << "  " << names.at(i) << std::endl;
      }
      throw Exception(message.str());
   }
}

bool Source::fixedSpectrum() const {
   return m_spectrum->getNumFreeParams() == 0;
}

double Source::pixelCounts(double emin, double emax,
                           double wtMin, double wtMax) const {
   optimizers::Function & spectrum = *m_spectrum;
   optimizers::dArg eminArg(emin);
   optimizers::dArg emaxArg(emax);

   double f1(spectrum(eminArg));
   double f2(spectrum(emaxArg));

   double y1(f1*wtMin);
   double y2(f2*wtMax);
   if (::getenv("USE_LINEAR_QUADRATURE")) {
      return (y1 + y2)*(emax - emin)/2.;
   }
// Quadrature in log-log space, suggested by J. Ballet.   
   return (y1*emin + y2*emax)/2.*std::log(emax/emin);
}

double Source::pixelCountsDeriv(double emin, double emax,
                                double wtMin, double wtMax,
                                const std::string & paramName) const {
   optimizers::Function & spectrum = *m_spectrum;
   optimizers::dArg eminArg(emin);
   optimizers::dArg emaxArg(emax);
   // double y1(spectrum(eminArg)*wtMin);
   // double y2(spectrum(emaxArg)*wtMax);

   double f1(spectrum.derivByParam(eminArg, paramName));
   double f2(spectrum.derivByParam(emaxArg, paramName));

   double dy1dp(f1*wtMin);
   double dy2dp(f2*wtMax);
   if (::getenv("USE_LINEAR_QUADRATURE")) {
      return (dy1dp + dy2dp)*(emax - emin)/2.;
   }
// Quadrature in log-log space, suggested by J. Ballet.   
   return (dy1dp*emin + dy2dp*emax)/2.*std::log(emax/emin);
}

void Source::getExposureArrays(double emin, double emax, 
                               std::vector<double> & energies,
                               std::vector<double> & exposures,
                               size_t nee) const {
   const std::vector<double> & roi_energies(m_energies);
   if (emin < roi_energies.front()) {
      emin = roi_energies.front();
   }
   if (emax > roi_energies.back()) {
      emax = roi_energies.back();
   }
   if (emin == roi_energies.front() && emax == roi_energies.back()) {
      energies = roi_energies;
      exposures = m_exposure;
      return;
   }
   if (nee == 0) {
      nee = roi_energies.size();
   }
   double roi_estep(std::log(roi_energies.back()/roi_energies.front())
                    /(roi_energies.size() - 1));
   double logestep(std::log(emax/emin)/(nee - 1));
   energies.clear();
   exposures.clear();
   for (size_t k(0); k < nee; k++) {
      double logE(k*logestep + std::log(emin));
      energies.push_back(std::exp(logE));
      size_t kk = static_cast<size_t>((logE - std::log(roi_energies.front()))
                                      /roi_estep);
      double exposure(0);
      if (kk == roi_energies.size()-1) { 
         exposure = m_exposure.back();
      } else {
         // log-linear interpolation
         exposure = m_exposure.at(kk) +
            (logE - std::log(roi_energies.at(kk)))/roi_estep
            *(m_exposure.at(kk+1) - m_exposure.at(kk));
      }
      exposures.push_back(exposure);
   }
}

void Source::getExposureSubArrays(double emin, double emax,
                                  std::vector<double> & energies,
                                  std::vector<double> & exposures) const {
   const std::vector<double> & roi_energies(m_energies);

   if (emin < roi_energies.front()) {
      emin = roi_energies.front();
   }
   if (emax > roi_energies.back()) {
      emax = roi_energies.back();
   }

   std::vector<double>::const_iterator first 
      = std::upper_bound(roi_energies.begin(), roi_energies.end(), emin);
   std::vector<double>::const_iterator last 
      = std::upper_bound(roi_energies.begin(), roi_energies.end(), emax);
   energies.resize(last - first);
   std::copy(first, last, energies.begin());
   size_t begin_offset = first - roi_energies.begin();
   size_t end_offset = last - roi_energies.begin();
   energies.insert(energies.begin(), emin);
   energies.push_back(emax);

   exposures.resize(last - first);
   std::copy(m_exposure.begin() + begin_offset,
             m_exposure.begin() + end_offset,
             exposures.begin());
   if (end_offset == roi_energies.size()) {
      end_offset = roi_energies.size() - 1;
   }
   double begin_exposure = (emin - roi_energies.at(begin_offset - 1))
      /(roi_energies.at(begin_offset) - roi_energies.at(begin_offset - 1))
      *(m_exposure.at(begin_offset) - m_exposure.at(begin_offset - 1))
      + m_exposure.at(begin_offset - 1);
   double end_exposure = (emax - roi_energies.at(end_offset - 1))
      /(roi_energies.at(end_offset) - roi_energies.at(end_offset - 1))
      *(m_exposure.at(end_offset) - m_exposure.at(end_offset - 1))
      + m_exposure.at(end_offset - 1);
   exposures.insert(exposures.begin(), begin_exposure);
   exposures.push_back(end_exposure);
}

void Source::computeGaussianParams(Response& gaussian_params,
				   const Response& resp,
				   const std::vector<double>& trueEnergies,
				   bool trueEnergiesUseLog) const {
   std::vector<double> y = resp;
   unsigned ne = trueEnergies.size();   
   double norm = TrapQuad(trueEnergies, y, trueEnergiesUseLog).integral();
   for(unsigned ie=0;ie<ne;ie++)y[ie] *= trueEnergies[ie];
   double int1 = TrapQuad(trueEnergies, y, trueEnergiesUseLog).integral();
   for(unsigned ie=0;ie<ne;ie++)y[ie] *= trueEnergies[ie];
   double int2 = TrapQuad(trueEnergies, y, trueEnergiesUseLog).integral();
   params.resize(3);
   params[0] = norm;
   params[1] = int1/norm;
   params[2] = std::sqrt(int2/norm - params[1]*params[1]);
} 

void Source::fillGaussianResp(Response& resp, const Response& gaussian_params,
			      const std::vector<double>& trueEnergies,
			      bool trueEnergiesUseLog) const {
   double norm = gaussian_params[0];
   double mean = gaussian_params[1];
   double sdev = gaussian_params[2];
   unsigned ne = trueEnergies.size();
   resp.resize(ne);
   norm /= sdev*std::sqrt(2.0*M_PI);
   for(unsigned ie=0;ie<ne;ie++) {
      double x = (trueEnergies[ie] - mean)/sdev;
      resp[ie] = norm*std::exp(-0.5*x*x);
   }
}

double value(trapQuad.integral());
  return value;
 }

} // namespace Likelihood
