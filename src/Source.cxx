/** 
 * @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Source.cxx,v 1.22 2012/02/08 00:23:55 jchiang Exp $
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

double Source::Npred() {
   optimizers::Function * specFunc = m_functions["Spectrum"];
   const std::vector<double> & energies(m_energies);
//    std::cout << m_name << ", Npred(): " 
//              << energies.size() << "  "
//              << energies.front() << "  "
//              << energies.back() << std::endl;

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

} // namespace Likelihood
