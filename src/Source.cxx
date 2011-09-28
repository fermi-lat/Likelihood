/** 
 * @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Source.cxx,v 1.20 2011/09/26 01:35:49 jchiang Exp $
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

Source::Source(const Source &rhs) {
// The deep copy of m_functions must be handled by the subclasses.
// Need to refactor this.
   m_name = rhs.m_name;
   m_srcType = rhs.m_srcType;
   m_useEdisp = rhs.m_useEdisp;
}

double Source::Npred() {
   optimizers::Function * specFunc = m_functions["Spectrum"];

// Evaluate the Npred integrand at the abscissa points contained in
// RoiCuts::energies().
   const RoiCuts & roiCuts = m_observation->roiCuts();
   const std::vector<double> & energies = roiCuts.energies();

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
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   // if (emin < energies.front() || emax > energies.back()) {
   //    std::ostringstream message;
   //    message << "Source::Npred(emin, emax):\n"
   //            << "emin = " << emin << "\n"
   //            << "emax = " << emax << "\n";
   //    throw std::out_of_range(message.str());
   // }

   double tol(1e-7);
   if (emin < energies.front()) {
      emin = energies.front();
   }
   if (emax > energies.back()) {
      emax = energies.back();
   }

   std::vector<double>::const_iterator first 
      = std::upper_bound(energies.begin(), energies.end(), emin);
   std::vector<double>::const_iterator last 
      = std::upper_bound(energies.begin(), energies.end(), emax);
   std::vector<double> my_energies(last - first);
   std::copy(first, last, my_energies.begin());
   size_t begin_offset = first - energies.begin();
   size_t end_offset = last - energies.begin();
   my_energies.insert(my_energies.begin(), emin);
   my_energies.push_back(emax);
   std::vector<double> exposure(last - first);
   std::copy(m_exposure.begin() + begin_offset,
             m_exposure.begin() + end_offset,
             exposure.begin());
   if (end_offset == energies.size()) {
      end_offset = energies.size() - 1;
   }
   double begin_exposure = (emin - energies.at(begin_offset - 1))
      /(energies.at(begin_offset) - energies.at(begin_offset - 1))
      *(m_exposure.at(begin_offset) - m_exposure.at(begin_offset - 1))
      + m_exposure.at(begin_offset - 1);
   double end_exposure = (emax - energies.at(end_offset - 1))
      /(energies.at(end_offset) - energies.at(end_offset - 1))
      *(m_exposure.at(end_offset) - m_exposure.at(end_offset - 1))
      + m_exposure.at(end_offset - 1);
   exposure.insert(exposure.begin(), begin_exposure);
   exposure.push_back(end_exposure);

   FuncMap::const_iterator my_func = m_functions.find("Spectrum");
   const optimizers::Function & specFunc = 
      const_cast<optimizers::Function &>(*my_func->second);

   std::vector<double> integrand(my_energies.size());
   for (unsigned int k = 0; k < my_energies.size(); k++) {
      optimizers::dArg eArg(my_energies.at(k));
      integrand.at(k) = specFunc(eArg)*exposure.at(k);
   }
   bool useLog;
   TrapQuad trapQuad(my_energies, integrand, useLog=true);
   return trapQuad.integral();
}

double Source::NpredDeriv(const std::string &paramName) {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
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

double Source::powerlaw_integral_est(double x1, double x2, 
                                     double y1, double y2, 
                                     double wt1, double wt2) {
   double gam(std::log(y2/y1)/std::log(x2/x1));
   double y0(y2/std::pow(x2, gam));

   double xbar;
   if (gam == -1) {
      xbar = (x2 - x1)/std::log(x2/x1);
   } else if (gam == -2) {
      xbar = std::log(x2/x1)/(1./x1 - 1./x2);
   } else {
      xbar = (gam + 1.)/(gam + 2.)*(std::pow(x2, gam+2.) - std::pow(x1, gam+2.))
         /(std::pow(x2, gam+1.) - std::pow(x1, gam+1.));
   }

   double wtbar((xbar - x1)/(x2 - x1)*(wt2 - wt1) + wt1);

   if (gam == -1) {
      return wtbar*y0*std::log(x2/x1);
   }
   return wtbar*y0/(gam+1.)*(std::pow(x2, gam+1.) - std::pow(x1, gam+1.));
}

} // namespace Likelihood
