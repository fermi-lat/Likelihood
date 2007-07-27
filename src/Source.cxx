/** 
 * @file Source.cxx
 * @brief Source class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Source.cxx,v 1.8 2006/12/24 19:52:35 jchiang Exp $
 */

#include <algorithm>
#include <stdexcept>

#include "optimizers/dArg.h"

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
   TrapQuad trapQuad(energies, NpredIntegrand);
   double value(trapQuad.integral());
   return value;
}

double Source::Npred(double emin, double emax) const {
   const std::vector<double> & energies = m_observation->roiCuts().energies();
   if (fabs((emin - energies.front())/emin) < 1e-2) {
      emin = energies.front();
   }
   if (fabs((emax - energies.back())/emax) < 1e-2) {
      emax = energies.back();
   }
   if (emin < energies.front() || emax > energies.back()) {
      throw std::out_of_range("Source::Npred(emin, emax)");
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
   TrapQuad trapQuad(my_energies, integrand);
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
      TrapQuad trapQuad(energies, myIntegrand);
      return trapQuad.integral();
   }
}

} // namespace Likelihood
