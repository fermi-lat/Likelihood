/**
 * @file MultipleBrokenPowerLaw.cxx
 * @brief User configurable multiply broken power-law.
 * @author J. Chiang
 *
 * $Header: $
 */

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "optimizers/dArg.h"
#include "Likelihood/MultipleBrokenPowerLaw.h"

namespace Likelihood {

MultipleBrokenPowerLaw::MultipleBrokenPowerLaw() {
   setMaxNumParams(100);
}

void MultipleBrokenPowerLaw::
addParams(double normalization, const std::vector<double> & photonIndexes,
          const std::vector<double> & breakEnergies) {
   if (photonIndexes.size() != breakEnergies.size() + 1) {
      throw std::runtime_error("Inconsistent number of photon indexes and "
                               "break energies.");
   }
   m_breakEnergies = breakEnergies;
   addParam("Normalization", normalization, true);
   for (size_t i(0); i < photonIndexes.size(); i++) {
      std::ostringstream parname;
      parname << "Index" << i;
      addParam(parname.str(), photonIndexes[i], true);
      m_indexNames.push_back(parname.str());
   }

   // Keep the break energies permanently fixed.
   for (size_t i(0); i < breakEnergies.size(); i++) {
      std::ostringstream parname;
      parname << "Break" << i;
      addParam(parname.str(), breakEnergies[i], false);
      setParamAlwaysFixed(parname.str());
      if (i > 0 && !(breakEnergies[i] > breakEnergies[i-1])) {
         throw std::runtime_error("Break energies must be monotonically "
                                  "increasing.");
      }
   }
}

double MultipleBrokenPowerLaw::value(optimizers::Arg & xarg) const {
   double x(dynamic_cast<optimizers::dArg &>(xarg).getValue());

   // Set the scale to be the first break energy.
   double x0(m_breakEnergies[0]);
   
   // Find the bracketing energies.
   size_t k(std::upper_bound(m_breakEnergies.begin(), m_breakEnergies.end(), x)
            - m_breakEnergies.begin());
   
   // The vector index for the photon index parameter is offset by 1 to
   // account for the Normalization parameter.
   double index(m_parameter[k + 1].getTrueValue());
   double value(norm(k)*std::pow(x/x0, index));
   return value;
}

double MultipleBrokenPowerLaw::
derivByParam(optimizers::Arg & xarg, const std::string & paramName) const {
   double x(dynamic_cast<optimizers::dArg &>(xarg).getValue());
   double x0(m_breakEnergies[0]);
   if (paramName == "Normalization") {
      double normalization(m_parameter[0].getValue());
      return value(xarg)/normalization;
   }
   std::vector<std::string>::const_iterator it
      = std::find(m_indexNames.begin(), m_indexNames.end(), paramName);
   if (it == m_indexNames.end()) {
      std::ostringstream what;
      what << "Likelihood::MultipleBrokenPowerLaw: "
           << "parameter name not found: " << paramName;
      throw std::runtime_error(what.str());
   }
   // Find the photon index vector index value.
   // Normalization is the first parameter, so add 1 to the location
   // in the m_indexNames vector to obtain the location in the
   // m_parameter vector.
   size_t m(it - m_indexNames.begin() + 1);
   if (m == 1) {
      if (x < m_breakEnergies[0]) {
         return value(xarg)*std::log(x/x0);
      } else {
         return 0;
      }
   }
   // Decrement the index by 1 to account for Normalization in
   // m_parameter[0] and by another 1 since that is the offset from
   // the derivative expression.
   double xmm1(m_breakEnergies[m-2]);
   if (x <= xmm1) {
      // Energies below the interval in question are not affected
      // by changes in the photon index.
      return 0;
   }
   // Find the bracketing energy break values.
   size_t k(std::upper_bound(m_breakEnergies.begin(), m_breakEnergies.end(), x)
            - m_breakEnergies.begin() + 1);
   if (m < k) {
      // Desired photon index is below the photon index corresponding
      // to the interval with the target energy, so the only effect on
      // the partial derivative is on the effective normalization at
      // m_breakEnergies[k-1].
      return value(xarg)*std::log(m_breakEnergies[m-1]/xmm1);
   }
   // We have a photon with target energy in the same range as the
   // photon index parameter with respect to which we are perfoming.
   return value(xarg)*std::log(x/xmm1)*m_parameter[k].getScale();
}

double MultipleBrokenPowerLaw::norm(size_t k) const {
   double n0(m_parameter[0].getTrueValue());
   if (k <= 1) {
      return n0;
   }
   double x0(m_breakEnergies[0]);
   double my_norm(n0);
   for (size_t m(1); m < k; m++) {
      double xm(m_breakEnergies[m]);
      // Increment m_breakEnergies index by 1 to account for
      // Normalization as first element of m_parameter.
      double pm(m_parameter[m+1].getTrueValue());
      double pm_p_1(m_parameter[m+2].getTrueValue());
      my_norm *= std::pow(xm/x0, pm - pm_p_1);
   }
   return my_norm;
}

} // namespace Likelihood
