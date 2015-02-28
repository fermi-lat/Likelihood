/**
 * @file PiecewisePowerLaw.cxx
 * @brief User configurable multiply broken power-law.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/PiecewisePowerLaw.cxx,v 1.3 2015/02/26 17:32:31 jchiang Exp $
 */

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "optimizers/dArg.h"
#include "Likelihood/PiecewisePowerLaw.h"

namespace Likelihood {

PiecewisePowerLaw::PiecewisePowerLaw() {
   setMaxNumParams(100);
   m_funcType = Addend;
   m_argType = "dArg";
   m_genericName = "PiecewisePowerLaw";
   m_normParName = "";
}

void PiecewisePowerLaw::
addParams(double indexL, double indexH, 
          const std::vector<double> & dNdEs,
          const std::vector<double> & energies) {
   if (dNdEs.size() != energies.size()) {
      throw std::runtime_error("Inconsistent number of dNdE values and "
                               "energies.");
   }
   m_energies = energies;
   addParam("IndexL", indexL, true);
   addParam("IndexH", indexH, true);
   for (size_t k(0); k < dNdEs.size(); k++) {
      std::ostringstream parname;
      parname << "dNdE" << k;
      addParam(parname.str(), dNdEs[k], true);
      m_dNdENames.push_back(parname.str());
   }

   // Keep the energies permanently fixed.
   for (size_t k(0); k < energies.size(); k++) {
      std::ostringstream parname;
      parname << "Energy" << k;
      addParam(parname.str(), energies[k], false);
      setParamAlwaysFixed(parname.str());
      if (k > 0 && !(energies[k] > energies[k-1])) {
         throw std::runtime_error("Energy parameters must be monotonically "
                                  "increasing.");
      }
   }
}

double PiecewisePowerLaw::value(optimizers::Arg & xarg) const {
   double x(dynamic_cast<optimizers::dArg &>(xarg).getValue());

   // Handle low and high end cases first.
   if (x < m_energies.front()) { // IndexL
      return norm(0)*std::pow(x/m_energies.front(),
                              m_parameter.front().getTrueValue());
   } else if (x > m_energies.back()) { // IndexH
      return norm(m_dNdENames.size()-1)*std::pow(x/m_energies.back(),
                                                 m_parameter[1].getTrueValue());
   }
   
   // Find the bracketing energies.
   size_t k(std::upper_bound(m_energies.begin(), m_energies.end(), x)
            - m_energies.begin() - 1);
   
   double value(norm(k)*std::pow(x/m_energies[k], plIndex(k)));
   return value;
}

double PiecewisePowerLaw::
derivByParam(optimizers::Arg & xarg, const std::string & paramName) const {
   if (paramName.substr(0, 6) == "Energy") {
      throw std::runtime_error("MultipleBPL: Parameter " + paramName 
                               + " must be fixed in the xml model definition.");
   }
   double x(dynamic_cast<optimizers::dArg &>(xarg).getValue());

   if (paramName == "IndexL") {
      if (x < m_energies.front()) {
         return value(xarg)*std::log(x/m_energies.front());
      } else {
         return 0;
      }
   }
   if (paramName == "IndexH") {
      if (x > m_energies.back()) {
         return value(xarg)*std::log(x/m_energies.back());
      } else {
         return 0;
      }
   }
   std::vector<std::string>::const_iterator it
      = std::find(m_dNdENames.begin(), m_dNdENames.end(), paramName);
   if (it == m_dNdENames.end()) {
      std::ostringstream what;
      what << "Likelihood::PiecewisePowerLaw: "
           << "parameter name not found: " << paramName;
      throw std::runtime_error(what.str());
   }
   // Extract the index of dNdE parameter from the name.
   if (paramName.substr(0, 4) != "dNdE") {
      throw std::runtime_error("PiecewisePowerLaw: Parameter not recognized: " 
                               + paramName);
   }
   size_t k(std::atoi(paramName.substr(4).c_str()));

   // End points where target energy is outside m_energies.
   if ((k == 0 && x < m_energies.front()) ||
       (k == m_dNdENames.size()-1 && x > m_energies.back())) {
      return value(xarg)/norm(k);
   }
   // End points where target energy is inside m_energies.
   if (k == 0 && x <= m_energies[1]) {
      return value(xarg)/norm(k)*(1. - std::log(x/m_energies[k])
                                  /std::log(m_energies[k+1]/m_energies[k]));
   }
   if (k == m_dNdENames.size()-1 && x > m_energies[k-1]) {
      return value(xarg)/norm(k)*(std::log(x/m_energies[k-1])
                                  /std::log(m_energies[k]/m_energies[k-1]));
   }

   // All other non-zero cases
   if (m_energies[k] <= x && (k < m_energies.size()-1 && x < m_energies[k+1])) {
      return value(xarg)/norm(k)*(1. - std::log(x/m_energies[k])
                                  /std::log(m_energies[k+1]/m_energies[k]));
   } else if ((k > 0 && x > m_energies[k-1]) && x < m_energies[k]) {
      return value(xarg)/norm(k)*(std::log(x/m_energies[k-1])
                                  /std::log(m_energies[k]/m_energies[k-1]));
   }
   return 0;
}

double PiecewisePowerLaw::norm(size_t k) const {
   if (k < 0 || k >= m_dNdENames.size()) {
      throw std::out_of_range("PiecewisePowerLaw::norm");
   }
   return m_parameter[2 + k].getTrueValue();
}

double PiecewisePowerLaw::plIndex(size_t k) const {
   if (k < 0 || (k+1 >= m_energies.size())) {
      throw std::out_of_range("PiecewisePowerLaw::plIndex: "
                              "out-of-range energy index.");
   }
   double dNdE0(m_parameter[k+2].getTrueValue());
   double dNdE1(m_parameter[k+3].getTrueValue());
   double gamma(std::log(dNdE1/dNdE0)/std::log(m_energies[k+1]/m_energies[k]));
   return gamma;
}

} // namespace Likelihood
