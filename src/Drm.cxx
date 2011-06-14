/**
 * @file Drm.cxx
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Drm.cxx,v 1.1 2011/06/14 06:31:53 jchiang Exp $
 */

#include "Likelihood/Drm.h"

namespace {
   double integrate(const std::vector<double> & x, 
                    const std::vector<double> & y) {
      
   }
} // anonymous namespace

//namespace Likelihood {

Drm::Drm(double ra, double dec, Observation & observation, 
         const std::vector<double> & ebounds, size_t npts) 
   : m_dir(ra, dec), m_observation(observation), m_npts(npts) {
   // Prepare the energy bounds array to be used for both true and
   // measured counts bins.
   m_ebounds.resize(ebounds.size());
   std::copy(ebounds.begin(), ebounds.end(), m_ebounds.begin());

   // Insert an interval at the beginning of m_ebounds for the
   // convolution.  This interval is used only by the true energy
   // bins.
   double de(std::log(m_ebounds.at(1)/m_ebounds.at(0)));
   m_ebounds.push_front(std::exp(std::log(m_ebounds[0]) - de));

   compute_drm();
}

void Drm::convolve(const std::vector<double> & true_counts,
                   std::vector<double> & meas_counts) const {
   if (true_counts.size() != m_ebounds.size() - 2) {
      throw std::runtime_error("Drm::convolve: Size of true_counts "
                               "does not equal size of energy grid.");
   }
   std::deque<double> counts(true_counts.size());
   std::copy(true_counts.begin(), true_counts.end(), counts.begin());
   double value = 
      true_counts[0]*std::exp(std::log(m_ebounds[0]/m_ebounds[2])/
                              std::log(m_ebounds[1]/m_ebounds[2])*
                              std::log(true_counts[0]/true_counts[1]));
   counts.push_front(value);

   meas_counts.resize(m_ebounds.size() - 2);
   for (size_t kp(0); kp < meas_counts.size(); kp++) {
      meas_counts[kp] = 0;
      for (size_t k(0); k < counts.size(); k++) {
         meas_counts[kp] += counts[k]*m_drm[kp][k];
      }
   }
}

void Drm::compute_drm() {
   m_drm.clear();
   for (size_t k(0); k < m_ebounds.size()-1; k++) {
      std::vector<double> row(m_ebounds.size() - 2);
      for (size_t kp(0); kp < row.size(); kp++) {
         std::vector<double> emeas(npts);
         get_emeas(kp, emeas);
         std::vector<double> disp(npts);
         get_disp(std::sqrt(m_ebounds[k]*m_ebounds[k+1]), emeas, disp);
         row[k] = integrate(emeas, disp);
      }
      m_drm.push_back(row);
   }
}

void Drm::get_emeas(size_t kp, std::vector<double> & emeas) const {
   double estep((m_ebounds[kp+2] - m_ebounds[kp+1])/(m_npts-1));
   emeas[0] = m_ebounds[kp+1];
   for (size_t k(1); k < m_npts; k++) {
      emeas[k] = emeas[k-1] + estep;
   }
}

void Drm::get_disp(double etrue, const std::vector<double> & emeas,
                   std::vector<double> & disp) const {
   disp.clear();
   const ResponseFunctions & resps(m_observation.respFuncs());
   const ExposureCube & expcube(m_observation.expCube());

   // Get the event types (usually just the list of conversion_type's)
   std::vector<int> evtTypes;
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator it;
   for (it = resps.begin(); it != resps.end(); ++it) {
      evtTypes.push_back(it->second->irfID());
   }

   size_t nmu(20);
   std::vector<double> mu_vals;
   double dmu(1./(nmu-1.));
   for (size_t i(0); i < nmu; i++) {
      mu_vals.push_back(1 - dmu*i);
   }

   for (size_t k(0); k < emeas.size(); k++) {
      std::vector<double> exposr(nmu, 0);
      std::vector<double> top(nmu, 0);
      size_t j(0);
      for (std::vector<double>::const_iterator mu(mu_vals.begin());
           mu != mu_vals.end(); ++mu, j++) {
         for (size_t i(0); i < evtTypes.size(); i++) {
            double theta(std::acos(*mu));
            double aeff(resps.aeff(etrue, theta, phi, evtTypes[i]));
            double livetime(expcube.livetime(m_dir, *mu));
            double edisp(resps.edisp(emeas, etrue, theta, phi, evtTypes[i]));
            exposr[j] += aeff*livetime;
            top[j] += edisp*exposr[j];
         }
      }
      disp.push_back(integrate(mu_vals, top)/integrate(mu_vals, exposr));
   }
}

} // namespace Likelihood
