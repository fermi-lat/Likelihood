/**
 * @file Drm.cxx
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/Drm.cxx,v 1.3 2011/06/14 22:41:50 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <stdexcept>

#include "irfInterface/Irfs.h"

#include "Likelihood/Drm.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"

namespace {
   double integrate(const std::vector<double> & x, 
                    const std::vector<double> & y) {
      if (x.size() != y.size()) {
         throw std::runtime_error("::integrate: x and y sizes do not match.");
      }
      double value(0);
      for (size_t i(0); i < x.size()-1; i++) {
         value += (y[i+1] + y[i])/2.*(x[i+1] - x[i]);
      }
      return value;
   }
} // anonymous namespace

namespace Likelihood {

Drm::Drm(double ra, double dec, const Observation & observation, 
         const std::vector<double> & ebounds, size_t npts) 
   : m_dir(ra, dec), m_observation(observation), m_npts(npts) {
   // Prepare the energy bounds array to be used for both true and
   // measured counts bins.
   m_ebounds.resize(ebounds.size());
   std::copy(ebounds.begin(), ebounds.end(), m_ebounds.begin());

   // Pad m_ebounds with extra energy bins at the beginning and at the
   // end for the convolution.  These intervals are used only by the
   // true energy bins.
   double de(std::log(m_ebounds.at(1)/m_ebounds.at(0)));
   m_ebounds.push_front(std::exp(std::log(m_ebounds.at(0)) - de));
   m_ebounds.push_back(std::exp(std::log(m_ebounds.back()) + de));

   compute_drm();
}

void Drm::convolve(const std::vector<double> & true_counts,
                   std::vector<double> & meas_counts) const {
   if (true_counts.size() != m_ebounds.size() - 3) {
      throw std::runtime_error("Drm::convolve: Size of true_counts "
                               "does not equal size of energy grid.");
   }
   std::deque<double> counts(true_counts.size());
   std::copy(true_counts.begin(), true_counts.end(), counts.begin());
   double value = 
      counts.at(0)*std::exp(std::log(m_ebounds.at(0)/m_ebounds.at(2))/
                            std::log(m_ebounds.at(1)/m_ebounds.at(2))*
                            std::log(counts.at(0)/counts.at(1)));
   counts.push_front(value);

   size_t nee(m_ebounds.size() - 1);
   size_t ncc(counts.size() - 1);
   value = counts.at(ncc-1)
      *std::exp(std::log(m_ebounds.at(nee)/m_ebounds.at(nee-2))/
                std::log(m_ebounds.at(nee-1)/m_ebounds.at(nee-2))*
                std::log(counts.at(ncc)/counts.at(ncc-1)));
   counts.push_back(value);

   meas_counts.resize(m_ebounds.size() - 3);
   for (size_t kp(0); kp < meas_counts.size(); kp++) {
      meas_counts[kp] = 0;
      for (size_t k(0); k < counts.size(); k++) {
         meas_counts[kp] += counts.at(k)*m_drm.at(k).at(kp);
      }
   }
}

void Drm::compute_drm() {
   m_drm.clear();
   for (size_t k(0); k < m_ebounds.size()-1; k++) {
      std::vector<double> row;
      for (size_t kp(0); kp < m_ebounds.size() - 3; kp++) {
         std::vector<double> emeas;
         get_emeas(kp, emeas);
         std::vector<double> disp(m_npts);
         get_disp(std::sqrt(m_ebounds.at(k)*m_ebounds.at(k+1)), emeas, disp);
         row.push_back(::integrate(emeas, disp));
      }
      m_drm.push_back(row);
   }
}

void Drm::get_emeas(size_t kp, std::vector<double> & emeas) const {
   emeas.clear();
   double estep((m_ebounds.at(kp+2) - m_ebounds.at(kp+1))/(m_npts-1));
   emeas.push_back(m_ebounds[kp+1]);
   for (size_t k(1); k < m_npts; k++) {
      emeas.push_back(emeas.back() + estep);
   }
}

void Drm::get_disp(double etrue, const std::vector<double> & emeas,
                   std::vector<double> & disp) const {
   disp.clear();
   const ResponseFunctions & resps(m_observation.respFuncs());
   const ExposureCube & expcube(m_observation.expCube());

   // Use phi-averged exposure
   double phi(-1);

   // Get the event types (usually just the list of conversion_type's)
   // and turn off phi-dependence temporarily.
   std::vector<bool> phideps;
   std::vector<int> evtTypes;
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator it;
   for (it = resps.begin(); it != resps.end(); ++it) {
      phideps.push_back(it->second->aeff()->usePhiDependence());
      it->second->aeff()->setPhiDependence(false);
      evtTypes.push_back(it->second->irfID());
   }

   size_t nmu(20);
   std::vector<double> mu_vals;
   double dmu(0.99/(nmu-1.));
   for (size_t i(0); i < nmu; i++) {
      mu_vals.push_back(1 - dmu*i);
   }

   for (size_t k(0); k < emeas.size(); k++) {
      std::vector<double> exposr(nmu, 0);
      std::vector<double> top(nmu, 0);
      size_t j(0);
      for (std::vector<double>::const_iterator mu(mu_vals.begin());
           mu != mu_vals.end(); ++mu, j++) {
         double theta(std::acos(*mu)*180./M_PI);
         double livetime(expcube.livetime(m_dir, *mu, phi));
         for (size_t i(0); i < evtTypes.size(); i++) {
            double aeff(resps.aeff(etrue, theta, phi, evtTypes[i]));
            double edisp(resps.edisp(emeas[k], etrue, theta, phi, evtTypes[i]));
            exposr[j] += aeff*livetime;
            top[j] += edisp*aeff*livetime;
         }
      }
      double numerator(::integrate(mu_vals, top));
      double denominator(::integrate(mu_vals, exposr));
      if (numerator == 0) {
         disp.push_back(0);
      } else {
         disp.push_back(numerator/denominator);
      }
   }

   // Restore status of phi-dependence.
   size_t i(0);
   for (it = resps.begin(); it != resps.end(); ++it, i++) {
      it->second->aeff()->setPhiDependence(phideps[i]);
   }
}

} // namespace Likelihood
