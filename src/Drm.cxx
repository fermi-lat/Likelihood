/**
 * @file Drm.cxx
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Drm.cxx,v 1.11 2014/09/22 21:12:47 jchiang Exp $
 */

#include <cmath>

#include <algorithm>
#include <stdexcept>

#include "irfInterface/Irfs.h"

#include "Likelihood/Drm.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/Source.h"


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
   double de(std::log(m_ebounds[1]/m_ebounds[0]));
   m_ebounds.push_front(std::exp(std::log(m_ebounds[0]) - de));
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
   double value(0);
   double min_counts_value(1e-10);
   if (counts[0] > min_counts_value && counts[1] > min_counts_value) {
      value = counts[1]*std::exp(std::log(m_ebounds[0]/m_ebounds[2])/
                                 std::log(m_ebounds[1]/m_ebounds[2])*
                                 std::log(counts[0]/counts[1]));
   } else {
      value = ((m_ebounds[0] - m_ebounds[2])/
               (m_ebounds[1] - m_ebounds[2])*
               (counts[0] - counts[1]) + counts[1]);
   }
   counts.push_front(value);

   size_t nee(m_ebounds.size() - 1);
   size_t ncc(counts.size() - 1);
   if (counts[ncc] > 0 && counts[ncc-1] > 0) {
      value = counts[ncc-1]
         *std::exp(std::log(m_ebounds[nee]/m_ebounds[nee-2])/
                   std::log(m_ebounds[nee-1]/m_ebounds[nee-2])*
                   std::log(counts[ncc]/counts[ncc-1]));
   } else {
      value = ((m_ebounds[nee] - m_ebounds[nee-2])/
               (m_ebounds[nee-1] - m_ebounds[nee-2])*
               (counts[ncc] - counts[ncc-1]) + counts[ncc-1]);
   }
   counts.push_back(value);

   meas_counts.resize(m_ebounds.size() - 3);
   for (size_t kp(0); kp < meas_counts.size(); kp++) {
      meas_counts[kp] = 0;
      for (size_t k(0); k < counts.size(); k++) {
         meas_counts[kp] += counts[k]*m_drm[k][kp];
      }
   }
}

void Drm::compute_drm() {
   m_drm.clear();
   for (size_t k(0); k < m_ebounds.size()-1; k++) {
      std::vector<double> row;
      for (size_t kp(0); kp < m_ebounds.size() - 3; kp++) {
         double emeas_min(m_ebounds[kp+1]);
         double emeas_max(m_ebounds[kp+2]);
         row.push_back(matrix_element(std::sqrt(m_ebounds[k]*m_ebounds[k+1]),
                                      emeas_min, emeas_max));
      }
      m_drm.push_back(row);
   }
}

double Drm::
matrix_element(double etrue, double emeas_min, double emeas_max) const {
   const ResponseFunctions & resps(m_observation.respFuncs());
   const ExposureCube & expcube(m_observation.expCube());
   double met((expcube.tstart() + expcube.tstop())/2.);

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

   std::vector<double> exposr(nmu, 0);
   std::vector<double> top(nmu, 0);
   size_t j(0);
   for (std::vector<double>::const_iterator mu(mu_vals.begin());
        mu != mu_vals.end(); ++mu, j++) {
      double theta(std::acos(*mu)*180./M_PI);
      double livetime(expcube.livetime(m_dir, *mu, phi));
      for (size_t i(0); i < evtTypes.size(); i++) {
         double aeff(resps.aeff(etrue, theta, phi, evtTypes[i], met));
         double edisp(resps.edisp(evtTypes[i]).integral(emeas_min, emeas_max, 
                                                        etrue, theta, phi,
                                                        met));
         exposr[j] += aeff*livetime;
         top[j] += edisp*aeff*livetime;
      }
   }
   double numerator(::integrate(mu_vals, top));
   double denominator(::integrate(mu_vals, exposr));
   double my_disp(0);
   if (numerator != 0) {
      my_disp = numerator/denominator;
   }

   // Restore status of phi-dependence.
   size_t i(0);
   for (it = resps.begin(); it != resps.end(); ++it, i++) {
      it->second->aeff()->setPhiDependence(phideps[i]);
   }

   return my_disp;
}


Drm_Cache::Drm_Cache(const Drm& drm,
		     const SourceMap & sourceMap,
		     const Source& src,
		     const std::vector<double>& energies,
		     bool use_edisp)
  :m_true_counts(energies.size()-1),
   m_meas_counts(energies.size()-1),
   m_xi(energies.size()-1),
   m_kref(energies.size()-1),
   m_use_edisp(use_edisp){
  update(drm,sourceMap,src,energies,use_edisp);
}

Drm_Cache::Drm_Cache(const Drm_Cache& other) 
  : m_true_counts(other.m_true_counts),
    m_meas_counts(other.m_meas_counts),
    m_xi(other.m_xi),
    m_kref(other.m_kref),
    m_use_edisp(other.m_use_edisp){
}

void Drm_Cache::update(const Drm& drm,
		       const SourceMap & sourceMap,
		       const Source& src,
		       const std::vector<double>& energies,
		       bool use_edisp) {
  const std::vector<double> & npreds = sourceMap.npreds();
  const std::vector<std::pair<double,double> >& npred_weights = sourceMap.npred_weights();
  size_t k(0);
  for (k = 0; k < energies.size()-1; k++) {
    m_true_counts[k] = src.pixelCounts(energies.at(k), energies.at(k+1),
				       npreds.at(k)*npred_weights[k].first, npreds.at(k+1)*npred_weights[k].second);
  }
  if ( use_edisp ) {
    drm.convolve(m_true_counts, m_meas_counts);
  } else {
    m_meas_counts.resize(m_true_counts.size());
    std::copy(m_true_counts.begin(),m_true_counts.end(),m_meas_counts.begin());
  }
  int kref(-1);
  for (k = 0; k < energies.size()-1; k++) {
    if ( m_true_counts[k] > 0 ) {
      // Still have counts in this true energy bin, so it can 
      // be used as a reference
      m_xi[k] = m_meas_counts[k] / m_true_counts[k];
      kref = k;
      m_kref[k] = -1;
    } else {
      // Don't have counts in this true energy bin.  
      m_xi[k] = m_meas_counts[k] / m_true_counts[kref];
      m_kref[k] = kref;
    }
  }
}

} // namespace Likelihood
