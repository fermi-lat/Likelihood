/**
 * @file Drm.cxx
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Drm.cxx,v 1.23 2017/10/12 21:54:14 echarles Exp $
 */

#include <cmath>

#include <algorithm>
#include <stdexcept>

#include "irfInterface/Irfs.h"

#include "Likelihood/Drm.h"
#include "Likelihood/FitUtils.h"
#include "Likelihood/ExposureCube.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/Source.h"
#include "Likelihood/WeightMap.h"


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
  : m_dir(ra, dec), m_observation(&observation), m_npts(npts){
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

   size_t nee(m_ebounds.size() - 1);

   m_log_ratio_lo = std::log(m_ebounds[0]/m_ebounds[2]) / std::log(m_ebounds[1]/m_ebounds[2]);
   m_diff_ratio_lo = (m_ebounds[0] - m_ebounds[2])/ (m_ebounds[1] - m_ebounds[2]);
   m_log_ratio_hi = std::log(m_ebounds[nee]/m_ebounds[nee-2]) / std::log(m_ebounds[nee-1]/m_ebounds[nee-2]);
   m_diff_ratio_hi = (m_ebounds[nee] - m_ebounds[nee-2]) / (m_ebounds[nee-1] - m_ebounds[nee-2]);

   compute_drm();
}


Drm::Drm(double ra, double dec, 
         const std::vector<double> & ebounds, 
	 const std::vector< std::vector<double> >& values)
  : m_dir(ra, dec), m_observation(0), m_npts(0), m_drm(values) {
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

   size_t nee(m_ebounds.size() - 1);

   m_log_ratio_lo = std::log(m_ebounds[0]/m_ebounds[2]) / std::log(m_ebounds[1]/m_ebounds[2]);
   m_diff_ratio_lo = (m_ebounds[0] - m_ebounds[2])/ (m_ebounds[1] - m_ebounds[2]);
   m_log_ratio_hi = std::log(m_ebounds[nee]/m_ebounds[nee-2]) / std::log(m_ebounds[nee-1]/m_ebounds[nee-2]);
   m_diff_ratio_hi = (m_ebounds[nee] - m_ebounds[nee-2]) / (m_ebounds[nee-1] - m_ebounds[nee-2]);
   compute_transpose();
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

   value = extrapolate_lo(counts);
   counts.push_front(value);

   value = extrapolate_hi(counts);
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
   if ( m_observation == 0 ) {
     throw std::runtime_error("Called compute_drm on a Drm that was build from a FITs table and does not have obsrvation data");
   }
   m_drm.clear();
   
   // EAC, Fix this code to use the binning from the livetime cube instead of sampling in cos theta.
   // Sampling in cos theta can fail for pointed observations where all of the livetime comes 
   // in a particular cos theta bin.
   compute_livetime();

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
   compute_transpose();
}


void Drm::compute_transpose() {
  size_t n_true = m_ebounds.size() - 1;
  size_t n_meas = n_true - 2;
  m_drmT.clear();
  m_drmT.resize(n_meas, std::vector<double>(n_true, 0));
  for ( size_t k(0); k < n_true; k++ ) {
    for ( size_t kp(0); kp < n_meas; kp++ ) {
       m_drmT[kp][k] = m_drm[k][kp];
    }
  }
}


void Drm::compute_livetime() {
  // This it can be done once for the entire matrix
  const ExposureCube & expcube = m_observation->expCube();
  const healpix::CosineBinner& cos_binner = expcube.get_cosine_binner(m_dir);
  size_t nmu = cos_binner.size();
  
  m_costheta_vals.resize(nmu, 0);
  m_theta_vals.resize(nmu, 0);
  m_livetime.resize(nmu, 0);

  size_t j(0);
  for ( std::vector<float>::const_iterator itrcos = cos_binner.begin(); 
	itrcos != cos_binner.end_costh(); itrcos++, j++) {
    double cos_theta = cos_binner.costheta(itrcos);
    m_costheta_vals[j] = cos_theta;
    m_theta_vals[j] = std::acos(cos_theta)*180./M_PI;
    m_livetime[j] = *itrcos;
  }

}


double Drm::
matrix_element(double etrue, double emeas_min, double emeas_max) const {
   const ResponseFunctions & resps(m_observation->respFuncs());
   const ExposureCube & expcube(m_observation->expCube());
   double met((expcube.tstart() + expcube.tstop())/2.);

   size_t nmu = m_livetime.size();
   if ( nmu == 0 ) {
     throw std::runtime_error("Drm::matrix_element() called before Drm::compute_livetime()");
   }

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

   std::vector<double> exposr(nmu, 0);
   std::vector<double> top(nmu, 0);

   for ( size_t j(0); j < nmu; j++ ) {
      double theta = m_theta_vals[j];
      double livetime = m_livetime[j];
      for (size_t i(0); i < evtTypes.size(); i++) {
         double aeff(resps.aeff(etrue, theta, phi, evtTypes[i], met));
         double edisp(resps.edisp(evtTypes[i]).integral(emeas_min, emeas_max, 
                                                        etrue, theta, phi,
                                                        met));
	 exposr[j] += aeff*livetime;
         top[j] += edisp*aeff*livetime;
      }
   }
   double numerator(::integrate(m_costheta_vals, top));
   double denominator(::integrate(m_costheta_vals, exposr));
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


double Drm::extrapolate_lo(const std::deque<double>& counts) const{
  // EAC, this can be wildly off for value near zero
  // double min_counts_value(1e-10);
  return extrapolate_lo(counts[0], counts[1]);
}

double Drm::extrapolate_lo(const double& c0, const double & c1) const{
  // EAC, this can be wildly off for value near zero
  // double min_counts_value(1e-10);
  static double min_counts_value(1e-2);

  double value(0.);
  if (c0 > min_counts_value && c1 > min_counts_value) {
    value = c1*std::exp(m_log_ratio_lo*std::log(c0/c1));
  } else {
    value = (m_diff_ratio_lo * (c0 - c1) + c1);
  }
  return value;
}

  
  
double Drm::extrapolate_hi(const std::deque<double>& counts) const{
  size_t ncc(counts.size() - 1);
  return extrapolate_hi(counts[ncc-1], counts[ncc]);
}


double Drm::extrapolate_hi(const double& c0, const double & c1) const{
  // EAC, this can be wildly off for value near zero
  // double min_counts_value(1e-10);
  static double min_counts_value(1e-2);

  double value(0.);
  if (c1 > min_counts_value && c0 > min_counts_value) {
    value = c0 * std::exp(m_log_ratio_hi * std::log(c1/c0));
  } else {
    value = (m_diff_ratio_hi*(c1 - c0) + c0);
  }
  return value;
}



Drm_Cache::Drm_Cache(const Drm* drm,
		     SourceMap & sourceMap,
		     const std::vector<double>& energies,
		     int edisp_val)
  :m_true_counts(energies.size()-1),
   m_meas_counts(energies.size()-1),
   m_xi(energies.size()-1),
   m_kref(energies.size()-1, -1),
   m_true_counts_wt(energies.size()-1),
   m_meas_counts_wt(energies.size()-1),
   m_edisp_val(drm!=0 ? edisp_val : -1 ){
  update(drm,sourceMap,energies,edisp_val);
}

Drm_Cache::Drm_Cache(const Drm_Cache& other) 
  : m_true_counts(other.m_true_counts),
    m_meas_counts(other.m_meas_counts),
    m_xi(other.m_xi),
    m_kref(other.m_kref),
    m_true_counts_wt(other.m_true_counts_wt),
    m_meas_counts_wt(other.m_meas_counts_wt),
    m_edisp_val(other.m_edisp_val){
}

void Drm_Cache::update(const Drm* drm,
		       SourceMap & sourceMap,
		       const std::vector<double>& energies, 
		       int edisp_flag){

  sourceMap.setSpectralValues(energies);  
  const std::vector<double>& npreds = sourceMap.npreds();
  const std::vector<std::vector<std::pair<double,double> > >& weighted_npreds = sourceMap.weighted_npreds();
  const std::vector<std::pair<double,double> > & spec_wts = sourceMap.specWts();
  const BinnedCountsCache& dataCache = *(sourceMap.dataCache());

  // These are the weights to be applied to the measured counts.
  size_t k(0);
  bool has_weights = sourceMap.weights() != 0;

  m_edisp_val = drm != 0 ? edisp_flag : -1;

  std::vector<double> edisp_col;
  for (k = 0; k < energies.size()-1; k++) {
    double counts(0.);
    double counts_wt(0.);
    try {
      FitUtils::npred_contribution(npreds, weighted_npreds.at(k).at(0), spec_wts, 1., k, counts, counts_wt);
    } catch (...) {
      std::cout << k << ' ' << weighted_npreds.size();
      if ( weighted_npreds.size() > k ) {
	std::cout << ' ' << weighted_npreds.at(k).size() << std::endl;
      }
    }
    m_true_counts[k] = counts;
    m_true_counts_wt[k] = counts_wt;

    if ( m_edisp_val > 0 ) {
      size_t kmin_edisp(0.);
      size_t kmax_edisp(0.);
      FitUtils::get_edisp_constants(sourceMap, dataCache, m_edisp_val, k, kmin_edisp, kmax_edisp, edisp_col);
      counts = 0.;
      counts_wt = 0.;
      FitUtils::npred_edisp(npreds, weighted_npreds.at(k), spec_wts, edisp_col, kmin_edisp, kmax_edisp, counts, counts_wt);      
    }
    m_meas_counts[k] = counts;
    m_meas_counts_wt[k] = counts_wt;   
  }
   
  if ( m_edisp_val == 0 ) {
    drm->convolve(m_true_counts, m_meas_counts);
  } 

  // EAC FIXME, Not sure why this can happen, but it is NOT a good thing.
  static bool first(true);
  for (k = 0; k < energies.size()-1; k++) {
    if ( m_meas_counts[k] < 0 ) {
      if (first) {
	first = false;
	std::cout << "Drm_Cache::update Measured counts < 0 " << sourceMap.name() << ' ' << k << ' '
		  << m_meas_counts[k] << ' ' << m_true_counts[k] << std::endl;
	for (size_t kk(0); kk < energies.size()-1; kk++) {
	  std::cout << m_true_counts[kk] << ' ';
	}
	std::cout << std::endl;
      }
      m_meas_counts[k] = 0.;
    }
  }
    

  int kref(-1);
  int kref_wt(-1);
  int idx(0);

  for (k = 0; k < energies.size()-1; k++) {
    if ( m_true_counts[k] > 0 ) {
      // Still have counts in this true energy bin, so it can 
      // be used as a reference
      m_xi[k] = m_meas_counts[k] / m_true_counts[k];
      kref = k;
      m_kref[k] = -1;
    } else {
      // Don't have counts in this true energy bin.  
      // Check to see if there are any true counts (i.e., if kref > 0)
      if ( kref >= 0 ) {
	// If so, use the kref bin
	m_xi[k] = m_meas_counts[k] / m_true_counts[kref];
	m_kref[k] = kref;
      } else {
	// No true counts, set xi to 1 and kref to -1
	m_xi[k] = 1.;
	m_kref[k] = -1.;
      }
    }
    
    if ( m_edisp_val == 0 ) {
      m_meas_counts_wt[k] = m_xi[k]*m_true_counts_wt[kref];
    }
    
    if ( m_xi[k] < 0 ) {
      std::cout << "Xi < 0 " << sourceMap.name() << ' ' << k << ' ' << kref << ' '
		<< m_xi[k] << ' ' 
		<< m_meas_counts[k] << ' ' << m_true_counts[k] << std::endl;	
    }

  }

}


size_t Drm_Cache::memory_size() const {
  size_t retVal(0);
  retVal += sizeof(*this);
  retVal += sizeof(double)*m_true_counts.capacity();
  retVal += sizeof(double)*m_meas_counts.capacity();
  retVal += sizeof(double)*m_xi.capacity();
  retVal += sizeof(int)*m_kref.capacity();
  retVal += sizeof(double)*m_true_counts_wt.capacity();
  retVal += sizeof(double)*m_meas_counts_wt.capacity();

}

} // namespace Likelihood
