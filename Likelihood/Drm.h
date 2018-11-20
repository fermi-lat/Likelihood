/**
 * @file Drm.h
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Drm.h,v 1.13 2017/08/17 23:43:54 echarles Exp $
 */

#ifndef Likelihood_Drm_h
#define Likelihood_Drm_h

#include <deque>
#include <vector>

#include "astro/SkyDir.h"

namespace Likelihood {

class Observation;
class SourceMap;
class Source;

class Drm {

public:
   
   Drm(double ra, double dec, const Observation & observation, 
       const std::vector<double> & ebounds, size_t npts=30);

   void convolve(const std::vector<double> & true_counts,
                 std::vector<double> & meas_counts) const;

   const std::vector<double> & row(size_t k) const {
      return m_drm.at(k);
   }

   const std::vector<double> & col(size_t k) const {
     return m_drmT.at(k);
   }
      
  inline const Observation& observation() const { return m_observation; }

  double matrix_element(double etrue, double emeas_min, 
			double emeas_max) const;

  double extrapolate_lo(const std::deque<double> & counts) const;

  double extrapolate_lo(const double& c0, const double & c1) const;
  
  double extrapolate_hi(const std::deque<double> & counts) const;

  double extrapolate_hi(const double& c0, const double & c1) const;


protected: 

  void compute_livetime();

  void compute_transpose();

private:

   astro::SkyDir m_dir;
   const Observation & m_observation;
   std::deque<double> m_ebounds;

   double m_log_ratio_lo;
   double m_diff_ratio_lo;
   
   double m_log_ratio_hi;
   double m_diff_ratio_hi;
   
   size_t m_npts;

   std::vector< double > m_costheta_vals;
   std::vector< double > m_theta_vals;
   std::vector< double > m_livetime;

   std::vector< std::vector<double> > m_drm;
   std::vector< std::vector<double> > m_drmT;

   void compute_drm();
  
    
                         
};


class Drm_Cache {
public:

  Drm_Cache(const Drm* drm,
	    SourceMap & sourceMap,
	    const std::vector<double>& energies,
	    int edisp_val);
  
  Drm_Cache(const Drm_Cache& other);

  Drm_Cache* clone() const {
    return new Drm_Cache(*this);
  }

  virtual ~Drm_Cache(){;}

  void update(const Drm* drm,
	      SourceMap & sourceMap,
	      const std::vector<double>& energies,
	      int edisp_flag);

  inline double get_correction(size_t k, int& kref, bool weight=false) const {
    kref = m_kref[k];
    return weight ? m_xi_wt[k] : m_xi[k];
  }

  inline const std::vector<double>& true_counts() const { return m_true_counts; }  
  inline const std::vector<double>& meas_counts() const { return m_meas_counts; }  
  inline const std::vector<double>& xi() const { return m_xi; }  
  inline const std::vector<double>& xi_wt() const { return m_xi_wt; }  

  inline const std::vector<double>& true_counts_wt() const { return m_true_counts_wt; }  
  inline const std::vector<double>& meas_counts_wt() const { return m_meas_counts_wt; }  

  inline const std::vector<int> kref() const { return m_kref; }

  inline int edisp_val() const { return m_edisp_val; }

  /* --------------------- Debugging -------------------- */
  size_t memory_size() const;

private:
  
  std::vector<double> m_true_counts;
  std::vector<double> m_meas_counts;  
  std::vector<double> m_xi;  
  std::vector<double> m_xi_wt;  
  std::vector<int> m_kref;

  std::vector<double> m_true_counts_wt;
  std::vector<double> m_meas_counts_wt;

  int m_edisp_val;

};



} // namespace Likelihood

#endif // Likelihood_Drm_h
