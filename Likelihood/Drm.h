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

    /* Standard c'tor,        
       Build from a position, and observation and vector of energies  */       
    Drm(double ra, double dec, const Observation & observation, 
	const std::vector<double> & ebounds, size_t edisp_bins=0);

    /* Recreation c'tor, 
       Build from a stored Drm */
    Drm(double ra, double dec, 
	const std::vector<double> & ebounds, 
	const std::vector< std::vector<double> >& values, size_t edisp_bins=0);

    /* The Reference direction */
    const astro::SkyDir& refDir() const { return m_dir; }

    /* The underlying observation */
    inline const Observation* observation() const { return m_observation; }

    /* The number of true energy bins in the response matrix
       This is equal to energies().size() + 2*edisp_bins() -1 */
    inline size_t ntrue() const { return m_ntrue; }

    /* The number of measured energy bins in the response matrix
       This is equal to energies().size() - 1 */
    inline size_t nmeas() const { return m_nmeas; }

    /* The number of extra bins of energy to consider on either side of the spectrum */
    inline size_t edisp_bins() const { return m_edisp_bins; }

    /* The vector of energy bin edges used in the analysis */
    inline const std::vector<double>& energies() const { return m_energies; }

    /* The extended vector of energy bin edges used for energy dispersion */
    inline const std::vector<double>& full_energies() const { return m_full_energies; }
    
    /* These are all of the responses for a True energy bin       
       Note that the numbering of true and measured energy bins 
       differs by edisp_bins() */
    const std::vector<double> & row(size_t k) const {
      return m_drm.at(k);
    }

    /* These are all of the responses for a Measured Energy bin */       
    const std::vector<double> & col(size_t k) const {
      return m_drmT.at(k);
    }

    /* Convonve a true counts disribution to get a measured counts distribution

       Note: this version assumes that you don't have a better way to
       extrapolate the true_counts spectrum and uses this
       class's extraploation methods.

       This requires that the length of the true_counts vector is equal
       to the nmeas() */
    void convolve(const std::vector<double> & true_counts,
		  std::vector<double> & meas_counts) const;

    /* Compute a matrix element */
    double matrix_element(double etrue, double emeas_min, 
			  double emeas_max) const;

    /* extend the true counts distribution on the low side */
    double extrapolate_lo(const std::deque<double> & counts) const;

    /* compute the true counts distribution for the next bin on the low side */
    double extrapolate_lo(const double& c0, const double & c1) const;
  
    /* extend the true counts distribution on the high side */
    double extrapolate_hi(const std::deque<double> & counts) const;

    /* compute the true counts distribution for the next bin on the high side */
    double extrapolate_hi(const double& c0, const double & c1) const;

  protected: 
    
    /* utility function to compute the livetime as a function of theta */
    void compute_livetime();

    /* utility function to compute the transverse of the drm */
    void compute_transpose();

  private:
    
    // The reference direction
    astro::SkyDir m_dir;
    
    // This is needed to compute the martix elements
    const Observation * m_observation;
    
    // The orignal vector of energy bin edges
    std::vector<double> m_energies;

    // The full vector of energy bin edges, including the extra bins tacked on
    std::vector<double> m_full_energies;
    
    // The number of true energy bins in the DRM
    size_t m_ntrue;
    
    // The number of measured energy bins in the DRM
    size_t m_nmeas;
    
    // The number of extra bins tacked on either end of the energy range
    size_t m_edisp_bins;
    
    // These are need for the low energy counts extrapolation
    double m_log_ratio_lo;
    double m_diff_ratio_lo;
    
    // These are needed for the high energy counts extrapolation  
    double m_log_ratio_hi;
    double m_diff_ratio_hi;
   
    // These are needed for doing the integrals over the instrument frame
    std::vector< double > m_costheta_vals;
    std::vector< double > m_theta_vals;
    std::vector< double > m_livetime;

    // This matrix is indexed m_drm[kTrue][kMeasured]
    std::vector< std::vector<double> > m_drm;
    
    // This matrix is indexed m_drm[kMeasured][kTrue]
    std::vector< std::vector<double> > m_drmT;
    
    void compute_drm();
           
  };


  class Drm_Cache {
  public:
    
    /* Standard c'tor, construct from a Drm and a SourceMap */
    Drm_Cache(const Drm& drm,
	      SourceMap & sourceMap);
    
    /* Copy c'tor */
    Drm_Cache(const Drm_Cache& other);
    
    /* Clone operator, uses copy c'tor */
    Drm_Cache* clone() const {
      return new Drm_Cache(*this);
    }
    
    /* D'tor */
    virtual ~Drm_Cache(){;}
    
    /* Spectrum of true counts */
    inline const std::vector<double>& true_counts() const { return m_true_counts; }  

    /* Spectrum of measured counts */
    inline const std::vector<double>& meas_counts() const { return m_meas_counts; }  

    /* Correction factors for rescaling method */
    inline const std::vector<double>& xi() const { return m_xi; }  

    /* Reference bins for rescaling method */
    inline const std::vector<int>& kref() const { return m_kref; }
    
    /* Weighted true counts */
    inline const std::vector<double>& true_counts_wt() const { return m_true_counts_wt; }  
    
    /* Weighted measured counts */
    inline const std::vector<double>& meas_counts_wt() const { return m_meas_counts_wt; }  
    
    /* Get the correction factor and reference bin for rescaling method */
    inline double get_correction(size_t k, int& kref) const {
      kref = m_kref[k];
      return m_xi[k];
    }

    /* Get the method to use for energy dispersion 
       
       edisp_val < 0 : no energy dispersion 
       edisp_val = 0 : rescaling method
       edisp_val > 0 : standard method with edisp_val bins
    */
    inline int edisp_val() const { return m_edisp_val; }

    /* Updated the cached values */
    void update(const Drm& drm,
		SourceMap & sourceMap);
    
    /* --------------------- Debugging -------------------- */
    size_t memory_size() const;
    
  private:
    
    /* Spectrum of true counts */
    std::vector<double> m_true_counts;
    /* Spectrum of measured counts */
    std::vector<double> m_meas_counts;  
    /* Correction factors for rescaling method */
    std::vector<double> m_xi;  
    /* Reference bins for rescaling method */
    std::vector<int> m_kref;

    /* Weighted true counts */   
    std::vector<double> m_true_counts_wt;
    /* Weighted measured counts */
    std::vector<double> m_meas_counts_wt;

    /* Method and number of bins to use for energy dispersion */
    int m_edisp_val;
    
  };

} // namespace Likelihood

#endif // Likelihood_Drm_h
