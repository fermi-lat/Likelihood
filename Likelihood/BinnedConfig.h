/**
 * @file BinnedConfig.h
 * @brief Small helper classes to deal with BinnedLikelihood configuration
 * @author E. Charles

 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FitUtils.h,v 1.14 2016/07/11 23:44:04 mdwood Exp $
 */

#ifndef Likelihood_BinnedConfig_h
#define Likelihood_BinnedConfig_h


namespace Likelihood {

  /* This is a small helper class to encapsulate information about
     how to treat the PSF Integration */
  class PsfIntegConfig {

  public:
    
    typedef enum { 
      //! Addaptive integration
      adaptive = 0,
      //! Annular integration
      annular = 1,
      //! Simple integration using pixel centers
      pixel_center = 2 } PSFIntegType;

    /* Get value of configuration paramters from the enviroment */
    static void get_envars(PSFIntegType& integ_type,
			   double& estimatorFtol,
			   double& estimatorPeakTh);

  public:

    /* Standard c'tor, all arguments can take default value */
    PsfIntegConfig(bool applyPsfCorrections = true,
		   bool performConvolution = true,
		   bool resample = true,
		   double resamp_factor = 2,
		   double minbinsz = 0.1,
		   PSFIntegType integ_type = adaptive,
		   double psfEstimatorFtol = 1e-3,
		   double psfEstimatorPeakTh = 1e-6,
		   bool verbose = true)
      :m_applyPsfCorrections(applyPsfCorrections),
       m_performConvolution(performConvolution),
       m_resample(resample),
       m_resamp_factor(resamp_factor),
       m_minbinsz(minbinsz),
       m_integ_type(integ_type),
       m_psfEstimatorFtol(psfEstimatorFtol),
       m_psfEstimatorPeakTh(psfEstimatorPeakTh),
       m_verbose(verbose){
    }

    /* Copy c'tor */
    PsfIntegConfig(const PsfIntegConfig& other)
      :m_applyPsfCorrections(other.m_applyPsfCorrections),
       m_performConvolution(other.m_performConvolution),
       m_resample(other.m_resample),
       m_resamp_factor(other.m_resamp_factor),
       m_minbinsz(other.m_minbinsz),
       m_integ_type(other.m_integ_type),
       m_psfEstimatorFtol(other.m_psfEstimatorFtol),
       m_psfEstimatorPeakTh(other.m_psfEstimatorPeakTh),
       m_verbose(other.m_verbose){
    }

    /* D'tor, trivial */
    ~PsfIntegConfig(){;}

    /* Simple access to all members */
    inline bool& applyPsfCorrections() { return m_applyPsfCorrections; }
    inline const bool& applyPsfCorrections() const { return m_applyPsfCorrections; }

    inline bool& performConvolution() { return m_performConvolution; }
    inline const bool& performConvolution() const { return m_performConvolution; }

    inline bool& resample() { return m_resample; }  
    inline const bool& resample() const { return m_resample; }  

    inline double& resamp_factor() { return m_resamp_factor; } 
    inline const double& resamp_factor() const { return m_resamp_factor; } 

    inline double& minbinsz() { return m_minbinsz; }
    inline const double& minbinsz() const { return m_minbinsz; }
 
    inline PSFIntegType& integ_type() { return m_integ_type; }
    inline const PSFIntegType& integ_type() const { return m_integ_type; }
    
    inline double& psfEstimatorFtol() { return m_psfEstimatorFtol; }
    inline const double& psfEstimatorFtol() const { return m_psfEstimatorFtol; }
    
    inline double& psfEstimatorPeakTh() { return m_psfEstimatorPeakTh; }
    inline const double& psfEstimatorPeakTh() const { return m_psfEstimatorPeakTh; }

    inline bool& verbose() { return m_verbose; }
    inline const bool& verbose() const { return m_verbose; }

  private:
    
    bool m_applyPsfCorrections;  //! Apply PSF integral corrections
    bool m_performConvolution;   //! Do PSF convolution when making SourceMaps
    bool m_resample;             //! turn on Resample when projecting Diffuse source maps
    double m_resamp_factor;      //! Resample factor for projecting Diffuse source maps
    double m_minbinsz;           //! Minimum pixel size for rebinning fine maps 
    PSFIntegType m_integ_type;   //! Integration type: "adaptive" or "annular" or "pixel_center"
    double m_psfEstimatorFtol;   //! Fractional tolerance on adaptive PSF Integration
    double m_psfEstimatorPeakTh; //! Peak threshold on adaptive PSF Integration
    bool m_verbose;              //! Turn on verbose output

  };
    

  /* This is a small helper class to encapsulate information about
     the BinnedLikelihood configuration */
  class BinnedLikeConfig {

  public:

    /* Standard c'tor, all arguments can take default value */
    BinnedLikeConfig(bool computePointSources = true,		     
		     bool applyPsfCorrections = true,
		     bool performConvolution = true,
		     bool resample = true,
		     double resamp_factor = 2,
		     double minbinsz = 0.1,
		     PsfIntegConfig::PSFIntegType integ_type = adaptive,
		     double psfEstimatorFtol = 1e-3,
		     double psfEstimatorPeakTh = 1e-6,
		     bool verbose = true,
		     bool use_edisp = false, 
		     bool use_single_fixed_map = true) 
      :m_computePointSources(computePointSources),       
       m_psf_integ_config(applyPsfCorrections,performConvolution,resample,resamp_factor,minbinsz,
			  integ_type,psfEstimatorFtol,psfEstimatorPeakTh,verbose),
       m_use_edisp(use_edisp),
       m_use_single_fixed_map(use_single_fixed_map) { 
    }
    
    BinnedLikeConfig(const BinnedLikeConfig& other)
      :m_computePointSources(other.m_computePointSources),
       m_psf_integ_config(other.m_psf_integ_config),
       m_use_edisp(other.m_use_edisp),
       m_use_single_fixed_map(other.m_use_single_fixed_map) { 
    }
    
    inline PsfIntegConfig& psf_integ_config() { return m_psf_integ_config; }
    inline const PsfIntegConfig& psf_integ_config() const { return m_psf_integ_config; }

    inline bool& computePointSources() { return m_computePointSources; } 
    inline bool& use_edisp() { return m_use_edisp; }
    inline bool& use_single_fixed_map() { return m_use_single_fixed_map; }
    
    inline const bool& computePointSources() const { return m_computePointSources; } 
    inline const bool& use_edisp() const { return m_use_edisp; }
    inline const bool& use_single_fixed_map() const { return m_use_single_fixed_map; }
    
  private:
    
    bool m_computePointSources;  //! Pre-compute and store SourceMaps for point sources
    PsfIntegConfig m_psf_integ_config; //! Parameters for PSF integration
    bool m_use_edisp;            //! Use energy dispersion
    bool m_use_single_fixed_map; //! Use a single model for all fixed components

  };

} // namespace Likelihood

#endif // Likelihood_PSFUtils_h
