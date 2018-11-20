/**
 * @file BinnedConfig.h
 * @brief Small helper classes to deal with BinnedLikelihood configuration
 * @author E. Charles

 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedConfig.h,v 1.6 2016/10/20 23:07:32 echarles Exp $
 */

#ifndef Likelihood_BinnedConfig_h
#define Likelihood_BinnedConfig_h


namespace Likelihood {

  class BinnedLikeConfig;

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


  public:

    /* Standard c'tor, all arguments can take default value */
    PsfIntegConfig(bool applyPsfCorrections = true,
		   bool performConvolution = true,
		   bool resample = true,
		   double resamp_factor = 2,
		   double minbinsz = 0.1,
		   PSFIntegType integ_type = PsfIntegConfig::adaptive,
		   double psfEstimatorFtol = 1e-3,
		   double psfEstimatorPeakTh = 1e-6,
		   bool verbose = true,
		   bool use_single_psf = false)
      :m_applyPsfCorrections(applyPsfCorrections),
       m_performConvolution(performConvolution),
       m_resample(resample),
       m_resamp_factor(resamp_factor),
       m_minbinsz(minbinsz),
       m_integ_type(integ_type),
       m_psfEstimatorFtol(psfEstimatorFtol),
       m_psfEstimatorPeakTh(psfEstimatorPeakTh),
       m_verbose(verbose),
       m_use_single_psf(use_single_psf){
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
       m_verbose(other.m_verbose),
       m_use_single_psf(other.m_use_single_psf){
    }

    /* D'tor, trivial */
    ~PsfIntegConfig(){;}

    /* Simple access to all members */
    inline void set_applyPsfCorrections(bool val) { m_applyPsfCorrections = val; }
    inline bool applyPsfCorrections() const { return m_applyPsfCorrections; }

    inline void set_performConvolution(bool val) { m_performConvolution = val; }
    inline bool performConvolution() const { return m_performConvolution; }

    inline void set_resample(bool val) { m_resample = val; }  
    inline bool resample() const { return m_resample; }  

    inline void set_resamp_factor(double val) { m_resamp_factor = val; } 
    inline double resamp_factor() const { return m_resamp_factor; } 

    inline void minbinsz(double val) { m_minbinsz = val; }
    inline double minbinsz() const { return m_minbinsz; }
 
    inline void set_integ_type(PSFIntegType val) { m_integ_type = val; }
    inline PSFIntegType integ_type() const { return m_integ_type; }
    
    inline void set_psfEstimatorFtol(double val) { m_psfEstimatorFtol = val; }
    inline double psfEstimatorFtol() const { return m_psfEstimatorFtol; }
    
    inline void set_psfEstimatorPeakTh(double val) { m_psfEstimatorPeakTh = val; }
    inline double psfEstimatorPeakTh() const { return m_psfEstimatorPeakTh; }

    inline void set_verbose(bool val) { m_verbose = val; }
    inline bool verbose() const { return m_verbose; }

    inline void set_use_single_psf(bool val) { m_use_single_psf = val; }
    inline bool use_single_psf() const { return m_use_single_psf; }

  private:
    
    friend class BinnedLikeConfig;
    
    bool m_applyPsfCorrections;  //! Apply PSF integral corrections
    bool m_performConvolution;   //! Do PSF convolution when making SourceMaps
    bool m_resample;             //! turn on Resample when projecting Diffuse source maps
    double m_resamp_factor;      //! Resample factor for projecting Diffuse source maps
    double m_minbinsz;           //! Minimum pixel size for rebinning fine maps 
    PSFIntegType m_integ_type;   //! Integration type: "adaptive" or "annular" or "pixel_center"
    double m_psfEstimatorFtol;   //! Fractional tolerance on adaptive PSF Integration
    double m_psfEstimatorPeakTh; //! Peak threshold on adaptive PSF Integration
    bool m_verbose;              //! Turn on verbose output
    bool m_use_single_psf;       //! Use a single PSF for all sources

  };
    

  /* This is a small helper class to encapsulate information about
     the BinnedLikelihood configuration */
  class BinnedLikeConfig {

  public:

    /* Get value of configuration paramters from the enviroment 

       Note that if the ENV vars are not set, this will _NOT_ overwrite
       the input values.
     */
    static void get_envars(PsfIntegConfig::PSFIntegType& integ_type,
			   double& estimatorFtol,
			   double& estimatorPeakTh,
			   int& edisp_val,
			   bool& use_linear_quadrature,
			   bool& save_all_srcmaps,
			   bool& use_single_psf);

  public:

    /* Standard c'tor, all arguments can take default value */
    BinnedLikeConfig(bool computePointSources = true,		     
		     bool applyPsfCorrections = true,
		     bool performConvolution = true,
		     bool resample = true,
		     double resamp_factor = 2,
		     double minbinsz = 0.1,
		     PsfIntegConfig::PSFIntegType integ_type = PsfIntegConfig::adaptive,
		     double psfEstimatorFtol = 1e-3,
		     double psfEstimatorPeakTh = 1e-6,
		     bool verbose = true,
		     int edisp_val = -1, 
		     bool use_single_fixed_map = true,
		     bool use_linear_quadrature = false,
		     bool save_all_srcmaps = false,
		     bool use_single_psf = false,
		     bool load_existing_srcmaps = true)
      :m_computePointSources(computePointSources),       
       m_psf_integ_config(applyPsfCorrections,performConvolution,resample,resamp_factor,minbinsz,
			  integ_type,psfEstimatorFtol,psfEstimatorPeakTh,verbose,use_single_psf),
       m_edisp_val(edisp_val),
       m_use_single_fixed_map(use_single_fixed_map),
       m_use_linear_quadrature(use_linear_quadrature),
       m_save_all_srcmaps(save_all_srcmaps),
       m_load_existing_srcmaps(load_existing_srcmaps) {
      get_envars(m_psf_integ_config.m_integ_type,
		 m_psf_integ_config.m_psfEstimatorFtol,
		 m_psf_integ_config.m_psfEstimatorPeakTh,
		 m_edisp_val,
		 m_use_linear_quadrature,
		 m_save_all_srcmaps,
		 m_psf_integ_config.m_use_single_psf);
    }
    
    BinnedLikeConfig(const BinnedLikeConfig& other)
      :m_computePointSources(other.m_computePointSources),
       m_psf_integ_config(other.m_psf_integ_config),
       m_edisp_val(other.m_edisp_val),
       m_use_single_fixed_map(other.m_use_single_fixed_map),
       m_use_linear_quadrature(other.m_use_linear_quadrature),
       m_save_all_srcmaps(other.m_save_all_srcmaps),
       m_load_existing_srcmaps(other.m_load_existing_srcmaps){
    }
    
    inline PsfIntegConfig& psf_integ_config() { return m_psf_integ_config; }
    inline const PsfIntegConfig& psf_integ_config() const { return m_psf_integ_config; }

    inline void set_computePointSources(bool val) {  m_computePointSources = val; } 
    inline void set_edisp_val(int val) {  m_edisp_val = val; }
    inline void set_use_single_fixed_map(bool val) {  m_use_single_fixed_map = val; }
    inline void set_use_linear_quadrature(bool val) {  m_use_linear_quadrature = val; }
    inline void set_save_all_srcmaps(bool val) {  m_save_all_srcmaps = val; }
    inline void set_load_existing_srcmaps(bool val) {  m_load_existing_srcmaps = val; }
   
    inline bool computePointSources() const { return m_computePointSources; } 
    inline int edisp_val() const { return m_edisp_val; }
    inline bool use_single_fixed_map() const { return m_use_single_fixed_map; }
    inline bool use_linear_quadrature() const { return m_use_linear_quadrature; }
    inline bool save_all_srcmaps() const { return m_save_all_srcmaps; }
    inline bool load_existing_srcmaps() const { return m_load_existing_srcmaps; }

  private:
    
    bool m_computePointSources;    //! Pre-compute and store SourceMaps for point sources
    PsfIntegConfig m_psf_integ_config; //! Parameters for PSF integration
    int m_edisp_val;               //! Use energy dispersion
    bool m_use_single_fixed_map;   //! Use a single model for all fixed components
    bool m_use_linear_quadrature;  //! Use linear quadrature for counts integration
    bool m_save_all_srcmaps;       //! Save the source maps for all sources
    bool m_load_existing_srcmaps;  //! Load existing source maps from the srcmaps file

  };

} // namespace Likelihood

#endif // Likelihood_PSFUtils_h
