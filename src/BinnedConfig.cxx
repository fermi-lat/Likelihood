/**
 * @file BinnedConfig.cxx
 * @brief Classes to store configuration for BinnedLikelihood fitting
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/BinnedConfig.cxx,v 1.3 2016/09/13 19:26:23 echarles Exp $
 */

#include "Likelihood/BinnedConfig.h"
#include <cstdlib>

namespace Likelihood {

  void BinnedLikeConfig::get_envars(PsfIntegConfig::PSFIntegType& estimatorMethod,
				    double& estimatorFtol,
				    double& estimatorPeakTh,
				    bool& use_edisp,
				    bool& use_linear_quadrature,
				    bool& save_all_srcmaps) {
         
    if(::getenv("USE_ADAPTIVE_PSF_ESTIMATOR")) {
      estimatorMethod = PsfIntegConfig::adaptive;
    } else if(::getenv("USE_ANNULAR_PSF_ESTIMATOR")) {
      estimatorMethod = PsfIntegConfig::annular;
    } else if(::getenv("USE_OLD_PSF_ESTIMATOR")) {
      estimatorMethod = PsfIntegConfig::pixel_center;
    } 
    
    if (::getenv("PSF_ADAPTIVE_ESTIMATOR_FTOL"))
      estimatorFtol = atof(::getenv("PSF_ADAPTIVE_ESTIMATOR_FTOL"));
    
    if (::getenv("PSF_ADAPTIVE_ESTIMATOR_PEAK_TH"))
      estimatorPeakTh  = atof(::getenv("PSF_ADAPTIVE_ESTIMATOR_PEAK_TH"));    

    if (::getenv("USE_BL_EDISP") ) {
      use_edisp = true;
    } 

    if (::getenv("USE_LINEAR_QUADRATURE") ) {
      use_linear_quadrature = true;
    }

    if (::getenv("SAVE_ALL_SRCMAPS") ) {
      save_all_srcmaps = true;
    }

  }
 
} // namespace Likelihood
