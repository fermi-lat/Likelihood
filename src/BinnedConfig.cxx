/**
 * @file BinnedConfig.cxx
 * @brief Classes to store configuration for BinnedLikelihood fitting
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FitUtils.cxx,v 1.22 2016/07/11 23:44:01 mdwood Exp $
 */

#include "Likelihood/BinnedConfig.h"
#include <cstdlib>

namespace Likelihood {

  void PsfIntegConfig::get_envars(PSFIntegType& estimatorMethod,
				  double& estimatorFtol,
				  double& estimatorPeakTh) {
         
    if(::getenv("USE_ADAPTIVE_PSF_ESTIMATOR")) {
      estimatorMethod = adaptive;
    } else if(::getenv("USE_ANNULAR_PSF_ESTIMATOR")) {
      estimatorMethod = annular;
    } else if(::getenv("USE_OLD_PSF_ESTIMATOR")) {
      estimatorMethod = pixel_center;
    } else {
      estimatorMethod = adaptive;
    }
    
    if (::getenv("PSF_ADAPTIVE_ESTIMATOR_FTOL"))
      estimatorFtol = atof(::getenv("PSF_ADAPTIVE_ESTIMATOR_FTOL"));
    
    if (::getenv("PSF_ADAPTIVE_ESTIMATOR_PEAK_TH"))
      estimatorPeakTh  = atof(::getenv("PSF_ADAPTIVE_ESTIMATOR_PEAK_TH"));    
  }
 
} // namespace Likelihood
