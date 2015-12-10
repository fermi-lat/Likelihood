/**
 * @file Convolve.h
 * @brief Functions to perform convolutions of HEALPix maps
 * @author E. Charles
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/ConvolveHealpix.h,v 1.3 2015/11/30 19:38:30 echarles Exp $
 */

#ifndef Likelihood_ConvolveHealpix_h
#define Likelihood_ConvolveHealpix_h

#include "healpix_map.h"
#include "Likelihood/MeanPsf.h"

namespace Likelihood {

  namespace ConvolveHealpix {
    
    void fillMapWithPSF_pole(const MeanPsf & psf, 
			     const double& energy, 
			     Healpix_Map<float>& outMap);
    
    void fillMapWithPSF_refDir(const MeanPsf& psf, 
			       const double& energy, 			       
			       const astro::SkyDir& refDir,
			       bool use_lb,
			       Healpix_Map<float>& outMap);

    void convolve(const Healpix_Map<float>& inMap,
		  const Healpix_Map<float>& psf,
		  Healpix_Map<float>& outMap);

    double psfMinPixSize(double energy);
    
  };

} // namespace Likelihood

#endif // Likelihood_ConvolveHealpix_h
