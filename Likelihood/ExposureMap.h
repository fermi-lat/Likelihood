/**
 * @file ExposureMap.h
 * @brief ExposureMap class declaration.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureMap.h,v 1.3 2003/03/25 23:22:02 jchiang Exp $
 */

#ifndef ExposureMap_h
#define ExposureMap_h

#include <vector>
#include <string>
#include <valarray>
#include "Likelihood/Function.h"
#include "Likelihood/FitsImage.h"

namespace Likelihood {

/**
 * @class ExposureMap 
 *
 * @brief This is a Singleton class that encapulates and provides
 * exposure map information, primarily for use by the DiffuseSource
 * class for integrating the response functions over the spatial
 * distributions of those Sources.
 *
 * The exposure map can be read in from an existing file (or perhaps 
 * computed ab initio given the ROI cuts and spacecraft data).
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureMap.h,v 1.3 2003/03/25 23:22:02 jchiang Exp $
 *
 */

class ExposureMap {

public:

   ~ExposureMap() {} //{delete s_mapData;}

   static ExposureMap * instance();

   //! Read exposure map FITS file and compute the static data members.
   static void readExposureFile(std::string exposureFile);

   void integrateSpatialDist(std::vector<double> &energies, 
                             Function * spatialDist, 
                             std::vector<double> &exposure);

   void fetchRA(std::valarray<double> &ra) 
      {ra.resize(s_ra.size()); ra = s_ra;}
   void fetchDec(std::valarray<double> &dec) 
      {dec.resize(s_ra.size()); dec = s_dec;}
   void fetchEnergies(std::vector<double> &energies) {energies = s_energies;}
   void fetchExposure(std::vector< std::valarray<double> > &exposure);

protected:

   ExposureMap() {}

private:

   static ExposureMap * s_instance;

   //! s_ra and s_dec are valarrays of size NAXIS1*NAXIS2.
   //! Traversing these valarrays in tandem yields all coordinate pairs
   //! of the image plane.
   static std::valarray<double> s_ra;
   static std::valarray<double> s_dec;

   //! True photon energies associated with each image plane.
   static std::vector<double> s_energies;

   //! s_exposure is a vector of size NAXIS3, corresponding to the
   //! number of true energy values identified with each plane in the
   //! exposure data cube.
   static std::vector< std::valarray<double> > s_exposure;

   static FitsImage *s_mapData;
};

} // namespace Likelihood

#endif  // ExposureMap_h
