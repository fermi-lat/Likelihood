/**
 * @file ExposureMap.h
 * @brief ExposureMap class declaration.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureMap.h,v 1.5 2003/04/25 18:32:18 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureMap.h,v 1.5 2003/04/25 18:32:18 jchiang Exp $
 *
 */

class ExposureMap {

public:

   ~ExposureMap() {} //{delete s_mapData;}

   static ExposureMap * instance();

   //! Read exposure map FITS file and compute the static data members.
   static void readExposureFile(std::string exposureFile);

   //! This method computes the energy-dependent coefficients for
   //! the predicted number of photons for this source.  The exposure 
   //! vector contains the integral of the exposure map times the 
   //! spatial distribution (spatialDist) over the Source Region
   //! as defined by the exposure map extent.  See 
   //! <a href="http://lheawww.gsfc.nasa.gov/~jchiang/SSC/like_3.ps>
   //! LikeMemo 3</a>, equations 20, 29, and 30 in particular.
   void integrateSpatialDist(std::vector<double> &energies, 
                             Function * spatialDist, 
                             std::vector<double> &exposure);

   void fetchRA(std::valarray<double> &ra) 
      {ra.resize(s_ra.size()); ra = s_ra;}
   void fetchDec(std::valarray<double> &dec) 
      {dec.resize(s_ra.size()); dec = s_dec;}
   void fetchEnergies(std::vector<double> &energies) {energies = s_energies;}
   void fetchExposure(std::vector< std::valarray<double> > &exposure);

   static void computeMap(const std::string &filename, double sr_radius = 30,
                          int nlong = 60, int nlat = 60, int nenergies = 10);

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

   //! write the FITS image file produced by computeMap()
   static void writeFitsFile(const std::string &filename,
                             std::vector<double> &lon,
                             std::vector<double> &lat,
                             std::vector<double> &energies,
                             std::vector< std::valarray<double> > &dataCube,
                             double ra0, double dec0);

};

} // namespace Likelihood

#endif  // ExposureMap_h
