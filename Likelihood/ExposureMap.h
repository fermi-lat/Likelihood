/**
 * @file ExposureMap.h
 * @brief ExposureMap class declaration.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureMap.h,v 1.18 2005/03/01 07:17:06 jchiang Exp $
 */

#ifndef Likelihood_ExposureMap_h
#define Likelihood_ExposureMap_h

#include <string>
#include <vector>

#include "optimizers/Function.h"

#include "Likelihood/FitsImage.h"

namespace Likelihood {

   class Observation;

/**
 * @class ExposureMap 
 *
 * @brief This Singleton class encapulates and provides exposure map
 * information, primarily for use by the DiffuseSource class for
 * integrating the response functions over the spatial distributions
 * of those Sources.
 *
 * The exposure map can be read in from an existing file (or perhaps 
 * computed ab initio given the ROI cuts and spacecraft data).
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ExposureMap.h,v 1.18 2005/03/01 07:17:06 jchiang Exp $
 *
 */

class ExposureMap {

public:

   ~ExposureMap() {}

   static ExposureMap * instance();

   /// Read exposure map FITS file and compute the static data members.
   static void readExposureFile(std::string exposureFile);

   /**
    * @brief This method computes the energy-dependent coefficients
    * for the predicted number of photons for this source.
    *
    * The exposure vector contains the integral of the exposure map
    * times the spatial distribution (spatialDist) over the Source
    * Region as defined by the exposure map extent.  See <a
    * href="http://lheawww.gsfc.nasa.gov/~jchiang/SSC/like_3.ps>
    * LikeMemo 3</a>, equations 20, 29, and 30 in particular.
    *
    * @param energies A vector of energies at which the DiffuseSource
    * spectrum is evaluated.
    * @param spatialDist A Function object that takes a SkyDirArg as 
    * its argument.
    * @param exposure A vector of exposure values characterizing the
    * DiffuseSource spectral response.
    */
   void integrateSpatialDist(const std::vector<double> &energies, 
                             optimizers::Function * spatialDist, 
                             std::vector<double> &exposure) const;

   /// Retrieve the RA of each pixel in the image plane
   void getRA(std::vector<double> & ra) const {
      ra.resize(s_ra.size());
      ra = s_ra;
   }

   /// Retrieve the Dec of each pixel in an image plane
   void getDec(std::vector<double> &dec) const {
      dec.resize(s_ra.size());
      dec = s_dec;
   }

   /// Retrieve the energies in MeV of each plane in the ExposureMap 
   /// frame stack
   void getEnergies(std::vector<double> &energies) const {
      energies = s_energies;
   }

   /**
    * @brief Compute the exposure map given the current set of
    * spacecraft data and write to a file.
    *
    * @param filename The output file
    * @param sr_radius The half of the range longitude and latitude in
    * degrees
    * @param nlong The number of pixels in the longitude coordinate
    * @param nlat The number of pixels in the latitude coordinate
    * @param nenergies The number of energies, which are in MeV.
    * These are logarithmically spaced with upper and lower bounds
    * given by the RoiCuts.
    */
   static void computeMap(std::string filename, 
                          const Observation & observation,
                          double sr_radius=30, int nlong=60, int nlat=60,
                          int nenergies=10);

   /// write the FITS image file produced by computeMap()
   static void writeFitsFile(const std::string &filename,
                             std::vector<double> &lon,
                             std::vector<double> &lat,
                             std::vector<double> &energies,
                             std::vector< std::vector<double> > &dataCube,
                             double ra0, double dec0);

protected:

   ExposureMap() {}

private:

   static ExposureMap * s_instance;

   static bool s_haveExposureMap;

   /// s_ra and s_dec are vectors of size NAXIS1*NAXIS2.
   /// Traversing these vectors in tandem yields all coordinate pairs
   /// of the image plane.
   static std::vector<double> s_ra;
   static std::vector<double> s_dec;

   /// True photon energies associated with each image plane.
   static std::vector<double> s_energies;

   /// s_exposure is a vector of size NAXIS3, corresponding to the
   /// number of true energy values identified with each plane in the
   /// exposure data cube.
   static std::vector< std::vector<double> > s_exposure;

   static FitsImage s_mapData;

};

} // namespace Likelihood

#endif  // Likelihood_ExposureMap_h
