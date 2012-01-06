/**
 * @file ExposureMap.h
 * @brief ExposureMap class declaration.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/ExposureMap.h,v 1.22 2006/10/03 20:52:09 jchiang Exp $
 */

#ifndef Likelihood_ExposureMap_h
#define Likelihood_ExposureMap_h

#include <string>
#include <vector>

#include "optimizers/Function.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

   class Observation;
   class WcsMap2;

/**
 * @class ExposureMap 
 *
 * @brief This class encapulates and provides exposure map
 * information, primarily for use by the DiffuseSource class for
 * integrating the response functions over the spatial distributions
 * of those Sources.
 *
 * The exposure map can be read in from an existing file (or perhaps 
 * computed ab initio given the ROI cuts and spacecraft data).
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/ExposureMap.h,v 1.22 2006/10/03 20:52:09 jchiang Exp $
 *
 */

class ExposureMap {

public:

   ExposureMap() : m_wcsmap(0), m_haveExposureMap(false) {}

   ~ExposureMap();

   /// Read exposure map FITS file and compute the static data members.
   void readExposureFile(std::string exposureFile);

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
    *        spectrum is evaluated.
    * @param spatialDist A Function object that takes a SkyDirArg as 
    *        its argument.
    * @param exposure A vector of exposure values characterizing the
    *        DiffuseSource spectral response.
    */
   void integrateSpatialDist(const std::vector<double> &energies, 
                             optimizers::Function * spatialDist, 
                             std::vector<double> &exposure) const;

   /// Retrieve the RA of each pixel in the image plane
   void getRA(std::vector<double> & ra) const {
      ra.resize(m_ra.size());
      ra = m_ra;
   }

   /// Retrieve the Dec of each pixel in an image plane
   void getDec(std::vector<double> &dec) const {
      dec.resize(m_ra.size());
      dec = m_dec;
   }

   /// Retrieve the energies in MeV of each plane in the ExposureMap 
   /// frame stack
   void getEnergies(std::vector<double> &energies) const {
      energies = m_energies;
   }

   bool haveMap() const {
      return m_haveExposureMap;
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
    * @param bool compute_submap If true, then compute just for the
    * range of longitude and latitudes pixels given by the following
    * parameters
    * @param nlongmin
    * @param nlongmax
    * @param nlatmin
    * @param nlatmax
    */
   void computeMap(std::string filename, 
                   const Observation & observation,
                   double sr_radius=30, int nlong=60, int nlat=60,
                   int nenergies=10, bool compute_submap=false,
                   int nlongmin=0, int nlongmax=0, int nlatmin=0, 
                   int nlatmax=0);

   static void readEnergyExtension(const std::string & filename,
                                   std::vector<double> & energies);

   double operator()(const astro::SkyDir & dir, double energy) const;

   double operator()(const astro::SkyDir & dir, int k) const;

   bool withinMapRadius(const astro::SkyDir & dir) const;

private:

   WcsMap2 * m_wcsmap;

   bool m_haveExposureMap;

   /// m_ra and m_dec are vectors of size NAXIS1*NAXIS2.
   /// Traversing these vectors in tandem yields all coordinate pairs
   /// of the image plane.
   std::vector<double> m_ra;
   std::vector<double> m_dec;
   std::vector<double> m_solidAngles;

   /// True photon energies associated with each image plane.
   std::vector<double> m_energies;

   /// m_exposure is a vector of size NAXIS3, corresponding to the
   /// number of true energy values identified with each plane in the
   /// exposure data cube.
   std::vector< std::vector<double> > m_exposure;

   void writeFitsFile(const std::string & filename,
                      const std::vector<long> & naxes,
                      double * crpix, double * crval, double * cdelt,
                      const std::vector<double> & energies, 
                      const std::vector<float> & expMap);
};

} // namespace Likelihood

#endif  // Likelihood_ExposureMap_h
