/**
 * @file ExposureMap.cxx
 * @brief Implementation for the Singleton ExposureMap class. This
 * class encapsulates exposure map information and makes it available
 * for use (primarily) by the DiffuseSource class.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ExposureMap.cxx,v 1.1 2003/03/25 23:22:03 jchiang Exp $
 */

#include "Likelihood/SkyDirArg.h"
#include "Likelihood/ExposureMap.h"

namespace Likelihood {

ExposureMap * ExposureMap::s_instance = 0;

std::vector<double> ExposureMap::s_energies;
std::valarray<double> ExposureMap::s_ra;
std::valarray<double> ExposureMap::s_dec;
std::vector< std::valarray<double> > ExposureMap::s_exposure;
FitsImage *ExposureMap::s_mapData = 0;

void ExposureMap::readExposureFile(std::string exposureFile) {

   s_mapData = new FitsImage(exposureFile);

   s_mapData->fetchCelestialArrays(s_ra, s_dec);

// Fetch the energy axis abscissa points. Here we assume that the
// exposure map has at least two image planes, and that the energies
// are along the third dimension so we set naxis = 2.
   int naxis = 2;
   s_mapData->fetchAxisVector(naxis, s_energies);

// pixel solid angles
   std::valarray<double> solidAngles;
   s_mapData->fetchSolidAngles(solidAngles);

// Fill the vector of the planes of s_exposure for each true photon
// energy.
   std::valarray<double> exposure;
   s_mapData->fetchImageData(exposure);

   int indx = 0;
   int npixels = solidAngles.size();
   for (unsigned int k = 0; k < s_energies.size(); k++) {
      std::valarray<double> expArray(npixels);
      for (int sa_indx = 0; sa_indx < npixels; sa_indx++) {
         expArray[sa_indx] = solidAngles[sa_indx]*exposure[indx];
         indx++;
      }
      s_exposure.push_back(expArray);
   }
}

void ExposureMap::integrateSpatialDist(std::vector<double> &energies,
                                       Function * spatialDist,
                                       std::vector<double> &exposure) {

// Fetch the exposure multiplied by the solid angle of the associated
// pixel
   std::vector< std::valarray<double> > my_exposure;
   fetchExposure(my_exposure);

   exposure.clear();
   exposure.reserve(energies.size());
   for (unsigned int k = 0; k < energies.size(); k++) {
      double srcExposure = 0;
// Find the index kk (of the energy array that describes the exposure
// map data) that will be used for interpolating at energies[k]; here
// we assume the ExposureMap energies are logarithmically spaced.
      std::vector<double>::const_iterator iterE;
      if (energies[k] < s_energies[0]) {
         iterE = s_energies.begin() + 1;
      } else if (energies[k] >= *(s_energies.end() - 1)) {
         iterE = s_energies.end() - 1;
      } else {
         iterE = upper_bound(s_energies.begin(), s_energies.end(), 
                             energies[k]);
      }
      int kk = iterE - s_energies.begin();
      for (unsigned int j = 0; j < s_ra.size(); j++) {
         double expsr = log(energies[k]/(*(iterE-1)))/log(*iterE/(*(iterE-1)));
         expsr = expsr*(my_exposure[kk][j] - my_exposure[kk-1][j])
            + my_exposure[kk-1][j];

         astro::SkyDir skyDir(s_ra[j], s_dec[j]);
         SkyDirArg dir(skyDir);
         srcExposure += expsr*(*spatialDist)(dir);
      }
      exposure.push_back(srcExposure);
   }
}

void ExposureMap::fetchExposure(std::vector< std::valarray<double> > 
                                &exposure) {
   if (!exposure.empty()) exposure.clear();

   exposure.reserve(s_exposure.size());
   for (unsigned int i = 0; i < s_exposure.size(); i++) {
      exposure.push_back(s_exposure[i]);
   }
}

ExposureMap * ExposureMap::instance() {
   if (s_instance == 0) {
      s_instance = new ExposureMap();
   }
   return s_instance;
}

} // namespace Likelihood
