/**
 * @file ExposureMap.cxx
 * @brief Implementation for the Singleton ExposureMap class. This
 * class encapsulates exposure map information and makes it available
 * for use (primarily) by the DiffuseSource class.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ExposureMap.cxx,v 1.7 2003/05/20 23:50:15 jchiang Exp $
 */

#include "Likelihood/SkyDirArg.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/RoiCuts.h"

//#define HAVE_CCFITS
#ifdef HAVE_CCFITS
#include <CCfits>
#endif

#include <utility>
#include <algorithm>

namespace Likelihood {

ExposureMap * ExposureMap::s_instance = 0;
bool ExposureMap::s_haveExposureMap = false;
std::vector<double> ExposureMap::s_energies;
std::valarray<double> ExposureMap::s_ra;
std::valarray<double> ExposureMap::s_dec;
std::vector< std::valarray<double> > ExposureMap::s_exposure;
FitsImage ExposureMap::s_mapData;

void ExposureMap::readExposureFile(std::string exposureFile) {

   s_mapData = FitsImage(exposureFile);

   s_mapData.fetchCelestialArrays(s_ra, s_dec);

// Fetch the energy axis abscissa points. Here we assume that the
// exposure map has at least two image planes, and that the energies
// are along the third dimension so we set naxis = 2.
   int naxis = 2;
   s_mapData.fetchAxisVector(naxis, s_energies);

// pixel solid angles
   std::valarray<double> solidAngles;
   s_mapData.fetchSolidAngles(solidAngles);

// Fill the vector of the planes of s_exposure for each true photon
// energy.
   s_exposure.clear();
   std::valarray<double> exposure;
   s_mapData.fetchImageData(exposure);

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
   s_haveExposureMap = true;
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
         iterE = std::upper_bound(s_energies.begin(), s_energies.end(), 
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
   if (s_haveExposureMap) {
      return s_instance;
   } else {
      return 0;
   }
}

void ExposureMap::computeMap(const std::string &filename, double sr_radius,
                             int nlon, int nlat, int nenergies) {
   RoiCuts *roi_cuts = RoiCuts::instance();
   std::pair<astro::SkyDir, double> roi = roi_cuts->getExtractionRegion();
   astro::SkyDir roiCenter = roi.first;
   FitsImage::EquinoxRotation eqRot(roiCenter.ra(), roiCenter.dec());

   double lonstep = 2.*sr_radius/(nlon-1);
   double latstep = 2.*sr_radius/(nlat-1);
   std::vector<double> lon;
   for (int i = 0; i < nlon; i++) lon.push_back(lonstep*i - sr_radius);
   std::vector<double> lat;
   for (int j = 0; j < nlat; j++) lat.push_back(latstep*j - sr_radius);

   std::pair<double, double> elims = roi_cuts->getEnergyCuts();
   double estep = log(elims.second/elims.first)/(nenergies - 1);
   std::vector<double> energies;
   for (int k = 0; k < nenergies; k++) {
      energies.push_back(elims.first*exp(estep*k));
   }

   std::vector< std::valarray<double> > exposureCube;
   exposureCube.resize(nenergies);
   for (int k = 0; k < nenergies; k++)
      exposureCube[k].resize(nlon*nlat);

   std::cerr << "Computing the ExposureMap";
   int indx = 0;
   for (int j = 0; j < nlat; j++) {
      for (int i = 0; i < nlon; i++) {
         if ((indx % ((nlon*nlat)/20)) == 0) std::cerr << ".";
         astro::SkyDir indir(lon[i], lat[j]);
         astro::SkyDir dir;
         eqRot.do_rotation(indir, dir);
         PointSource ptsrc;
         bool updateExposure = false;
         ptsrc.setDir(dir.ra(), dir.dec(), updateExposure);
         std::vector<double> exposure;
         int verbose = 0;
         ptsrc.computeExposure(energies, exposure, verbose);
         for (int k = 0; k < nenergies; k++)
            exposureCube[k][indx] = exposure[k];
         indx++;
      }
   }
   std::cerr << "!" << std::endl;
   writeFitsFile(filename, lon, lat, energies, exposureCube, 
                 roiCenter.ra(), roiCenter.dec());
}

void ExposureMap::writeFitsFile(const std::string &filename,
                                std::vector<double> &lon,
                                std::vector<double> &lat,
                                std::vector<double> &energies,
                                std::vector< std::valarray<double> > &dataCube,
                                double ra0, double dec0) {
// Use CCfits to produce the ExposureMap FITS file; this
// implementation is based on the writeImage() example in the
// CCfits-1.2 distribution (absent the superfluous use of the generic
// algorithms).

#ifdef HAVE_CCFITS
   int nlon = lon.size();
   double lonstep = lon[1] - lon[0];
   int nlat = lat.size();
   double latstep = lat[1] - lat[0];
   double emin = energies[0];
   double estep = log(energies[1]/energies[0]);
   int nenergies = energies.size();

   long naxis = 3;
   long naxes[3] = {nlon, nlat, nenergies};

   std::auto_ptr<FITS> pFits(0);

   try { // to overwrite a possibly existing file...
      pFits.reset( new FITS("!"+filename, DOUBLE_IMG, naxis, naxes) );
   } catch (FITS::CantCreate) {
      std::cerr << "ExposureMap::writeFitsFile: Can't create file "
                << filename << std::endl;
      return;
   }

   long nelements = naxes[0]*naxes[1]*naxes[2];

// add the keywords and values
   pFits->pHDU().addKey("CRVAL1", 0., 
                        "reference value for longitude coordinate");
   pFits->pHDU().addKey("CRVAL2", 0., 
                        "reference value for latitude coordinate");
   pFits->pHDU().addKey("CDELT1", lonstep, "step in longitude coordinate");
   pFits->pHDU().addKey("CDELT2", latstep, "step in latitude coordinate");
   pFits->pHDU().addKey("CRPIX1", static_cast<float>(nlon/2), 
                      "reference pixel for longitude coordinate");
   pFits->pHDU().addKey("CRPIX2", static_cast<float>(nlat/2), 
                      "reference pixel for latitude coordinate");
   pFits->pHDU().addKey("CTYPE1", "degrees", "units for longitude");
   pFits->pHDU().addKey("CTYPE2", "degrees", "units for latitude");
   pFits->pHDU().addKey("CRVAL3", log(emin), 
                      "reference value for log_energy coordinate");
   pFits->pHDU().addKey("CDELT3", estep, "step in log_energy coordinate");
   pFits->pHDU().addKey("CRPIX3", 1., 
                      "reference pixel for log_energy coordinate");
   pFits->pHDU().addKey("CTYPE3", "log_MeV", "units for log_energy");
   pFits->pHDU().addKey("LONPOLE", ra0, "RA of the ROI center");
   pFits->pHDU().addKey("LATPOLE", dec0, "DEC of the ROI center");

// flatten the exposureCube into a single valarray
   std::valarray<double> exposure(nelements);
   int indx = 0;
   for (int k = 0; k < nenergies; k++) {
      int cindx = 0;
      for (int i = 0; i < nlon; i++) {
         for (int j = 0; j < nlat; j++) {
            exposure[indx] = dataCube[k][cindx];
            indx++;
            cindx++;
         }
      }
   }

// write the primary HDU
   long fpixel(1);
   pFits->pHDU().write(fpixel, nelements, exposure);
#endif // HAVE_CCFITS
}

} // namespace Likelihood
