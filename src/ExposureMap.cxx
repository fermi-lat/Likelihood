/**
 * @file ExposureMap.cxx
 * @brief Implementation for the Singleton ExposureMap class. This
 * class encapsulates exposure map information and makes it available
 * for use (primarily) by the DiffuseSource class.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/ExposureMap.cxx,v 1.19 2004/04/07 20:18:05 jchiang Exp $
 */
#include <utility>
#include <algorithm>

#include "fitsio.h"

#include "facilities/Util.h"

#include "Likelihood/SkyDirArg.h"
#include "Likelihood/ExposureMap.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/RoiCuts.h"

namespace {

void fitsReportError(FILE *stream, int status) {
   fits_report_error(stream, status);
   if (status != 0) {
      throw std::string("writeExposureFile: cfitsio error.");
   }
}

} // unnamed namespace

namespace Likelihood {

ExposureMap * ExposureMap::s_instance = 0;
bool ExposureMap::s_haveExposureMap = false;
std::vector<double> ExposureMap::s_energies;
std::valarray<double> ExposureMap::s_ra;
std::valarray<double> ExposureMap::s_dec;
std::vector< std::valarray<double> > ExposureMap::s_exposure;
FitsImage ExposureMap::s_mapData;

void ExposureMap::readExposureFile(std::string exposureFile) {

// Expand any environment variables in the exposure file name.
   facilities::Util::expandEnvVar(&exposureFile);

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
                                       optimizers::Function * spatialDist,
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

void ExposureMap::computeMap(std::string filename, double sr_radius,
                             int nlon, int nlat, int nenergies) {
// Expand any environment variables in the map filename.
   facilities::Util::expandEnvVar(&filename);

   RoiCuts *roi_cuts = RoiCuts::instance();
   astro::SkyDir roiCenter = roi_cuts->extractionRegion().center();
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
         bool verbose(false);
         if (PointSource::s_exposure == 0) {
            ptsrc.computeExposure(energies, exposure, verbose);
         } else {
            ptsrc.computeExposureWithHyperCube(energies, exposure, verbose);
         }
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
                                std::vector<double> &glon,
                                std::vector<double> &glat,
                                std::vector<double> &energies,
                                std::vector< std::valarray<double> > &exposure,
                                double ra0, double dec0) {

   fitsfile *fptr;
   int status = 0;
   
// Always overwrite an existing file.
   remove(filename.c_str());
   fits_create_file(&fptr, filename.c_str(), &status);
   fitsReportError(stderr, status);

   int bitpix = DOUBLE_IMG;
   long naxis = 3;
   long naxes[] = {glon.size(), glat.size(), energies.size()};
   fits_create_img(fptr, bitpix, naxis, naxes, &status);
   fitsReportError(stderr, status);

// Write the exposure map data.
   long group = 0;
   long dim1 = glon.size();
   long dim2 = glat.size();

// Repack exposure into a C array.
   double *exp_array = new double[glon.size()*glat.size()*energies.size()];
   int indx = 0;
   for (unsigned int k = 0; k < energies.size(); k++) {
      for (unsigned int i = 0; i < exposure[k].size(); i++) {
         exp_array[indx] = exposure[k][i];
         indx++;
      }
   }

   fits_write_3d_dbl(fptr, group, dim1, dim2, 
                     glon.size(), glat.size(), energies.size(),
                     exp_array, &status);
   delete[] exp_array;
   fitsReportError(stderr, status);

// Write some keywords.
   double l0 = glon[0];
   fits_update_key(fptr, TDOUBLE, "CRVAL1", &l0, 
                   "longitude of reference pixel", &status);
   fitsReportError(stderr, status);
   double b0 = glat[0];
   fits_update_key(fptr, TDOUBLE, "CRVAL2", &b0, 
                   "latitude of reference pixel", &status);
   fitsReportError(stderr, status);
   
   double lstep = glon[1] - glon[0];
   fits_update_key(fptr, TDOUBLE, "CDELT1", &lstep, 
                   "longitude step at ref. pixel", &status);
   fitsReportError(stderr, status);
   double bstep = glat[1] - glat[0];
   fits_update_key(fptr, TDOUBLE, "CDELT2", &bstep, 
                   "latitude step at ref. pixel", &status);
   fitsReportError(stderr, status);
   
   float crpix1 = 1.0;
   fits_update_key(fptr, TFLOAT, "CRPIX1", &crpix1, 
                   "reference pixel for longitude coordinate", &status);
   fitsReportError(stderr, status);
   float crpix2 = 1.0;
   fits_update_key(fptr, TFLOAT, "CRPIX2", &crpix2, 
                   "reference pixel for latitude coordinate", &status);
   fitsReportError(stderr, status);
   
   char *ctype1 = "GLON-CAR";
   fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, 
                   "right ascension", &status);
   fitsReportError(stderr, status);
   char *ctype2 = "GLAT-CAR";
   fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, 
                   "declination", &status);
   fitsReportError(stderr, status);

   double logEmin = log(energies[0]);
   fits_update_key(fptr, TDOUBLE, "CRVAL3", &logEmin,
                   "reference value for log_energy coordinate", &status);
   fitsReportError(stderr, status);

   int nee = energies.size();
   double estep = log(energies[nee-1]/energies[0])/(nee-1);
   fits_update_key(fptr, TDOUBLE, "CDELT3", &estep, 
                   "step in log_energy coordinate", &status);
   fitsReportError(stderr, status);

   float crpix3 = 1.;
   fits_update_key(fptr, TFLOAT, "CRPIX3", &crpix3,
                   "reference pixel for log_energy coordinate", &status);
   fitsReportError(stderr, status);

   char *ctype3 = "log_MeV";
   fits_update_key(fptr, TSTRING, "CTYPE3", ctype3,
                   "units for log_energy", &status);
   fitsReportError(stderr, status);

   fits_update_key(fptr, TDOUBLE, "ROI_RA", &ra0, "RA of ROI center", 
                   &status);
   fitsReportError(stderr, status);
   fits_update_key(fptr, TDOUBLE, "ROI_DEC", &dec0, "DEC of ROI center",
                   &status);
   fitsReportError(stderr, status);
   
// Write the energy array as a binary table.
   int nrows = energies.size();
   int tfields = 1;
   char *ttype[] = {"Energy"};
   char *tform[] = {"1D"};
   char *tunit[] = {"MeV"};
   char extname[] = "True Photon Energy Array";
   
   int firstrow  = 1;
   int firstelem = 1;
   
   fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                   tunit, extname, &status);
   fitsReportError(stderr, status);
   
   fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, 
                  &energies[0], &status);
   fitsReportError(stderr, status);
   
   fits_close_file(fptr, &status);
   fitsReportError(stderr, status);
   
   return;
}

} // namespace Likelihood
