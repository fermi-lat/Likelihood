/** @file FitsImage.cxx
 * @brief Implementation of FitsImage member functions
 * @authors J. Chiang
 *
 * $Header$
 *
 */

#include "fitsio.h"
#include <iostream>
#include "Likelihood/FitsImage.h"

namespace Likelihood {

FitsImage::FitsImage(std::string &fitsfile) {
   m_filename = fitsfile;
   read_fits_image(fitsfile, m_axes, m_image);
   for (unsigned int i = 0; i < m_axes.size(); i++) {
      std::vector<double> axisVector;
      m_axes[i].computeAxisVector(axisVector);
      m_axisVectors.push_back(axisVector);
   }
   if (m_haveRefCoord) m_eqRot = new EquinoxRotation(m_lonpole, m_latpole);
}

void FitsImage::fetchAxisDims(std::vector<int> &axisDims) {
   if (!axisDims.empty()) axisDims.clear();

   for (unsigned int i = 0; i < m_axes.size(); i++)
      axisDims.push_back(m_axes[i].size);
}

void FitsImage::fetchAxisNames(std::vector<std::string> &axisNames) {   
   if (!axisNames.empty()) axisNames.clear();

   for (unsigned int i = 0; i < m_axes.size(); i++)
      axisNames.push_back(m_axes[i].axisType);
}

void FitsImage::fetchAxisVector(unsigned int naxis,
                                std::vector<double> &axisVector) {
   if (naxis >= m_axes.size()) {
      std::cerr << "FitsImage::fetchAxisVector: Invalid axis number "
                << naxis << std::endl;
      assert(false);
   }
   axisVector = m_axisVectors[naxis];
}

void FitsImage::fetchCelestialArrays(std::valarray<double> &lonArray,
                                     std::valarray<double> &latArray) {
   unsigned int npixels = m_axes[0].size*m_axes[1].size;
   lonArray.resize(npixels);
   latArray.resize(npixels);
   for (int j = 0; j < m_axes[1].size; j++) {
      for (int i = 0; i < m_axes[0].size; i++) {
         int indx = i + j*m_axes[0].size;
         if (m_haveRefCoord) {
            astro::SkyDir inVec(m_axisVectors[0][i], m_axisVectors[1][j]);
            astro::SkyDir outVec;
            m_eqRot->do_rotation(inVec, outVec);
            lonArray[indx] = outVec.ra();
            latArray[indx] = outVec.dec();
         } else {
            lonArray[indx] = m_axisVectors[0][i];
            latArray[indx] = m_axisVectors[1][j];
         }
      }
   }
}
         
void FitsImage::fetchSolidAngles(std::valarray<double> &solidAngles){
// This solid angle calculation *assumes* that m_axes[0] is a
// longitudinal coordinate and that m_axes[1] is a latitudinal one.
// Furthermore, the axis units are assumed to be *degrees*, while the
// solid angles are returned as steradians.

   solidAngles.resize(m_axes[0].size*m_axes[1].size);
   for (int i = 0; i < m_axes[0].size; i++) {
      for (int j = 0; j < m_axes[1].size; j++) {
         int indx = i + j*m_axes[0].size;
         double thetamin = (m_axisVectors[1][j] - m_axes[1].step/2.)*M_PI/180;
         double thetamax = (m_axisVectors[1][j] + m_axes[1].step/2.)*M_PI/180;
         solidAngles[indx] = m_axes[0].step*M_PI/180
            *(sin(thetamax) - sin(thetamin));
      }
   }
}

void FitsImage::fetchImageData(std::valarray<double> &imageData) {
   imageData.resize(m_image.size());
   imageData = m_image;
}

void FitsImage::AxisParams::computeAxisVector(std::vector<double> &axisVector) {
   if (!axisVector.empty()) axisVector.clear();
   axisVector.reserve(size);
   for (int i = 0; i < size; i++) {
      double value = step*(i - refPixel + 1) + refVal;
      if (logScale) value = exp(value);
      axisVector.push_back(value);
   }
}

void FitsImage::read_fits_image(std::string &filename, 
                                std::vector<AxisParams> &axes,
                                std::valarray<double> &image) {
   
   fitsfile * fptr = 0;
   char *file = const_cast<char *>(filename.c_str());
   int status = 0;

   fits_open_file(&fptr, file, READONLY, &status);
   fits_report_error(stderr, status);

// Get dimensions of the data cube
   long naxes;
   char comment[80];
   fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status);
   fits_report_error(stderr, status);

// Assume at least 1 image plane, but at most 3 dimensions...
   if (naxes != 2 && naxes != 3) {
      std::cerr << "FitsImage::read_fits_image: "
                << "FITS file " << filename 
                << " does not have the expected number of dimensions:"
                << " naxes = " << naxes
                << std::endl;
      assert(false);
   }

// prepare the axes vector 
   axes.clear();
   axes.resize(naxes);

// keyword names
   char *naxis[] = {"NAXIS1", "NAXIS2", "NAXIS3"};
   char *crval[] = {"CRVAL1", "CRVAL2", "CRVAL3"};
   char *cdelt[] = {"CDELT1", "CDELT2", "CDELT3"};
   char *crpix[] = {"CRPIX1", "CRPIX2", "CRPIX3"};
   char *ctype[] = {"CTYPE1", "CTYPE2", "CTYPE3"};

   long ivalue;
   double value;
   char svalue[40];

   long npixels = 1;
   for (int i = 0; i < naxes; i++) {
// axis size
      fits_read_key_lng(fptr, naxis[i], &ivalue, comment, &status);
      fits_report_error(stderr, status);
      axes[i].size = ivalue;

// Compute the number of pixels in the image.
      npixels *= ivalue;
   }

// account for degenerate case of NAXIS3 = 1

   if (naxes == 3 && axes[2].size == 1) {
      naxes = 2;
      axes.resize(naxes);
   }

   for (int i = 0; i < naxes; i++) {
// reference values
      fits_read_key_dbl(fptr, crval[i], &value, comment, &status);
      fits_report_error(stderr, status);
      axes[i].refVal = value;

// step sizes
      fits_read_key_dbl(fptr, cdelt[i], &value, comment, &status);
      fits_report_error(stderr, status);
      axes[i].step = value;

// reference pixels
      fits_read_key_lng(fptr, crpix[i], &ivalue, comment, &status);
      fits_report_error(stderr, status);
      axes[i].refPixel = ivalue;

// axis types and commentary
      fits_read_key_str(fptr, ctype[i], svalue, comment, &status);
      fits_report_error(stderr, status);
      axes[i].axisType = svalue;
      axes[i].comment = comment;
      
// Check for logarithmic scaling.
      int offset = axes[i].axisType.substr(0).find_first_of("log_");
      if (offset == 0) {
         axes[i].logScale = true;
      } else {
         axes[i].logScale = false;
      }
   } // naxes

// Check for LONPOLE and LATPOLE keywords.

   fits_read_key_dbl(fptr, "LONPOLE", &value, comment, &status);
   if (status == 0) {
      m_lonpole = value;
      m_haveRefCoord = true;
   } else if (status == KEY_NO_EXIST) {
      m_haveRefCoord = false;
   } else {
      std::cerr << KEY_NO_EXIST << " " << status << std::endl;
      fits_report_error(stderr, status);
   }
   status = 0;

   fits_read_key_dbl(fptr, "LATPOLE", &value, comment, &status);
   if (status == 0) {
      m_latpole = value;
      m_haveRefCoord = true;
   } else if (status == KEY_NO_EXIST) {
      m_haveRefCoord = false;
   } else {
      std::cerr << KEY_NO_EXIST << " " << status << std::endl;
      fits_report_error(stderr, status);
   }
   status = 0;

// Read in the image pixels.
   long group = 0;
   long fpixel = 1;
   double nullval = 0.;
   int anynull;
   double *tmpImage;
   tmpImage = new double[npixels];
   fits_read_img_dbl(fptr, group, fpixel, npixels, nullval, 
                     tmpImage, &anynull, &status);
   fits_report_error(stderr, status);

   image.resize(npixels);

   for (int i = 0; i < npixels; i++)
      image[i] = tmpImage[i];

   delete [] tmpImage;

   fits_close_file(fptr, &status);
   fits_report_error(stderr, status);
}

FitsImage::EquinoxRotation::EquinoxRotation(double alpha0, double delta0) {

   rotMatrix.clear();

// convert to radians
   alpha0 *= M_PI/180;
   delta0 *= M_PI/180;

// build the rotation matrix, using Fortran-like indexing, i.e.,
// rotMatrix[row][column].  Note that this is the *transpose* of the
// matrix given in LikeMemo 3.
   
   double ca = cos(alpha0);
   double sa = sin(alpha0);
   double cd = cos(delta0);
   double sd = sin(delta0);

   std::vector<double> row(3);

   row[0] = cd*ca;
   row[1] = -sa;
   row[2] = -sd*ca;
   rotMatrix.push_back(row);

   row[0] = cd*sa;
   row[1] = ca;
   row[2] = -sd*sa;
   rotMatrix.push_back(row);

   row[0] = sd;
   row[1] = 0;
   row[2] = cd;
   rotMatrix.push_back(row);
}

void FitsImage::EquinoxRotation::do_rotation(const astro::SkyDir &inDir,
                                             astro::SkyDir &outDir) {
// I'm sure there is a Hep3Vector way to do this, nonetheless....

   std::vector<double> inVec(3), outVec(3);

// Need to ensure that inDir is created using 
// (ra, dec) = (longitude, latitude)
   double alpha = inDir.ra()*M_PI/180;
   double delta = inDir.dec()*M_PI/180;
   inVec[0] = cos(delta)*cos(alpha);
   inVec[1] = cos(delta)*sin(alpha);
   inVec[2] = sin(delta);

// Apply the rotation

   for (int i = 0; i < 3; i++) {
      outVec[i] = 0;
      for (int j = 0; j < 3; j++)
         outVec[i] += rotMatrix[i][j]*inVec[j];
   }
   outDir = astro::SkyDir(Hep3Vector(outVec[0], outVec[1], outVec[2]));
}   

} // namespace Likelihood
