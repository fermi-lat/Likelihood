/** 
 * @file FitsImage.cxx
 * @brief Implementation of FitsImage member functions
 * @authors J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FitsImage.cxx,v 1.18 2004/09/28 14:45:37 jchiang Exp $
 *
 */

#include <iostream>
#include <memory>
#include <sstream>

#include "Likelihood/FitsImage.h"

namespace Likelihood {

#include "fitsio.h"

FitsImage::FitsImage(const std::string &fitsfile) 
   : st_facilities::FitsImage(fitsfile), m_eqRot(0) {
   if (haveRefCoord()) {
      m_eqRot = new EquinoxRotation(m_roiRa, m_roiDec);
   }
}

FitsImage::FitsImage(const FitsImage &rhs) : st_facilities::FitsImage(rhs) {
   m_roiRa = rhs.m_roiRa;
   m_roiDec = rhs.m_roiDec;
   if (rhs.m_eqRot) m_eqRot = rhs.m_eqRot->clone();
}

FitsImage& FitsImage::operator=(const FitsImage &rhs) {
   if (this != &rhs) {
      m_filename = rhs.m_filename;
      m_axes = rhs.m_axes;
      m_axisVectors = rhs.m_axisVectors;
      m_image = rhs.m_image;
      m_roiRa = rhs.m_roiRa;
      m_roiDec = rhs.m_roiDec;
      if (rhs.m_eqRot) {
         delete m_eqRot;
         m_eqRot = rhs.m_eqRot->clone();
      }
   }
   return *this;
}

void FitsImage::getCelestialArrays(std::vector<double> &lonArray,
                                   std::vector<double> &latArray) {
   unsigned int npixels = m_axes[0].size*m_axes[1].size;
   lonArray.resize(npixels);
   latArray.resize(npixels);
   for (int j = 0; j < m_axes[1].size; j++) {
      for (int i = 0; i < m_axes[0].size; i++) {
         int indx = i + j*m_axes[0].size;
         if (m_eqRot) {
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

bool FitsImage::haveRefCoord() {
   bool have_ref_coord;

   fitsfile * fptr = 0;
   char * file = const_cast<char *>(m_filename.c_str());
   int status = 0;

   fits_open_file(&fptr, file, READONLY, &status);
   fitsReportError(stderr, status);

   double value;
   char comment[80];
   fits_read_key_dbl(fptr, "ROI_RA", &value, comment, &status);
   if (status == 0) {
      m_roiRa = value;
      have_ref_coord = true;
   } else if (status == KEY_NO_EXIST) {
      have_ref_coord = false;
   } else {
      fitsReportError(stderr, status);
   }
   status = 0;

   fits_read_key_dbl(fptr, "ROI_DEC", &value, comment, &status);
   if (status == 0) {
      m_roiDec = value;
      have_ref_coord = true;
   } else if (status == KEY_NO_EXIST) {
      have_ref_coord = false;
   } else {
      fitsReportError(stderr, status);
   }
   status = 0;

   fits_close_file(fptr, &status);
   fitsReportError(stderr, status);

   return have_ref_coord;
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
      for (int j = 0; j < 3; j++) {
         outVec[i] += rotMatrix[i][j]*inVec[j];
      }
   }
   outDir = astro::SkyDir(Hep3Vector(outVec[0], outVec[1], outVec[2]));
}   

void FitsImage::fitsReportError(FILE *stream, int status) const {
   fits_report_error(stream, status);
   if (status != 0) {
      throw std::string("Likelihood::FitsImage: cfitsio error.");
   }
}

} // namespace Likelihood
