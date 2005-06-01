/** 
 * @file FitsImage.cxx
 * @brief Implementation of FitsImage member functions
 * @authors J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FitsImage.cxx,v 1.24 2005/05/21 23:39:02 jchiang Exp $
 *
 */

#include <iostream>
#include <memory>
#include <stdexcept>
#include <sstream>

#include "Likelihood/EquinoxRotation.h"
#include "Likelihood/FitsImage.h"

namespace Likelihood {

std::string FitsImage::s_routineName("");

FitsImage::FitsImage(const std::string &fitsfile) 
   : st_facilities::FitsImage(fitsfile), m_eqRot(0) {
   if (haveRefCoord()) {
      m_eqRot = new EquinoxRotation(m_roiRa, m_roiDec);
   }
   std::vector<std::string> axisNames;
   getAxisNames(axisNames);
   if (axisNames[0].find("RA") != std::string::npos) {
      m_coordSys = "Equatorial";
   } else if (axisNames[0].find("GLON") != std::string::npos) {
      m_coordSys = "Galactic";
   } else {
      std::ostringstream message;
      message << "Likelihood::FitsImage:\n"
              << "Unrecognized coordinate system in " << fitsfile << ".\n"
              << "Axis names: ";
      for (unsigned int i = 0; i < axisNames.size(); i++) {
         message << axisNames.at(i) << "  ";
      }
      throw std::runtime_error(message.str());
   }
}

FitsImage::FitsImage(const FitsImage &rhs) : st_facilities::FitsImage(rhs) {
   m_roiRa = rhs.m_roiRa;
   m_roiDec = rhs.m_roiDec;
   m_coordSys = rhs.m_coordSys;
   if (rhs.m_eqRot) {
      delete m_eqRot;
      m_eqRot = rhs.m_eqRot->clone();
   }
}

FitsImage::~FitsImage() {
   delete m_eqRot;
}

FitsImage& FitsImage::operator=(const FitsImage &rhs) {
   if (this != &rhs) {
      m_filename = rhs.m_filename;
      m_axes = rhs.m_axes;
      m_axisVectors = rhs.m_axisVectors;
      m_image = rhs.m_image;
      m_roiRa = rhs.m_roiRa;
      m_roiDec = rhs.m_roiDec;
      m_coordSys = rhs.m_coordSys;
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

void FitsImage::getPixelBounds(unsigned int naxis,
                               std::vector<double> & pixelBounds) const {
   getAxisVector(naxis, pixelBounds);
   double pixelSize = std::fabs(pixelBounds.at(1) - pixelBounds.at(0));
   if (pixelBounds.back() > pixelBounds.front()) {
      pixelBounds.push_back(pixelBounds.back() + pixelSize);
      for (unsigned int i = 0; i < pixelBounds.size(); i++) {
         pixelBounds.at(i) -= pixelSize;
      }
   } else {
      pixelBounds.insert(pixelBounds.begin(),
                         pixelBounds.front() + pixelSize);
   }
}

bool FitsImage::haveRefCoord() {
   s_routineName = "haveRefCoord";

   bool have_ref_coord;

   fitsfile * fptr = 0;
   char * file = const_cast<char *>(m_filename.c_str());
   int status = 0;

   fits_open_file(&fptr, file, READONLY, &status);
   fitsReportError(status);

   double value;
   char comment[80];
   fits_read_key_dbl(fptr, "ROI_RA", &value, comment, &status);
   if (status == 0) {
      m_roiRa = value;
      have_ref_coord = true;
   } else if (status == KEY_NO_EXIST) {
      have_ref_coord = false;
   } else {
      fitsReportError(status);
   }
   status = 0;

   fits_read_key_dbl(fptr, "ROI_DEC", &value, comment, &status);
   if (status == 0) {
      m_roiDec = value;
      have_ref_coord = true;
   } else if (status == KEY_NO_EXIST) {
      have_ref_coord = false;
   } else {
      fitsReportError(status);
   }
   status = 0;

   fits_close_file(fptr, &status);
   fitsReportError(status);

   return have_ref_coord;
}

int FitsImage::findHdu(const std::string & fitsFile, 
                       const std::string & extension) {
   s_routineName = "findHdu";
   
   int status(0);
   fitsfile * fptr = 0;

   fits_open_file(&fptr, fitsFile.c_str(), READONLY, &status);
   fitsReportError(status);
   
   int nhdus;
   fits_get_num_hdus(fptr, &nhdus, &status);
   fitsReportError(status);

   int hdutype(0);
   char extname[20];
   char comment[72];
   for (int hdu = 1; hdu < nhdus+1; hdu++) {
      fits_movabs_hdu(fptr, hdu, &hdutype, &status);
      fitsReportError(status);
      
      fits_read_key_str(fptr, "EXTNAME", extname, comment, &status);
      if (status == 202) {
         status = 0;
         continue;
      } else {
         fitsReportError(status);
      }
      
      if (extension == extname) {
         fits_close_file(fptr, &status);
         fitsReportError(status);
         return hdu;
      }
   }
   fits_close_file(fptr, &status);
   fitsReportError(status);

   std::ostringstream message;
   message << "FitsImage::findHdu: HDU number not found for file "
           << fitsFile << " and extension " << extension;
   throw std::runtime_error(message.str());
   return -1;
}

void FitsImage::readColumn(fitsfile * fptr, const std::string & colname,
                           std::vector<double> & coldata) {
   s_routineName = "readColumn";

   int status(0);
   int colnum(0);
   fits_get_colnum(fptr, CASEINSEN, const_cast<char *>(colname.c_str()),
                   &colnum, &status);
   fitsReportError(status);

   long nrows(0);
   fits_get_num_rows(fptr, &nrows, &status);
   fitsReportError(status);

   int anynul(0), nulval(0);
   coldata.resize(nrows);
   fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, &nulval, &coldata[0],
                 &anynul, &status);
   fitsReportError(status);
}

void FitsImage::fitsReportError(int status, std::string routine) {
   if (status == 0) {
      return;
   }
   if (routine == "") {
      routine = "FitsImage::" + s_routineName;
   }
   fits_report_error(stderr, status);
   std::ostringstream message;
   message << routine << ": CFITSIO error " << status;
   throw std::runtime_error(message.str());
}

} // namespace Likelihood
