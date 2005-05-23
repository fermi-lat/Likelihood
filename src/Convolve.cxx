/**
 * @file Convolve.cxx
 * @brief Static methods to perform 1D and 2D convolutions using the
 *        FFTW library
 * @author J. Chiang
 *
 * $Header$
 */

#include <iostream>
#include <stdexcept>
#include <vector>

#include <fftw3.h>

#include "Likelihood/Convolve.h"

namespace {
   fftw_complex * 
   complexVector(const std::vector<double> & input, bool fill=true) {
      unsigned int npts = input.size();
      fftw_complex * output;
      output = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*npts);
      if (fill) {
         for (unsigned int i = 0; i < npts; i++) {
            output[i][0] = input.at(i);
            output[i][1] = 0;
         }
      }
      return output;
   }

   fftw_complex * 
   complexArray(const std::vector< std::vector<double> > & input,
                bool fill=true, unsigned int nx=0, unsigned int ny=0) {
      unsigned int dx(0), dy(0);
      if (ny == 0 || ny < input.size()) {
         ny = input.size();
      } else {
         dy = (ny - input.size())/2;
      }
      if (nx == 0 || nx < input.at(0).size()) {
         nx = input.at(0).size();
      } else {
         dx = (nx - input.at(0).size())/2;
      }
      fftw_complex * output;
      output = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nx*ny);
      if (fill) {
         for (unsigned int i = 0; i < nx*ny; i++) {
            output[i][0] = 0;
            output[i][1] = 0;
         }
         for (unsigned int i = 0; i < input.size(); i++) {
            for (unsigned int j = 0; j < input.at(i).size(); j++) {
               unsigned int indx = (i + dx)*ny + j + dy;
               if (indx < nx*ny) {
                  output[indx][0] = input.at(i).at(j);
               }
            }
         }
      }
      return output;
   }
} // unnamed namespace

namespace Likelihood {

std::vector< std::vector<double> > 
Convolve::convolve2d(const std::vector< std::vector<double> > & signal,
                     const std::vector< std::vector<double> > & psf) {
   if (psf.size() > signal.size()) {
      throw std::runtime_error("Convolve::convolve2d: Psf size must "
                               "be smaller than the signal size.");
   }

   int nx = signal.at(0).size();
   int ny = signal.size();
   unsigned int npts = nx*ny;

   fftw_complex * in = ::complexArray(signal);
   fftw_complex * out = ::complexArray(signal, false);
   fftw_plan splan = fftw_plan_dft_2d(nx, ny, in, out,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(splan);

   fftw_complex * ipsf = ::complexArray(psf, true, nx, ny);
   fftw_complex * opsf = ::complexArray(psf, false, nx, ny);
   fftw_plan pplan = fftw_plan_dft_2d(nx, ny, ipsf, opsf,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(pplan);

   fftw_complex * iconv = ::complexArray(signal, false);
   fftw_complex * oconv = ::complexArray(signal, false);
   for (unsigned int i = 0; i < npts; i++) {
      iconv[i][0] = out[i][0]*opsf[i][0] - out[i][1]*opsf[i][1];
      iconv[i][1] = out[i][0]*opsf[i][1] + out[i][1]*opsf[i][0];
   }
   fftw_plan cplan = fftw_plan_dft_2d(nx, ny, iconv, oconv,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(cplan);

   unsigned int indx;
   std::vector< std::vector<double> > output;
   output.reserve(signal.size());
   for (unsigned int i = signal.size()/2 - 1; i < signal.size(); i++) {
      std::vector<double> local;
      for (unsigned int j = signal.at(i).size()/2 - 1;
           j < signal.at(i).size(); j++) {
         indx = i*signal.at(i).size() + j;
         local.push_back(oconv[indx][0]/npts);
      }
      for (unsigned int j = 0; j < signal.at(i).size()/2 - 1; j++) {
         indx = i*signal.at(i).size() + j;
         local.push_back(oconv[indx][0]/npts);
      }
      output.push_back(local);
   }
   for (unsigned int i = 0; i < signal.size()/2 - 1; i++) {
      std::vector<double> local;
      for (unsigned int j = signal.at(i).size()/2 - 1;
           j < signal.at(i).size(); j++) {
         indx = i*signal.at(i).size() + j;
         local.push_back(oconv[indx][0]/npts);
      }
      for (unsigned int j = 0; j < signal.at(i).size()/2 - 1; j++) {
         indx = i*signal.at(i).size() + j;
         local.push_back(oconv[indx][0]/npts);
      }
      output.push_back(local);
   }

   fftw_destroy_plan(splan);
   fftw_destroy_plan(pplan);
   fftw_destroy_plan(cplan);

   fftw_free(in);
   fftw_free(out);
   fftw_free(ipsf);
   fftw_free(opsf);
   fftw_free(iconv);
   fftw_free(oconv);

   return output;
}

std::vector<double> 
Convolve::convolve(const std::vector<double> & signal,
                   const std::vector<double> & psf) {
   if (psf.size() != signal.size()) {
      throw std::runtime_error("convolve: psf size must be the same "
                               "as the signal size.");
   }

   fftw_complex * in = ::complexVector(signal);
   fftw_complex * out = ::complexVector(signal, false);
   fftw_complex * ipsf = ::complexVector(psf);
   fftw_complex * opsf = ::complexVector(psf, false);

   fftw_plan splan = fftw_plan_dft_1d(signal.size(), in, out,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_plan pplan = fftw_plan_dft_1d(signal.size(), ipsf, opsf,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(splan);
   fftw_execute(pplan);

   fftw_complex * iconv = ::complexVector(signal, false);
   fftw_complex * oconv = ::complexVector(signal, false);
   for (unsigned int i = 0; i < signal.size(); i++) {
      iconv[i][0] = out[i][0]*opsf[i][0] - out[i][1]*opsf[i][1];
      iconv[i][1] = out[i][0]*opsf[i][1] + out[i][1]*opsf[i][0];
   }
   fftw_plan cplan = fftw_plan_dft_1d(signal.size(), iconv, oconv,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(cplan);

   std::vector<double> output;
   output.reserve(signal.size());
   for (unsigned int i = signal.size()/2-1; i < signal.size(); i++) {
      output.push_back(oconv[i][0]/signal.size());
   }
   for (unsigned int i = 0; i < signal.size()/2-1; i++) {
      output.push_back(oconv[i][0]/signal.size());
   }

   fftw_destroy_plan(splan);
   fftw_destroy_plan(pplan);
   fftw_destroy_plan(cplan);

   fftw_free(in);
   fftw_free(out);
   fftw_free(ipsf);
   fftw_free(opsf);
   fftw_free(iconv);
   fftw_free(oconv);

   return output;
}

} // namespace Likelihood
