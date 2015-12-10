/**
 * @file Convolve.cxx
 * @brief Static methods to perform 1D and 2D convolutions using the
 *        FFTW library
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/Convolve.cxx,v 1.3 2015/03/03 06:00:00 echarles Exp $
 */

#include <iostream>
#include <stdexcept>
#include <vector>

#include "fftw/fftw3.h"

#include "Likelihood/Convolve.h"

namespace {
   fftw_complex * 
   complexVector(const std::vector<double> & input, bool fill=true) {
      size_t npts = input.size();
      fftw_complex * output;
      output = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*npts);
      if (fill) {
         for (size_t i = 0; i < npts; i++) {
            output[i][0] = input.at(i);
            output[i][1] = 0;
         }
      }
      return output;
   }

   fftw_complex * 
   complexArray(const std::vector< std::vector<double> > & input,
                bool fill=true, size_t nx=0, size_t ny=0) {
      size_t dx(0), dy(0);
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
         for (size_t i = 0; i < nx*ny; i++) {
            output[i][0] = 0;
            output[i][1] = 0;
         }
         for (size_t j = 0; j < input.size(); j++) {
            for (size_t i = 0; i < input.at(j).size(); i++) {
               size_t indx = (j + dy)*nx + i + dx;
               if (indx < nx*ny) {
                  output[indx][0] = input.at(j).at(i);
               }
            }
         }
      }
      return output;
   }
   std::vector< std::vector<double> > 
   subArray(std::vector< std::vector<double> > & input,
            const std::vector< std::vector<double> > & signal) {
      size_t dx = (input.size() - signal.size())/2;
      size_t dy = (input.at(0).size() - signal.at(0).size())/2;
      std::vector< std::vector<double> > output(signal.size());
      for (size_t i = 0; i < signal.size(); i++) {
         output.at(i).resize(signal.at(i).size());
         for (size_t j = 0; j < signal.at(i).size(); j++) {
            output.at(i).at(j) = input.at(i + dx).at(j + dy);
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

   size_t nx = signal.at(0).size();
   size_t ny = signal.size();
// pad out the signal array so that it is n x n
   if (nx < ny) {
      nx = ny;
   } else if (ny < nx) {
      ny = nx;
   }
   nx*=2;
   ny*=2;
   size_t npts = nx*ny;

   fftw_complex * in = ::complexArray(signal, true, nx, ny);
   fftw_complex * out = ::complexArray(signal, false, nx, ny);
   fftw_plan splan = fftw_plan_dft_2d(nx, ny, in, out,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(splan);

   fftw_complex * ipsf = ::complexArray(psf, true, nx, ny);
   fftw_complex * opsf = ::complexArray(psf, false, nx, ny);
   fftw_plan pplan = fftw_plan_dft_2d(nx, ny, ipsf, opsf,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(pplan);

   fftw_complex * iconv = ::complexArray(signal, false, nx, ny);
   fftw_complex * oconv = ::complexArray(signal, false, nx, ny);
   for (size_t i = 0; i < npts; i++) {
      iconv[i][0] = out[i][0]*opsf[i][0] - out[i][1]*opsf[i][1];
      iconv[i][1] = out[i][0]*opsf[i][1] + out[i][1]*opsf[i][0];
   }
   fftw_plan cplan = fftw_plan_dft_2d(nx, ny, iconv, oconv,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(cplan);

   size_t indx;
   std::vector< std::vector<double> > output;
   output.reserve(ny);
   
   for (size_t i = ny/2 - 1; i < ny; i++) {
      std::vector<double> local;
      for (size_t j = nx/2 - 1; j < nx; j++) {
         indx = i*nx + j;
         local.push_back(oconv[indx][0]/npts);
      }
      for (size_t j = 0; j < nx/2 - 1; j++) {
         indx = i*nx + j;
         local.push_back(oconv[indx][0]/npts);
      }
      output.push_back(local);
   }
   for (size_t i = 0; i < ny/2 - 1; i++) {
      std::vector<double> local;
      for (size_t j = nx/2 - 1; j < nx; j++) {
         indx = i*nx + j;
         local.push_back(oconv[indx][0]/npts);
      }
      for (size_t j = 0; j < nx/2 - 1; j++) {
         indx = i*nx + j;
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

   return ::subArray(output, signal);
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
   for (size_t i = 0; i < signal.size(); i++) {
      iconv[i][0] = out[i][0]*opsf[i][0] - out[i][1]*opsf[i][1];
      iconv[i][1] = out[i][0]*opsf[i][1] + out[i][1]*opsf[i][0];
   }
   fftw_plan cplan = fftw_plan_dft_1d(signal.size(), iconv, oconv,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(cplan);

   std::vector<double> output;
   output.reserve(signal.size());
   for (size_t i = signal.size()/2-1; i < signal.size(); i++) {
      output.push_back(oconv[i][0]/signal.size());
   }
   for (size_t i = 0; i < signal.size()/2-1; i++) {
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

std::vector<float> 
Convolve::convolve(const std::vector<float> & signal,
                   const std::vector<float> & psf) {
   std::vector<double> dsignal(signal.size(), 0);
   std::vector<double> dpsf(psf.size(), 0);
   for (size_t i(0); i < signal.size(); i++) {
      dsignal[i] = signal[i];
   }
   for (size_t i(0); i < psf.size(); i++) {
      dpsf[i] = psf[i];
   }
   std::vector<double> foo(convolve(dsignal, dpsf));
   std::vector<float> my_result(foo.size(), 0);
   for (size_t i(0); i < foo.size(); i++) {
      my_result[i] = foo[i];
   }
   return my_result;
}

std::vector< std::vector<float> > 
Convolve::convolve2d(const std::vector< std::vector<float> > & signal,
                     const std::vector< std::vector<float> > & psf) {
   std::vector< std::vector<double> > dsignal(signal.size(),
                                              std::vector<double>(signal[0].size(), 0));
   std::vector< std::vector<double> > dpsf(psf.size(),
                                           std::vector<double>(psf[0].size(), 0));
   for (size_t i(0); i < signal.size(); i++) {
      for (size_t j(0); j < signal[i].size(); j++) {
         dsignal[i][j] = signal[i][j];
      }
   }
   for (size_t i(0); i < psf.size(); i++) {
      for (size_t j(0); j < psf[i].size(); j++) {
         dpsf[i][j] = psf[i][j];
      }
   }
   std::vector< std::vector<double> > foo(convolve2d(dsignal, dpsf));
   std::vector< std::vector<float> > my_result(foo.size(),
                                               std::vector<float>(foo[0].size(), 0));
   for (size_t i(0); i < foo.size(); i++) {
      for (size_t j(0); j < foo[i].size(); j++) {
         my_result[i][j] = foo[i][j];
      }
   }
   return my_result;
}


} // namespace Likelihood
