/**
 * @file Mcmc.cxx
 * @brief Implementation for generating a Markov Chain Monte Carlo of
 * a Statistic object using the Variable-at-a-time Metropolis-Hastings
 * update method.
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Mcmc.cxx,v 1.5 2003/06/04 16:24:40 jchiang Exp $
 */

#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "fitsio.h"
#include "Likelihood/Mcmc.h"
#include "LikelihoodException.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

namespace Likelihood {

Mcmc::Mcmc(Statistic &stat) {
   m_stat = &stat;
   estimateTransWidths();
}

void Mcmc::generateSamples(std::vector< std::vector<double> > &samples, 
                           unsigned long nsamp) {
// Get initial values of Parameters... 
   std::vector<double> paramValues;
   m_stat->getFreeParamValues(paramValues);

// and the Parameter objects themselves (for the bounds information)
   std::vector<Parameter> params;
   m_stat->getFreeParams(params);

   samples.clear();

   while (samples.size() < nsamp) {
// Loop over parameters, treating each update step as a trial
      for (unsigned int i = 0; i < paramValues.size(); i++) {
         std::vector<double> newParamValues = paramValues;
         double transProbRatio;
         newParamValues[i] = drawValue(params[i], m_transitionWidths[i],
                                       transProbRatio);
// Hastings ratio
         double alpha = transProbRatio*exp(m_stat->value(newParamValues)
                                           - m_stat->value(paramValues));
// Metropolis rejection criterion (Use CLHEP call here?)
//          double drand = static_cast<double>(rand())
//             /static_cast<double>(RAND_MAX);
         double drand = RandFlat::shoot();
         if (drand < alpha) {
// Accept the new point in Parameter space
            paramValues[i] = newParamValues[i];
            params[i].setValue(paramValues[i]);
         }
// We always append the current point after the update step
         samples.push_back(paramValues);
      }
   }
}

void Mcmc::writeSamples(std::string filename, 
                        std::vector< std::vector<double> > &samples) 
   throw(LikelihoodException) {

   fitsfile *fptr;
   int status = 0;

   std::cout << filename << std::endl;
// create the file
   remove(filename.c_str());
   fits_create_file(&fptr, filename.c_str(), &status);
   if (status != 0) {
      fits_report_error(stderr, status);
      throw LikelihoodException("Mcmc::writeSamples: cfitsio errors.");
   }

// create the binary table, with the labels identifying the content of
// each column

   int tfields = static_cast<int>(samples[0].size());
//   static int tfields_max = 100;
   char *ttype[100];
   char *tform[100];
   char *tunit[100];

   std::cout << "creating binary table" << std::endl;
   for (int i = 0; i < tfields; i++) {
      std::ostringstream type;
      type << "param" << i;
      ttype[i] = strdup(std::string(type.str()).c_str());
      tform[i] = strdup("1E");
      tunit[i] = strdup("None");
      std::cout << ttype[i] << "  "
                << tform[i] << "  "
                << tunit[i] << std::endl;
   }
   char *extname = strdup("Mcmc data");
   int nrows = samples.size();
   fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform,
		   tunit, extname, &status);
   if (status != 0) {
      fits_report_error(stderr, status);
      throw LikelihoodException("Mcmc::writeSamples: cfitsio errors.");
   }

   int firstrow  = 1;  /* first row in table to write   */
   int firstelem = 1;  /* first element in row  (ignored in ASCII tables) */

// write each column
   for (int i = 0; i < tfields; i++) {
// repack the data into a vector of floats
      std::vector<float> my_data;
      for (int j = 0; j < nrows; j++) {
         my_data.push_back(samples[j][i]);
      }
// Since the data in vectors are stored sequentially, one can pass the
// pointer to the first object as a C array.
      fits_write_col(fptr, TFLOAT, i+1, firstrow, firstelem, nrows, 
                     &my_data[0], &status);
      if (status != 0) {
         fits_report_error(stderr, status);
         throw LikelihoodException("Mcmc::writeSamples: cfitsio errors.");
      }
   }
   fits_close_file(fptr, &status);
   if (status != 0) {
      fits_report_error(stderr, status);
      throw LikelihoodException("Mcmc::writeSamples: cfitsio errors.");
   }
}

void Mcmc::estimateTransWidths() {
                               
   m_transitionWidths.clear();

// Use an approximate Hessian to get estimates of each Parameters'
// error bars.
   std::vector<double> params;
   m_stat->getFreeParamValues(params);
   std::vector<double> derivs;
   m_stat->getFreeDerivs(derivs);
   double eps = 1e-5;

// Estimate the diagonal elements of the Hessian
   for (unsigned int i = 0; i < params.size(); i++) {
      std::vector<double> new_params = params;
      double delta = eps*params[i];
      new_params[i] += delta;
      m_stat->value(new_params);
      std::vector<double> new_derivs;
      m_stat->getFreeDerivs(new_derivs);
      double hessian = (new_derivs[i] - derivs[i])/delta;
      m_transitionWidths.push_back(sqrt(abs(1./hessian)));
   }
}

double Mcmc::drawValue(Parameter &param, double dx, double &transProbRatio) {
// Normalize the top-hat function (i.e., find the width) at the input
// Parameter location given the Parameter bounds.

   double xl = param.getBounds().first;
   double xu = param.getBounds().second;
   double x0 = param.getValue();
   double width = std::min(xu, x0 + dx) - std::max(xl, x0 - dx);

// Draw the trial value...
//   double drand = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
   double drand = RandFlat::shoot();
   double y = drand*width + std::max(xl, x0 - dx);

// and compute the ratio of the transition probability densities
   transProbRatio = width/(std::min(xu, y + dx) - std::max(xl, y - dx));

   return y;
}

} // namespace Likelihood
