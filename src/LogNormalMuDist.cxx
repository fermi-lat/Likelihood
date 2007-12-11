/**
 * @file LogNormalMuDist.cxx
 * @brief Implementation of function to compute a log-Normal
 * distribution of sampling points in cos(theta) for diffuse response
 * integral.
 * 
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <deque>
#include <iostream>

#include "LogNormalMuDist.h"

namespace {
   double erf(double x) {
      double z = std::fabs(x);
      double t = 1.0/(1. + 0.5*z);
      double ans = t*std::exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
                   t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
                   t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
                    
      if (x >= 0) {
         return 1. - ans;
      }
      return  ans - 1.;
   }
}

namespace Likelihood {

LogNormalMuDist * LogNormalMuDist::s_instance(0);

LogNormalMuDist * LogNormalMuDist::instance() {
   if (!s_instance) {
      s_instance = new LogNormalMuDist();
   }
   return s_instance;
}

LogNormalMuDist::
LogNormalMuDist(double muSlope, double muIntercept, double sigma,
                double emin, double emax, size_t numEnergies)
   : m_sigma(sigma), m_logEmin(std::log(emin)), m_logEmax(std::log(emax)),
     m_numEnergies(numEnergies),
     m_estep((m_logEmax - m_logEmin)/(numEnergies - 1)) {

   m_muPoints.clear();
   for (size_t k(0); k < numEnergies - 1; k++) {
      double logE(m_estep*(k + 0.5) + m_logEmin);
      double mu(muSlope*logE + muIntercept);
      std::vector<double> mus;
      createSample(mu, mus);
      m_muPoints.push_back(mus);
   }
}

const std::vector<double> & LogNormalMuDist::muPoints(double energy) const {
   return m_muPoints.at(energyIndex(energy));
}

void LogNormalMuDist::createSample(double mu, std::vector<double> & mus) {
   double xmin(1e-4);
   double xmax(2.*M_PI);
   size_t nx(200);

   std::vector<double> xx;
   std::vector<double> cdf;
   double xstep(std::log(xmax/xmin)/(nx - 1));
   for (size_t i(0); i < nx; i++) {
      xx.push_back(xmin*std::exp(xstep*i));
      double arg((std::log(xx.back()) - mu)/m_sigma/std::sqrt(2.));
      cdf.push_back(0.5 + 0.5*::erf(arg));
   }
   
   std::deque<double> my_mus;
   double dxi(1./(nx-1));
   for (size_t i(1); i < nx; i++) {
      double xi(dxi*i);
      std::vector<double>::const_iterator it = 
         std::upper_bound(cdf.begin(), cdf.end(), xi);
      if (it != cdf.end()) {
         size_t indx(it - cdf.begin());
         double xval((xi - cdf.at(indx-1))/(cdf.at(indx) - cdf.at(indx+1))
                     *(xx.at(indx) - xx.at(indx-1)) + xx.at(indx-1));
         my_mus.push_front(std::cos(xval));
      }
   }
   mus.resize(my_mus.size());
   std::copy(my_mus.begin(), my_mus.end(), mus.begin());
}

size_t LogNormalMuDist::energyIndex(double energy) const {
   double logE(std::log(energy));
   if (logE <= m_logEmin) {
      return 1;
   } else if (logE >= m_logEmax) {
      return m_muPoints.size() - 1;
   }
   return static_cast<size_t>((logE - m_logEmin)/m_estep);
}

} // namespace Likelihood
