/** 
 * @file Psf.cxx
 * @brief Implementation for the LAT Point-Spread Function class
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Psf.cxx,v 1.11 2003/04/25 21:51:29 burnett Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <valarray>

#include "astro/SkyDir.h"
#include "Likelihood/Psf.h"

namespace Likelihood {

Psf * Psf::s_instance = 0;

void Psf::readPsfData(const std::string &file, int hdu) {
   switch (hdu) {
   case Front:
      m_psfData.add_columns("ENERGY THETA SIG1_F SIG2_F W");
      break;
   case Back:
      m_psfData.add_columns("ENERGY THETA SIG1_B SIG2_B W");
      break;
   case Combined:
      m_psfData.add_columns("ENERGY THETA SIG1_C SIG2_C W");
      break;
   default:
      std::cerr << "Invalid HDU for PSF data: " << hdu << std::endl;
      exit(0);
      break;
   }

   m_psfFile = file;
   m_psfHdu = hdu;
   m_psfData.read_FITS_table(m_psfFile, m_psfHdu);

   for (int i = 0; i < m_psfData[0].dim; i++) {
      m_energy.push_back(m_psfData[0].val[i]);
   }
   for (int i = 0; i < m_psfData[1].dim; i++) {
      m_theta.push_back(m_psfData[1].val[i]);
   }
//     m_sig1 = m_psfData[2].val;
//     m_sig2 = m_psfData[3].val;
   m_sig1.resize(m_psfData[2].dim);
   for (int i = 0; i < m_psfData[2].dim; i++) {
      m_sig1[i] = m_psfData[2].val[i];
   }
   m_sig2.resize(m_psfData[3].dim);
   for (int i = 0; i < m_psfData[3].dim; i++) {
      m_sig2[i] = m_psfData[3].val[i];
   }
   for (int i = 0; i < m_psfData[4].dim; i++) {
      m_wt.push_back(m_psfData[4].val[i]);
   }
}

double Psf::value(const astro::SkyDir &appDir, double energy, 
                  const astro::SkyDir &srcDir, double time) {
// angle between photon and source directions
   double separation = appDir.SkyDir::difference(srcDir);

// Compute the index corresponding to the desired time.
// Here we assume the scData->vec[].time are at regular intervals.

   double tstep = scData->vec[1].time - scData->vec[0].time;
   int indx;
   indx = static_cast<int>((time - scData->vec[0].time)/tstep);

// inclination wrt spacecraft z-axis in degrees
   double inc = srcDir.SkyDir::difference(scData->vec[indx].zAxis)*180./M_PI;

   if (inc < incMax()) {
      return value(separation, energy, inc);
   } else {
      return 0;
   }
}

double Psf::value(double separation, double energy, double inc) {
   std::vector<double> psf_params;
   
   fillPsfParams(energy, inc, psf_params);

   double sig1 = psf_params[0];
   double sig2 = psf_params[1];
   double wt = psf_params[2];

// compute the psf as the weighted sum of two Gaussians projected onto
// the Celestial sphere (so cannot use the Euclidean space Gaussian
// function)
   sig1 *= M_PI/180.;
   sig2 *= M_PI/180.;

   double mu = cos(separation);
   double part1 = wt*exp(-(1. - mu)/sig1/sig1)
      /sig1/sig1/(1. - exp(-2./sig1/sig1));
   double part2 = (1. - wt)*exp(-(1. - mu)/sig2/sig2)
      /sig2/sig2/(1. - exp(-2./sig2/sig2));

   return (part1 + part2)/2./M_PI;
}

void Psf::fillPsfParams(double energy, double inc, 
                        std::vector<double> &psf_params) {

// do a bilinear interpolation on the effective area data
   double sig1val = bilinear(m_energy, energy, m_theta, inc, m_sig1);
   double sig2val = bilinear(m_energy, energy, m_theta, inc, m_sig2);

// simply set the weight using the upper bound energy
   std::vector<double>::const_iterator ie;
   if (energy < *(m_energy.begin())) {
      ie = m_energy.begin();
   } else if (energy >= *(m_energy.end() - 1)) {
      ie = m_energy.end() - 1;
   } else {
	   ie = std::upper_bound(m_energy.begin(), m_energy.end(), energy);
   }
   double wt = m_wt[ie - m_energy.begin()];

   psf_params.push_back(sig1val);
   psf_params.push_back(sig2val);
   psf_params.push_back(wt);
}

Psf * Psf::instance() {
   if (s_instance == 0) {
      s_instance = new Psf();
   }
   return s_instance;
}

} // namespace Likelihood
