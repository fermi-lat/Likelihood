#include <vector>
#include <string>
#include <cmath>

#include "../Likelihood/Psf.h"

namespace Likelihood {

Psf * Psf::s_instance = 0;

void Psf::readPsfData(const std::string &file, int hdu) {
   enum {FRONT = 2, BACK, COMBINED};

   switch (hdu) {
   case FRONT:
      m_psfData.add_columns("ENERGY THETA SIG1_F SIG2_F W");
      break;
   case BACK:
      m_psfData.add_columns("ENERGY THETA SIG1_B SIG2_B W");
      break;
   case COMBINED:
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
   m_sig1 = m_psfData[2].val;
   m_sig2 = m_psfData[3].val;
   for (int i = 0; i < m_psfData[4].dim; i++) {
      m_wt.push_back(m_psfData[4].val[i]);
   }
}

double Psf::value(astro::SkyDir appDir, double energy, 
		  astro::SkyDir srcDir, double time) {
// angle between photon and source directions
   double separation = appDir.SkyDir::difference(srcDir);

// Compute the index corresponding to the desired time.
// Here we assume the m_scTimes are at regular intervals.
   double tstep = m_scData[1].time - m_scData[0].time;
   int indx;
   indx = static_cast<int>((time - m_scData[0].time)/tstep);

// inclination wrt spacecraft z-axis in degrees
   double inc = srcDir.SkyDir::difference(m_scData[indx].zAxis)*180./M_PI;

   if (inc < 70.) {
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
			
// convert energy to MeV
   energy *= 1e3;

// do a bilinear interpolation on the effective area data
// this is the ugly code from glean (uses unit-offset kludge of NR 1.2)

// find the energy index
   int ie;
   m_hunt(m_psfData[0].val-1, m_psfData[0].dim, energy, &ie);

// kludge to deal with energies outside of the nominal boundaries 
   if (ie == 0) { 
      ie = 1;
   } else if (ie == m_psfData[0].dim) {
      ie = m_psfData[0].dim - 1;
   }
   
// find the theta index
   int it;
   m_hunt(m_psfData[1].val-1, m_psfData[1].dim, inc, &it);
   if (it == 0) it = 1;

// do the interpolations
   double sig1val
      = m_bilinear(m_psfData[0].dim, m_psfData[0].val-1, ie, energy, 
                   m_psfData[1].dim, m_psfData[1].val-1, it, inc, 
                   m_sig1);
   double sig2val
      = m_bilinear(m_psfData[0].dim, m_psfData[0].val-1, ie, energy, 
                   m_psfData[1].dim, m_psfData[1].val-1, it, inc, 
                   m_sig2);

// simply set the weight using the energy index
   double wt = m_wt[ie];

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
