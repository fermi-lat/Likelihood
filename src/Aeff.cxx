/** 
 * @file Aeff.cxx
 * @brief Implementation for LAT effective area class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Aeff.cxx,v 1.10 2003/04/25 18:32:19 jchiang Exp $
 */

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "Likelihood/Aeff.h"
#include "astro/SkyDir.h"

namespace Likelihood {

Aeff * Aeff::s_instance = 0;

void Aeff::readAeffData(const std::string &file, int hdu) {

   switch (hdu) {
   case Front:
      m_aeffData.add_columns("ENERGY THETA AEFF_F");
      break;
   case Back:
      m_aeffData.add_columns("ENERGY THETA AEFF_B");
      break;
   case Combined:
      m_aeffData.add_columns("ENERGY THETA AEFF_C");
      break;
   default:
      std::cerr << "Invalid HDU for Aeff data: " << hdu << std::endl;
      exit(0);
      break;
   }

   m_aeffFile = file;
   m_aeffHdu = hdu;
   m_aeffData.read_FITS_table(m_aeffFile, m_aeffHdu);

   for (int i = 0; i < m_aeffData[0].dim; i++) {
      m_energy.push_back(m_aeffData[0].val[i]);
   }
   for (int i = 0; i < m_aeffData[1].dim; i++) {
      m_theta.push_back(m_aeffData[1].val[i]);
   }
   m_aeff.resize(m_aeffData[2].dim);
   for (int i = 0; i < m_aeffData[2].dim; i++)
      m_aeff[i] = m_aeffData[2].val[i];
}

double Aeff::value(double energy, const astro::SkyDir &dir, double time) {
// Compute the index corresponding to the desired time.
// Here we assume the scData->vec[].times are at regular intervals.
   double tstep = scData->vec[1].time - scData->vec[0].time;
   int indx;
   indx = static_cast<int>((time - scData->vec[0].time)/tstep);

// inclination wrt spacecraft z-axis in degrees
   double inc = dir.SkyDir::difference(scData->vec[indx].zAxis)*180./M_PI;

   if (inc < incMax()) {
      return value(energy, inc);
   } else {
      return 0;
   }
}

double Aeff::value(double energy, double inc) {
// do a bilinear interpolation on the effective area data
   double aeffval = bilinear(m_energy, energy, m_theta, inc, m_aeff);

   return aeffval;
}

Aeff * Aeff::instance() {
   if (s_instance == 0) {
      s_instance = new Aeff();
   }
   return s_instance;
}

} // namespace Likelihood
