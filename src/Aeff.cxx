#include <vector>
#include <string>
#include <cmath>

#include "../Likelihood/Aeff.h"

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
   m_aeff = m_aeffData[2].val;
}

double Aeff::value(double energy, astro::SkyDir dir, double time) {
// Compute the index corresponding to the desired time.
// Here we assume the scData->vec[].times are at regular intervals.
   double tstep = scData->vec[1].time - scData->vec[0].time;
   int indx;
   indx = static_cast<int>((time - scData->vec[0].time)/tstep);

// inclination wrt spacecraft z-axis in degrees
   double inc = dir.SkyDir::difference(scData->vec[indx].zAxis)*180./M_PI;

   if (inc < incMax) {
      return value(energy, inc);
   } else {
      return 0;
   }
}

double Aeff::value(double energy, double inc) {
// do a bilinear interpolation on the effective area data
// this is the ugly code from glean (uses unit-offset kludge of NR 1.2)

// find the energy index
   int ie;
   m_hunt(m_aeffData[0].val-1, m_aeffData[0].dim, energy, &ie);

// kludge to deal with energies outside of the nominal boundaries 
   if (ie == 0) { 
//      ie = 1;
      ie = 0;
   } else if (ie == m_aeffData[0].dim) {
      ie = m_aeffData[0].dim - 1;
   }
   
// find the theta index
   int it;
   m_hunt(m_aeffData[1].val-1, m_aeffData[1].dim, inc, &it);
   if (it == 0) it = 1;

   double aeffval 
      = m_bilinear(m_aeffData[0].dim, m_aeffData[0].val-1, ie, energy, 
                   m_aeffData[1].dim, m_aeffData[1].val-1, it, inc, 
                   m_aeff);
   
   return aeffval;
}

Aeff * Aeff::instance() {
   if (s_instance == 0) {
      s_instance = new Aeff();
   }
   return s_instance;
}

} // namespace Likelihood
