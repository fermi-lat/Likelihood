/** @file PointSource.cxx
 * @brief PointSource class implementation
 *
 * $Header:
 */

#include <vector>
#include <string>
#include <cmath>

#include "astro/SkyDir.h"
#include "../Likelihood/PointSource.h"
#include "../Likelihood/Psf.h"
#include "../Likelihood/Aeff.h"
#include "../Likelihood/ScData.h"
#include "../Likelihood/RoiCuts.h"
#include "../Likelihood/dArg.h"

namespace Likelihood {

void PointSource::m_makeEnergyVector(int nee) {
   RoiCuts *roiCuts = RoiCuts::instance();
   
// set up a logrithmic grid of energies for doing the integral over 
// the spectrum
   double emin = (roiCuts->getEnergyCuts()).first;
   double emax = (roiCuts->getEnergyCuts()).second;
   double estep = log(emax/emin)/(nee-1);
   
   m_energies.reserve(nee);
   for (int i = 0; i < nee; i++)
      m_energies.push_back(emin*exp(i*estep));
}

PointSource::PointSource(const PointSource &rhs) : Source(rhs) {
   m_dir = rhs.m_dir;
   m_spectrum = rhs.m_spectrum;
   m_energies = rhs.m_energies;
}

double PointSource::fluxDensity(double energy, double time,
                                const astro::SkyDir &dir) const {

// Scale the energy spectrum by the psf value and the effective area
// and convolve with the energy dispersion (now a delta-function in
// energy), all of which are functions of time and spacecraft attitude
// and orbital position.

   Psf *psf = Psf::instance();
   Aeff *aeff = Aeff::instance();

   dArg energy_arg(energy);
   double spectrum = (*m_spectrum)(energy_arg);
   double psf_val = (*psf)(dir, energy, m_dir.getDir(), time);
   double aeff_val = (*aeff)(energy, dir, time);
   return spectrum*psf_val*aeff_val;      
}

double PointSource::Npred() {
// presently implemented for no psf cuts....

   Aeff *aeff = Aeff::instance();
   ScData *scData = ScData::instance();

   Function *specFunc = m_functions["Spectrum"];

// compute time integral ordinates
   std::vector<double> time_ordinates;
   time_ordinates.reserve(scData->vec.size());
   for (unsigned int k = 0; k < scData->vec.size(); k++) {
// energy integral using trapezoidal rule
      double time_ordinate = 0;
      for (unsigned int j = 0; j < m_energies.size()-1; j++) {
	 dArg ejp1arg(m_energies[j+1]);
	 dArg ejarg(m_energies[j]);
	 double inclination = 
	    scData->vec[k].zAxis.SkyDir::difference(getDir())*180/M_PI;
	 if (inclination < 70.) {
	    time_ordinate += 
	       (specFunc->value(ejp1arg)
		*(*aeff)(m_energies[j+1], inclination)
		+ specFunc->value(ejarg)
		*(*aeff)(m_energies[j], inclination))/2.
		*(m_energies[j+1] - m_energies[j]);
	 } else {
	    time_ordinate = 0.;
	 }
      }
      time_ordinates.push_back(time_ordinate);
   }

// now perform the time integral
   double my_value = 0;
   for (unsigned int k = 0; k < scData->vec.size()-1; k++) {
      my_value += (time_ordinates[k+1] + time_ordinates[k])/2.
	 *(scData->vec[k+1].time - scData->vec[k].time);
   }

   return my_value;   
}

} // namespace Likelihood
