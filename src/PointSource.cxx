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

namespace Likelihood {

PointSource::PointSource(const PointSource &rhs) : Source(rhs) {
   m_dir = rhs.m_dir;
   m_spectrum = rhs.m_spectrum;
}

double PointSource::fluxDensity(double energy, double time,
                                const astro::SkyDir &dir) const {

// Scale the energy spectrum by the psf value and the effective area
// and convolve with the energy dispersion (now a delta-function in
// energy), all of which are functions of time and spacecraft attitude
// and orbital position.

   Psf *psf = Psf::instance();
   Aeff *aeff = Aeff::instance();

   return (*m_spectrum)(energy)*(*psf)(dir, energy, m_dir.getDir(), time)
      *(*aeff)(energy, dir, time);
}

} // namespace Likelihood
