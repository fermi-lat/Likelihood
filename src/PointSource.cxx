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
#include "Gaussian.h"

namespace Likelihood {

PointSource::PointSource(const PointSource &rhs) : Source(rhs) {
   m_dir = rhs.m_dir;
   m_spectrum = rhs.m_spectrum;
}

double PointSource::fluxDensity(double energy, double time,
                                const astro::SkyDir &dir) const {

   double separation = getSeparation(dir);

// One should scale the energy spectrum by the psf value and the
// effective area and convolve with the energy dispersion, all of
// which are functions of time and spacecraft attitude and orbital
// position.

// For now, we just scale by a jury-rigged Gaussian psf with an ad-hoc
// energy scaling.

   const double sigma0 = 9.*M_PI/180.;       // Gaussian width in radians
   const double E0 = 0.030;                  // at 30 MeV
   const double psf_index = -0.75;
   Gaussian my_psf(1., 0., pow(sigma0*(energy/E0), psf_index));

   return (*m_spectrum)(energy)*my_psf(separation);
}

void PointSource::setDir(double ra, double dec) {
   astro::SkyDir dir(ra, dec);
   m_dir = dir;
}

} // namespace Likelihood
