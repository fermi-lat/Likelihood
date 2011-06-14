/**
 * @file Drm.h
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Drm.h,v 1.1 2011/06/14 06:31:53 jchiang Exp $
 */

#ifndef Likelihood_Drm_h
#define Likelihood_Drm_h

#include <deque>
#include <vector>

#include "astro/SkyDir.h"

namespace Likelihood {

class Observation;

class Drm {

public:
   
   Drm(Observation & observation, const std::vector<double> & ebounds,
       size_t npts=30);

   void convolve(const std::vector<double> & true_counts,
                 std::vector<double> & meas_counts) const;
       
private:

   astro::SkyDir m_dir;
   Observation & m_observation;
   std::deque<double> m_ebounds;
   size_t m_npts;

   std::vector< std::vector<double> > m_drm;

   void compute_drm()

};

} // namespace Likelihood

#endif // Likelihood_Drm_h
