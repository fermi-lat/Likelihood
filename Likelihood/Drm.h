/**
 * @file Drm.h
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Drm.h,v 1.5 2011/10/18 04:56:55 jchiang Exp $
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
   
   Drm(double ra, double dec, const Observation & observation, 
       const std::vector<double> & ebounds, size_t npts=30);

   void convolve(const std::vector<double> & true_counts,
                 std::vector<double> & meas_counts) const;

   const std::vector<double> & row(size_t k) const {
      return m_drm.at(k);
   }
       
private:

   astro::SkyDir m_dir;
   const Observation & m_observation;
   std::deque<double> m_ebounds;
   size_t m_npts;

   std::vector< std::vector<double> > m_drm;

   void compute_drm();
      
   double matrix_element(double etrue, double emeas_min, 
                         double emeas_max) const;
                         
};

} // namespace Likelihood

#endif // Likelihood_Drm_h
