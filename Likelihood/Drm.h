/**
 * @file Drm.h
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/Drm.h,v 1.3 2011/06/14 22:41:49 jchiang Exp $
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
       
private:

   astro::SkyDir m_dir;
   const Observation & m_observation;
   std::deque<double> m_ebounds;
   size_t m_npts;

   std::vector< std::vector<double> > m_drm;

   void compute_drm();

   void get_emeas(size_t kp, std::vector<double> & emeas) const;

   void get_disp(double etrue, const std::vector<double> & emeas,
                 std::vector<double> & disp) const;
      

};

} // namespace Likelihood

#endif // Likelihood_Drm_h
