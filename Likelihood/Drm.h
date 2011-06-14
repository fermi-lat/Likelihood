/**
 * @file Drm.h
 * @brief Detector response matrix for use in convolving model cubes
 * with mean energy dispersion.
 *
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_Drm_h
#define Likelihood_Drm_h

#include <deque>
#include <vector>

namespace Likelihood {

class Observation;

class Drm {

public:
   
   Drm(Observation & observation, const std::vector<double> & ebounds,
       size_t npts=30);

   void convolve(const std::vector<double> & true_counts,
                 std::vector<double> & meas_counts) const;
       
private:

   Observation & m_observation;
   std::deque<double> m_ebounds;
   size_t m_npts;

   std::vector< std::vector<double> > m_drm;

   void compute_drm()

};

} // namespace Likelihood

#endif // Likelihood_Drm_h
