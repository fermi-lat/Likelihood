/**
 * @file EquinoxRotation.h
 * @brief Given a direction on the Sky, rotate to a coordinate
 * system that has its equinox at that location.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_EquinoxRotation_h
#define Likelihood_EquinoxRotation_h

#include <vector>

#include "astro/SkyDir.h"

namespace Likelihood {

/**
 * @class EquinoxRotation
 * @brief Class to perform the "Equinox Rotation" described in
 * <a href="http://lheawww.gsfc.nasa.gov/~jchiang/SSC/like_3.ps">
 * LikeMemo 3</a>.
 * @author J. Chiang
 *
 * $Header$
 */

class EquinoxRotation {

public:
     
   EquinoxRotation() {}
     
   EquinoxRotation(double alpha0, double delta0);

   ~EquinoxRotation() {}
  
   void do_rotation(const astro::SkyDir & inDir, astro::SkyDir & outDir,
                    bool reverse=false);
  
   EquinoxRotation * clone() const {
      return new EquinoxRotation(*this);
   }

private:

   std::vector< std::vector<double> > m_rotMatrix;

};                       

} // namespace Likelihood

#endif // Likelihood_EquinoxRotation_h
