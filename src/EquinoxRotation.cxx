/**
 * @file EquinoxRotation.cxx
 * @brief Rotate to a coordinate system with equinox at a specified
 * direction.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/EquinoxRotation.cxx,v 1.2 2005/05/23 05:51:26 jchiang Exp $
 */

#include <cmath>

#include "Likelihood/EquinoxRotation.h"

namespace Likelihood {

EquinoxRotation::EquinoxRotation(double alpha0, double delta0) {
     
   m_rotMatrix.clear();
     
// Convert to radians.
   alpha0 *= M_PI/180;
   delta0 *= M_PI/180;

// Build the rotation matrix, using Fortran-like indexing, i.e.,
// m_rotMatrix[row][column].  Note that this is the *transpose* of the
// matrix given in LikeMemo 3.
   
   double ca = cos(alpha0);
   double sa = sin(alpha0);
   double cd = cos(delta0);
   double sd = sin(delta0);

   std::vector<double> row(3);

   row[0] = cd*ca;
   row[1] = -sa;
   row[2] = -sd*ca;
   m_rotMatrix.push_back(row);

   row[0] = cd*sa;
   row[1] = ca;
   row[2] = -sd*sa;
   m_rotMatrix.push_back(row);

   row[0] = sd;
   row[1] = 0;
   row[2] = cd;
   m_rotMatrix.push_back(row);
}

void EquinoxRotation::do_rotation(const astro::SkyDir & inDir,
                                  astro::SkyDir & outDir, bool reverse) const {
   std::vector<double> inVec(3), outVec(3);

   double alpha = inDir.ra()*M_PI/180;
   double delta = inDir.dec()*M_PI/180;
   inVec[0] = cos(delta)*cos(alpha);
   inVec[1] = cos(delta)*sin(alpha);
   inVec[2] = sin(delta);

// Apply the rotation
   for (int i = 0; i < 3; i++) {
      outVec[i] = 0;
      for (int j = 0; j < 3; j++) {
         if (reverse) {
            outVec[i] += m_rotMatrix[j][i]*inVec[j];
         } else {
            outVec[i] += m_rotMatrix[i][j]*inVec[j];
         }
      }
   }
   outDir = astro::SkyDir(CLHEP::Hep3Vector(outVec[0], outVec[1], outVec[2]));
}

} // namespace Likelihood 
