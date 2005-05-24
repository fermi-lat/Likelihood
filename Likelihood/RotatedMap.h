/**
 * @file RotatedMap.h
 * @brief Map-like representation of a diffuse source that is
 * rotated to an equinox-centered coordinate system.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RotatedMap.h,v 1.2 2005/05/23 19:12:51 jchiang Exp $
 */

#ifndef Likelihood_RotatedMap_h
#define Likelihood_RotatedMap_h

#include <vector>

#include "astro/SkyDir.h"

#include "Likelihood/EquinoxRotation.h"

namespace Likelihood {

class BinnedExposure;
class DiffuseSource;
class MeanPsf;

/**
 * @class RotatedMap
 * @brief Map-like representation of a diffuse source that is
 * rotated to an equinox-centered coordinate system.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/RotatedMap.h,v 1.2 2005/05/23 19:12:51 jchiang Exp $
 */

class RotatedMap {

public:

   RotatedMap(const DiffuseSource & diffuseSource, double ra, double dec,
              double radius, int npts);

   double operator()(const astro::SkyDir & dir) const;

   RotatedMap convolve(double energy, const MeanPsf & psf,
                       const BinnedExposure & exposure) const;

   void getUnrotatedMap(std::vector< std::vector<double> > & map,
                        const std::vector<double> & ras,
                        const std::vector<double> & decs,
                        astro::SkyDir::CoordSystem 
                        coordSys=astro::SkyDir::EQUATORIAL) const;
private:

   std::vector< std::vector<double> > m_image;
   std::vector<double> m_lons;
   std::vector<double> m_lats;
   EquinoxRotation m_rot;

   RotatedMap(const std::vector< std::vector<double> > & rotatedImage,
              const std::vector<double> & lons, 
              const std::vector<double> & lats, 
              const EquinoxRotation & rot) 
      : m_image(rotatedImage), m_lons(lons), m_lats(lats), m_rot(rot) {}

   void linearArray(double xmin, double xmax, int npts, 
                    std::vector<double> & x) const;

};

} // namespace Likelihood

#endif // Likelihood_RotatedMap_h
