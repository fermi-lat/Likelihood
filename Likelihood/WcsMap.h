/**
 * @file WcsMap.h
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/WcsMap.h,v 1.3 2006/03/23 00:21:25 jchiang Exp $
 */

#ifndef Likelihood_WcsMap_h
#define Likelihood_WcsMap_h

#include <vector>

#include "astro/SkyDir.h"

namespace Likelihood {

class BinnedExposure;
class DiffuseSource;
class MeanPsf;

/**
 * @class WcsMap
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/WcsMap.h,v 1.3 2006/03/23 00:21:25 jchiang Exp $
 */

class WcsMap {

public:

   WcsMap(const std::string & filename, const std::string & extension="");

   WcsMap(const DiffuseSource & diffuseSource, double ra, double dec,
          double radius, int npts, double energy=100.,
          const std::string & proj_name="STG", bool use_lb=false);

   ~WcsMap();

   WcsMap(const WcsMap &);

   WcsMap & operator=(const WcsMap &);

   double operator()(const astro::SkyDir & dir) const;

   WcsMap convolve(double energy, const MeanPsf & psf,
                   const BinnedExposure & exposure,
                   bool performConvolution=true) const;

   const std::vector< std::vector<double> > & image() const {
      return m_image;
   }

private:

   astro::SkyDir m_refDir;

   std::vector< std::vector<double> > m_image;

   int m_naxis1;
   int m_naxis2;

   astro::SkyProj * m_proj;

   WcsMap() : m_proj(0) {}

};

} // namespace Likelihood

#endif // Likelihood_WcsMap_h
