/**
 * @file WcsMap.h
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/WcsMap.h,v 1.8 2009/02/03 07:24:36 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/WcsMap.h,v 1.8 2009/02/03 07:24:36 jchiang Exp $
 */

class WcsMap {

public:

   WcsMap(const std::string & filename, const std::string & extension="",
          bool interpolate=true);

   WcsMap(const DiffuseSource & diffuseSource, double ra, double dec,
          double pix_size, int npts, double energy=100.,
          const std::string & proj_name="STG", bool use_lb=false,
          bool interpolate=false);

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

   /// @return Solid angle of the (ilon, ilat) pixel
   static double solidAngle(const astro::SkyProj & proj, 
                            double ilon, double ilat);

   double solidAngle(double ilon, double ilat) const;

   /// @return Pixel value as a function index
   double pixelValue(double ilon, double ilat) const;
   
   /// @return SkyDir corresponding to the pixel indices
   astro::SkyDir skyDir(double ilon, double ilat) const;

   int nxpix() const {
      return m_naxis1;
   }

   int nypix() const {
      return m_naxis2;
   }

   bool insideMap(const astro::SkyDir & dir) const;

private:

   astro::SkyDir m_refDir;

   std::vector< std::vector<double> > m_image;

   int m_naxis1;
   int m_naxis2;

   astro::SkyProj * m_proj;
   
   bool m_interpolate;

   bool m_isPeriodic;

   astro::SkyDir::CoordSystem m_coordSys;

   WcsMap();

};

} // namespace Likelihood

#endif // Likelihood_WcsMap_h
