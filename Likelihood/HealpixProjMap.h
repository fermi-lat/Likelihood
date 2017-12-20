/**
 * @file HealpixProjMap.h
 * @brief A map with reference point centered on the image and that
 * uses astro::ProjBase for indexing its internal representation.
 * @author E. Charles, from J. Chiang WcsMap2
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/HealpixProjMap.h,v 1.1 2015/12/10 00:57:58 echarles Exp $
 */

#ifndef Likelihood_HealpixProjMap_h
#define Likelihood_HealpixProjMap_h

#include <vector>
#include <utility>

#include "astro/SkyDir.h"
#include "astro/ProjBase.h"
#include "astro/HealpixProj.h"
#include "Likelihood/ProjMap.h"
#include "Likelihood/CountsMapBase.h"
#include "healpix_map.h"

namespace Likelihood {

class BinnedExposureBase;
class DiffuseSource;
class MeanPsf;
class CountsMapBase;
class CountsMapHealpix;

/**
 * @class HealpixProjMap
 *
 */

class HealpixProjMap : public ProjMap {

public:

  static void foldVector(const std::vector<float>& vector_in,
			 const astro::HealpixProj& hp_proj, 
			 int npix, int nebins,
			 std::vector<Healpix_Map<float> >& image_out);
   
  static void convertToDifferential(std::vector<Healpix_Map<float> >& image_out,
				    const std::vector<double>& energy_bin_widths,
				    int npix, const double& solid_angle);

  static void convertToIntegral(std::vector<float> image_out,
				const std::vector<Healpix_Map<float> >& image_in,
				const std::vector<double>& energy_bin_widths,
				int npix, const double& solid_angle);

public:

   HealpixProjMap(const std::string & filename, 
		  const std::string & extension,
		  bool interpolate=true, bool enforceEnergyRange=false);
   
   HealpixProjMap(const DiffuseSource & diffuseSource, 
		  int order, Healpix_Ordering_Scheme scheme, 
		  double energy=100.,bool use_lb=false,
		  double radius=180, double ra=0, double dec=0.,
		  bool interpolate=false, bool enforceEnergyRange=false);

   HealpixProjMap(const DiffuseSource & diffuseSource, 
		  int nside, Healpix_Ordering_Scheme scheme, const nside_dummy,
		  double energy=100.,bool use_lb=false,
		  double radius=180, double ra=0, double dec=0.,
		  bool interpolate=false, bool enforceEnergyRange=false);

   HealpixProjMap(const CountsMapHealpix& theMap, 
		  CountsMapBase::ConversionType cType);

   HealpixProjMap(const HealpixProjMap &, bool copy_image=true);

   HealpixProjMap(const HealpixProjMap &, const double& energy, const Healpix_Map<float>& image);

   virtual ~HealpixProjMap();

   virtual HealpixProjMap & operator=(const HealpixProjMap &);

   virtual HealpixProjMap* cast_healpix() { return this; }

   virtual double operator()(const astro::SkyDir & dir, double energy=-1) const;

   virtual double operator()(const astro::SkyDir & dir, int k) const;

   virtual ProjMap* convolveAll(const MeanPsf & psf,
				const BinnedExposureBase * exposure=0,
				bool performConvolution=true) const;

   virtual ProjMap* convolve(double energy, const MeanPsf & psf,
			     const BinnedExposureBase * exposure=0,
			     bool performConvolution=true,
			     int k=0) const;
   
   virtual CountsMapBase* makeCountsMap(const CountsMapBase& counts_map) const;

   virtual double solidAngle(double ilon, double ilat) const { return m_solidAngle; }

   inline double solidAngleHealpix() const { return m_solidAngle; }

   inline const std::vector< Healpix_Map<float> >& image() const { return m_image; }

   /// @return Pixel value as a function of index
   virtual double pixelValue(double ilon, double ilat, int k=0) const;
   
   virtual std::pair<astro::SkyDir, astro::SkyDir> minMaxDistPixels(const astro::SkyDir & dir) const;

   virtual bool insideMap(const astro::SkyDir & dir) const;

   virtual double pixelSize() const { return m_pixelSize; }
   
   /// Rebin the map data by the specified factor.  By default, this
   /// will average over combined pixels, as would be appropriate for
   /// intensity or exposure maps.  For counts maps, set
   /// average=false, so that combined bins are summed.  Return a new
   /// ProjMap object with the new geometry.  The reference direction
   /// will be unchanged from the original (and so will not generally
   /// point to the center of the map.)
   virtual ProjMap* rebin(unsigned int factor, bool average=true);

   // In case our map is of less than the whole sky we will want to 
   // a compact vector representation of the data and will need to 
   // map from global HEALPix index to the index in the compact representation
   int globalToLocal(int glo) const;   
   int localToGlobal(int loc) const;
   int nPixels() const;

protected:

   HealpixProjMap();

   virtual void computeMapIntegrals();

   void latchCacheData();

private:

   typedef Healpix_Map<float> ImagePlane_t;
   std::vector<ImagePlane_t> m_image;

   double m_solidAngle;

   double m_pixelSize;

   void check_negative_pixels(const ImagePlane_t &) const;

   // it is actually the same as m_proj
   // it is just here from convinience
   astro::HealpixProj* m_healpixProj;

};

} // namespace Likelihood

#endif // Likelihood_HealpixProjMap_h
