/**
 * @file WcsMap2.h
 * @brief A map with reference point centered on the image and that
 * uses WCS projections for indexing its internal representation.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/WcsMap2.h,v 1.12 2016/09/14 21:47:00 echarles Exp $
 */

#ifndef Likelihood_WcsMap2_h
#define Likelihood_WcsMap2_h

#include <vector>
#include <utility>

// EAC, make a base class for WcsMap2
#include "Likelihood/ProjMap.h"
#include "astro/SkyDir.h"
#include "Likelihood/CountsMapBase.h"

namespace Likelihood {

class BinnedExposureBase;
class DiffuseSource;
class MeanPsf;
class SpatialFunction;
class CountsMapBase;
class CountsMap;
class Pixel;

/**
 * @class WcsMap2
 *
 */

class WcsMap2 : public ProjMap {

public:

   static void foldVector(const std::vector<float>& vector_in,
			  int naxis1, int naxis2, int naxis3, 
			  std::vector<std::vector<std::vector<float> > >& image_out);
			  
   static void fillSolidAngles(const std::vector<Pixel> & pixels,
			       int naxis1, int naxis2, 
			       std::vector<std::vector<float> >& image_out);

  static void convertToDifferential(std::vector<std::vector<std::vector<float> > >& image_out,
				    const std::vector<double>& energy_bin_widths,
				    const std::vector<std::vector<float> >& solid_angles);

  static void convertToIntegral(std::vector<float>& image_out,
				const std::vector<std::vector<std::vector<float> > >& image_in,
				const std::vector<double>& energy_bin_widths,
				const std::vector<std::vector<float> >& solid_angles);

public:

   WcsMap2(const std::string & filename, const std::string & extension="",
           bool interpolate=true, bool enforceEnergyRange=false,
	   bool computeIntegrals=true);

   WcsMap2(const DiffuseSource & diffuseSource, double ra, double dec,
           double pix_size, int npts, double energy=100.,
           const std::string & proj_name="STG", bool use_lb=false,
           bool interpolate=false, bool enforceEnergyRange=false,
	   bool computeIntegrals=true);

   WcsMap2(const DiffuseSource & diffuseSource, double ra, double dec,
           double crpix1, double crpix2, double cdelt1, double cdelt2,
           int naxis1, int naxis2, double energy=100.,
           const std::string & proj_name="STG", bool use_lb=false,
           bool interpolate=false, bool enforceEnergyRange=false,
	   bool computeIntegrals=true);

   WcsMap2(const CountsMap& theMap, CountsMapBase::ConversionType cType);
 
   WcsMap2(const WcsMap2 &, bool copy_image = true);

   WcsMap2(const WcsMap2 &, const double & energy,
	   const std::vector< std::vector<float> >& image);

   virtual ~WcsMap2();

   virtual WcsMap2 & operator=(const WcsMap2 &);

   virtual WcsMap2* cast_wcs() { return this; }

   virtual double operator()(const astro::SkyDir & dir, double energy=-1) const;

   virtual double operator()(const astro::SkyDir & dir, int k) const;

   virtual ProjMap* convolveAll(const MeanPsf & psf,
				const BinnedExposureBase * exposure=0,
				bool performConvolution=true) const;

   virtual ProjMap* convolve(double energy, const MeanPsf & psf,
			     const BinnedExposureBase * exposure=0,
			     bool performConvolution=true,
			     int k=0) const;

   virtual ProjMap* convolve(double energy, const MeanPsf & psf,
			     const BinnedExposureBase * exposure,
			     const SpatialFunction& fn,
			     int k=0) const;   

   virtual CountsMapBase* makeCountsMap(const CountsMapBase& counts_map) const;

   const std::vector< std::vector< std::vector<float> > > & image() const {
      return m_image;
   }


   /// @return Solid angle of the (ilon, ilat) pixel
   static double solidAngle(const astro::ProjBase & proj, 
                            double ilon, double ilat);

   virtual double solidAngle(double ilon, double ilat) const;

   const std::vector< std::vector<float> > & solidAngles() const;

   /// @return Pixel value as a function of index
   virtual double pixelValue(double ilon, double ilat, int k=0) const;
   
   int nxpix() const {
      return m_naxis1;
   }

   int nypix() const {
      return m_naxis2;
   }


   virtual std::pair<astro::SkyDir, astro::SkyDir> minMaxDistPixels(const astro::SkyDir & dir) const;

   virtual bool insideMap(const astro::SkyDir & dir) const;

   virtual double pixelSize() const;

   void getCorners(std::vector<astro::SkyDir> & corners) const;

   double cdelt1() const {
      return m_cdelt1;
   }

   double cdelt2() const {
      return m_cdelt2;
   }

   /// Rebin the map data by the specified factor.  By default, this
   /// will average over combined pixels, as would be appropriate for
   /// intensity or exposure maps.  For counts maps, set
   /// average=false, so that combined bins are summed.  Return a new
   /// WcsMap2 object with the new geometry.  The reference direction
   /// will be unchanged from the original (and so will not generally
   /// point to the center of the map.)
   ProjMap* rebin(unsigned int factor, bool average=true);

   inline double crval1() const { return m_crval1; }
   inline double crval2() const { return m_crval2; }

   inline bool periodic() const { return m_isPeriodic; }

protected:

  virtual void computeMapIntegrals();

private:

//   typedef std::vector< std::vector<double> > Imageplane_t;
   typedef std::vector< std::vector<float> > ImagePlane_t;

   bool m_isPeriodic;

   std::vector<ImagePlane_t> m_image;

   mutable ImagePlane_t m_solidAngles;

   int m_naxes;
   int m_naxis1;
   int m_naxis2;

   /// astro::SkyProj provides almost no introspection, so we store the
   /// map projection locally.
   double m_crpix1, m_crpix2;
   double m_crval1, m_crval2;
   double m_cdelt1, m_cdelt2;
   double m_crota2;
   
   WcsMap2();

   void check_negative_pixels(const ImagePlane_t &) const;

};

} // namespace Likelihood

#endif // Likelihood_WcsMap2_h
