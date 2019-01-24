/**
 * @file PSFUtils.h
 * @brief Functions to deal with PSF Integration and convolution
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 *  This file contains a number of functions useful for PSF Integration and convolution.
 *
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PSFUtils.h,v 1.8 2017/09/29 01:44:14 echarles Exp $
 */

#ifndef Likelihood_PSFUtils_h
#define Likelihood_PSFUtils_h

#include <vector>
#include "healpix_map.h"

#include "Likelihood/FileUtils.h"

namespace st_stream {
  class StreamFormatter;
}

namespace astro {
  class ProjBase;
  class HealpixProj;
  class SkyDir;
}

namespace Likelihood {

  class Source;
  class CompositeSource;
  class DiffuseSource;
  class PointSource;
  class BinnedCountsCache;
  class BinnedLikeConfig;
  class CountsMapBase;
  class CountsMap;
  class CountsMapHealpix;
  class Drm;
  class ProjMap;
  class WcsMap2;
  class MeanPsf;
  class Observation;
  class BinnedExposureBase;
  class Pixel;
  class PsfIntegConfig;

  namespace PSFUtils {
   
    /* Split at double into the integer and remainder parts, and deal with boundry conditions 

       val  : the input value
       n    : the size of the array
       i    : set to the integer index for val, this will be between 0 and n-2
       rem  : set to the floating point remainder for val ( i.e., rem = val - i )      
     */
    void get_index_and_remainder(double val, size_t n, size_t& i, double& rem); 


    /* A specialized version of bilinear interpolator that works when 
       the x and y arguments are the indicies into a 2D array 
       
       x,y  :  Coordinate indices
       z    :  Grid values

       returns the interpolated value.  
       any calls for values outside the array will return the extrapolated value       
    */
    double bilinear_on_grid(double x, double y,
			    const std::vector< std::vector<double> > &z);    

    /* This is inefficient, it is only here for testing, don't use it */
    double test_bilinear(double x, double y,
			 const std::vector< std::vector<double> > &z);  

    /* Get the maximum distance pixel from a particular direction 

       pixel   : Set of pixels to consider
       dir     : Direction to consider 
       
       returns the maximum angular separation in degrees
    */
    double maxRadius(const std::vector<Likelihood::Pixel> & pixels,
		     const astro::SkyDir & dir);


    /* Get the maximum distance pixel in a counts map from a particular direction.
       This version used the boundry pixels of the map to improve the speed.      

       src     : A point source
       dataMap : A counts map
       
       returns the maximum angular separation in degrees
    */    
    double maxPsfRadius(const PointSource& src,const CountsMap& dataMap);

    
    /* Get the maximum map size needed for a diffuse map.

       src     : A point source
       dataMap : A counts map
       
       returns the smaller of the Source size and the data map size
    */
    double neededMapSize(const DiffuseSource& src,const CountsMapBase& dataMap);



    /* A utility function to convert from WCS to HEALPix.   
       FIXME: should this live somewhere else?
       
       inputMap     : The map we are converting
       proj         : The HEALPix projection we are converting into
       hpm          : Filled with coverted map data
       
    */
    void fillHealpixFromWcsMap(const WcsMap2& inputMap,
			       const astro::HealpixProj& proj,
			       Healpix_Map<float>& hpm);
    

   /* A utility function to convert from WCS to HEALPix.   
      FIXME: should this live somewhere else?

      This version is for partial-sky HEALPix maps.

      inputMap     : The map we are converting
      proj         : The HEALPix projection we are converting into
      pixList      : The list of HEALPix pixels to fill
      hpm          : Filled with coverted map data
 
    */
    void fillHealpixFromWcsMap(const WcsMap2& inputMap,
			       const astro::HealpixProj& proj,
			       const std::vector<int>& pixList,
			       Healpix_Map<float>& hpm);

    /* Build the PSF for a particular direction */
    MeanPsf* build_psf(const Source& src, 
		       const CountsMapBase& dataMap,
		       const std::vector<double>& energies,
		       const Observation& obs);


    /* Test to see if a diffuse source has a MapCubeFuction */
    bool haveMapCubeFunction(DiffuseSource& src);
    
    /* Compute the resampling factor for a diffuse source 

       This compares the diffuse source map pixel size with the counts map pixel size.
     */
    double computeResampFactor(const DiffuseSource & src,
			       const CountsMapBase & dataMap,
			       const PsfIntegConfig& config);
    
    /* Compute the non-cartesian corrections angular separation between the pixels in a 
       counts map and a point source and put them on a vector

       src          : The point source in question.
       dataMap      : The counts map in question
       pixelOffsets : Filled with the non-cartesian corrections.   
           I.e,. the actual angular separations are ( 1 + pixelOffsets ) * ref_pixel_size
    */
    void createOffsetMap(const Source& src,
			 const CountsMap& dataMap,
			 std::vector< std::vector< double > >& pixelOffsets);


    /* Generic function to make a SourceMap model 

       src         : The source in question
       dataMap     : The counts map in question
       energies    : The vector of energies to consider
       meanpsf     : The average PSF across the ROI
       bexpmap     : The Binned Exposure map (actually a cube )
       config      : Parameters for PSF integration 
       srcMapsFile : The file with the source maps for the Composite source
       drm         : For energy dispersion, NULL -> no energy dispersion       
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       kmin        : Minimum energy layer
       kmax        : Maximum energy layer
     */
    int makeModelMap(const Source& src, 
		     const BinnedCountsCache& dataCache,
		     const std::vector<double>& energies,
		     const MeanPsf& meanpsf,
		     const BinnedExposureBase* bexpmap,
		     const PsfIntegConfig& config,	
		     const std::string & srcMapsFile,
		     const Drm* drm,
		     st_stream::StreamFormatter& formatter,
		     std::vector<float>& modelmap,
		     FileUtils::SrcMapType& mapType,
		     int kmin=0, int kmax=-1);
		      
    /* Compute the SourceMap model values for a Diffuse source.  
       This version just calls one of the versions below.

       src         : The source in question
       dataMap     : The counts map in question
       energies    : The vector of energies to consider
       meanpsf     : The average PSF across the ROI
       bexpmap     : The Binned Exposure map (actually a cube )
       config      : Parameters for PSF integration 
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       kmin        : Minimum energy layer
       kmax        : Maximum energy layer
      
       return 0 for success, error code otherwise 
    */
    int makeDiffuseMap(const DiffuseSource& src, 
		       const CountsMapBase& dataMap,
		       const std::vector<double>& energies,
		       const MeanPsf& meanpsf,
		       const BinnedExposureBase* bexpmap,
		       const PsfIntegConfig& config,	
		       st_stream::StreamFormatter& formatter,
		       std::vector<float>& modelmap,
		       FileUtils::SrcMapType& mapType,
		       int kmin=0, int kmax=-1);
		      

    /* Compute the SourceMap model values for a Diffuse source.       
       This version is called for WCS-based Counts Maps

       src         : The source in question
       dataMap     : The counts map in question
       energies    : The vector of energies to consider
       meanpsf     : The average PSF across the ROI
       bexpmap     : The Binned Exposure map (actually a cube )
       config      : Parameters for PSF integration 
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       kmin        : Minimum energy layer
       kmax        : Maximum energy layer
      
       return 0 for success, error code otherwise 
    */
    int makeDiffuseMap_wcs(const DiffuseSource& src, 
			   const CountsMap& dataMap,
			   const std::vector<double>& energies,
			   const MeanPsf& meanpsf,
			   const BinnedExposureBase* bexpmap,
			   const PsfIntegConfig& config,
			   st_stream::StreamFormatter& formatter,
			   std::vector<float>& modelmap,
			   FileUtils::SrcMapType& mapType,
			   int kmin=0, int kmax=-1);
    
    /* Compute the SourceMap model values for a Diffuse source.       
       This version is called for all-sky HEALPix-based Counts Maps

       src         : The source in question
       dataMap     : The counts map in question
       energies    : The vector of energies to consider
       meanpsf     : The average PSF across the ROI
       bexpmap     : The Binned Exposure map (actually a cube )
       config      : Parameters for PSF integration 
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       kmin        : Minimum energy layer
       kmax        : Maximum energy layer
     
       return 0 for success, error code otherwise 
    */
    int makeDiffuseMap_healpix(const DiffuseSource& src, 
			       const CountsMapHealpix& dataMap,
			       const std::vector<double>& energies,
			       const MeanPsf& meanpsf,
			       const BinnedExposureBase* bexpmap,
			       const PsfIntegConfig& config,
			       st_stream::StreamFormatter& formatter,
			       std::vector<float>& modelmap,
			       FileUtils::SrcMapType& mapType,
			       int kmin=0, int kmax=-1);
    
    /* Compute the SourceMap model values for a Diffuse source.       
       This version is called partial-sky HEALPix-based Counts Maps

       src         : The source in question
       dataMap     : The counts map in question
       energies    : The vector of energies to consider
       meanpsf     : The average PSF across the ROI
       bexpmap     : The Binned Exposure map (actually a cube )
       config      : Parameters for PSF integration 
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       kmin        : Minimum energy layer
       kmax        : Maximum energy layer
     
       return 0 for success, error code otherwise 
    */
    int makeDiffuseMap_native(const DiffuseSource& src, 
			      const CountsMapHealpix& dataMap,
			      const std::vector<double>& energies,
			      const MeanPsf& meanpsf,
			      const BinnedExposureBase* bexpmap,
			      const PsfIntegConfig& config,
			      st_stream::StreamFormatter& formatter,
			      std::vector<float>& modelmap,
			      FileUtils::SrcMapType& mapType,
			      int kmin=0, int kmax=-1);
   
    /* Compute the SourceMap model values for a Point source.  
       This version just calls one of the versions below.

       src         : The source in question
       dataMap     : The counts map in question
       energies    : The vector of energies to consider
       config      : Parameters for PSF integration 
       meanpsf     : The average PSF across the ROI
       bexpmap     : The Binned Exposure map (actually a cube )
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       kmin        : Minimum energy layer
       kmax        : Maximum energy layer
       
       return 0 for success, error code otherwise 
    */
    int makePointSourceMap(const PointSource& src, 
			   const CountsMapBase& dataMap,
			   const std::vector<double>& energies,
			   const PsfIntegConfig& config,
			   const MeanPsf& meanpsf,
			   const BinnedExposureBase* bexpmap,
			   st_stream::StreamFormatter& formatter,
			   std::vector<float>& modelmap,
			   FileUtils::SrcMapType& mapType,
			   int kmin=0, int kmax=-1);
    
    /* Compute the SourceMap model values for a Point source.  
       This version is called for WCS-based Counts Maps

       pointSrc    : The source in question
       dataMap     : The counts map in question
       energies    : The vector of energies to consider
       config      : Parameters for PSF integration 
       meanpsf     : The average PSF across the ROI
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       kmin        : Minimum energy layer
       kmax        : Maximum energy layer
    
       return 0 for success, error code otherwise 
    */
    int makePointSourceMap_wcs(const PointSource& pointSrc, 
			       const CountsMap& dataMap,
			       const std::vector<double>& energies,
			       const PsfIntegConfig& config,
			       const MeanPsf& meanpsf,
			       const BinnedExposureBase* bexpmap,
			       st_stream::StreamFormatter& formatter,
			       std::vector<float>& modelmap,
			       FileUtils::SrcMapType& mapType,
			       int kmin=0, int kmax=-1);
    
    /* Compute the SourceMap model values for a Point source.  
       This version is called for HEALPix-based Counts Maps

       pointSrc    : The source in question
       dataMap     : The counts map in question
       energies    : The vector of energies to consider
       config      : Parameters for PSF integration 
       meanpsf     : The average PSF across the ROI
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       kmin        : Minimum energy layer
       kmax        : Maximum energy layer
      
       return 0 for success, error code otherwise 
    */
    int makePointSourceMap_healpix(const PointSource& pointSrc, 
				   const CountsMapHealpix& dataMap,
				   const std::vector<double>& energies,
				   const PsfIntegConfig& config,
				   const MeanPsf& meanpsf,
				   const BinnedExposureBase* bexpmap,
				   st_stream::StreamFormatter& formatter,
				   std::vector<float>& modelmap,
				   FileUtils::SrcMapType& mapType,
				   int kmin=0, int kmax=-1);


    /* Compute the SourceMap model values for a Composite source.  
 
       compSrc     : The source in question
       dataCache   : Cache with the CountsMap to use as template
       energies    : The vector of energies to consider
       srcMapsFile : The file with the source maps for the Composite source
       drm         : For energy dispersion, NULL -> no energy dispersion
       formatter   : Stream for writting progress messages
       modelmap    : Filled with the model values
       mapType     : Enum that specifies how to store the map
       
       return 0 for success, error code otherwise 
    */
    int makeCompositeMap(const CompositeSource& compSrc, 
			 const BinnedCountsCache& dataCache,
			 const std::vector<double>& energies,
			 const std::string & srcMapsFile,
			 const Drm* drm,
			 st_stream::StreamFormatter& formatter,
			 std::vector<float>& modelmap,
			 FileUtils::SrcMapType& mapType,
			 int kmin=0, int kmax=-1);

    /* Get the PSF value for a point source at a particular pixel and energy 

       meanpsf     : The average PSF at the src direction
       energy      : The energy in question
       srcDir      : The location of the poitn source
       pixel       : The pixel in question
       srcCoord    : The coordinates of the source
       pixCoord    : The coordinates of the pixel
       pixelOffsets  : Correction factors for projection
       ref_pixel_size : Size of reference pixel
       config      : Parameters for PSF integration
       
       returns integral of PSF over the pixel in question at the energy in question
    */
    double psfValueEstimate(const MeanPsf & meanPsf, double energy, 
			    const astro::SkyDir & srcDir, 
			    const Pixel & pixel,
			    const std::pair<double, double>& srcCoord,
			    const std::pair<double, double>& pixCoord,
			    const std::vector< std::vector< double > >& pixelOffsets,
			    double ref_pixel_size,
			    const PsfIntegConfig& config);

    /* Get the PSF value for a point source at a particular pixel and energy 
       by taking the value at the pixel center

       meanpsf     : The average PSF at the src direction
       energy      : The energy in question
       srcDir      : The location of the poitn source
       pixel       : The pixel in question

       returns integral of PSF over the pixel in question at the energy in question
    */
    double psfValue_pixelCenter(const MeanPsf & meanPsf, double energy, 
				const astro::SkyDir & srcDir, 
				const Pixel & pixel);
    
    /* Get the PSF value for a point source at a particular pixel and energy 
       by approximating the integral over the pixel using and equal area annulus

       meanpsf     : The average PSF at the src direction
       energy      : The energy in question
       srcDir      : The location of the poitn source
       pixel       : The pixel in question

       returns integral of PSF over the pixel in question at the energy in question
    */
    double psfValue_annular(const MeanPsf & meanPsf, double energy, 
			    const astro::SkyDir & srcDir, 
			    const Pixel & pixel);

    /* Get the PSF value for a point source at a particular pixel and energy 
       by using addaptive integration

       meanpsf     : The average PSF at the src direction
       energy      : The energy in question
       srcDir      : The location of the poitn source
       pixel       : The pixel in question
       srcCoord    : The coordinates of the source
       pixCoord    : The coordinates of the pixel
       pixelOffsets  : Correction factors for projection
       ref_pixel_size : Size of reference pixel
       config      : Parameters for PSF integration

       returns integral of PSF over the pixel in question at the energy in question
    */
    double integrate_psf_adaptive(const MeanPsf & meanPsf, double energy, 
				  const astro::SkyDir & srcDir, 
				  const Pixel & pixel,
				  const std::pair<double, double>& srcCoord,
				  const std::pair<double, double>& pixCoord,
				  const std::vector< std::vector< double > >& pixelOffsets,
				  double ref_pixel_size,
				  const PsfIntegConfig& config);

    /* Get the PSF value for a point source at a particular pixel and energy 
       by integrating over sub-pixels

       meanpsf     : The average PSF at the src direction
       energy      : The energy in question
       srcDir      : The location of the point source
       pixel       : The pixel in question
       npts        : Number for intergration poitns

       returns integral of PSF over the pixel in question at the energy in question
    */  
    double integrate_psf(const MeanPsf & meanPsf, double energy,
			 const astro::SkyDir & srcDir, 
			 const Pixel & pixel, 
			 size_t npts = 11);
    
    /* Get the PSF value for a point source at a particular pixel and energy 
       by integrating over sub-pixels and interpolating projection correction factors

       meanpsf     : The average PSF at the src direction
       energy      : The energy in question
       srcDir      : The location of the point source
       pixel       : The pixel in question
       srcCoord    : The coordinates of the source
       pixCoord    : The coordinates of the pixel
       pixelOffsets  : Correction factors for projection
       ref_pixel_size : Size of reference pixel
       npts        : Number for intergration poitns

       returns integral of PSF over the pixel in question at the energy in question
    */  
    double integrate_psf_fast(const MeanPsf & meanPsf, double energy, 
			      const astro::SkyDir & srcDir, 
			      const Pixel & pixel, 
			      const std::pair<double, double>& srcCoord,
			      const std::pair<double, double>& pixCoord,
			      const std::vector< std::vector< double > >& pixelOffsets,
			      double ref_pixel_size,
			      size_t npts);

  } // namespace PSFUtils

} // namespace Likelihood

#endif // Likelihood_PSFUtils_h
