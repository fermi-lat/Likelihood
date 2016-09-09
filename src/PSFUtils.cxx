/**
 * @file Convolve.cxx
 * @brief Functions to  deal with PSF Integration and convolution
 * @author E. Charles, from code in SourceMap by J. Chiang and M. Wood.
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FitUtils.cxx,v 1.22 2016/07/11 23:44:01 mdwood Exp $
 */

#include "Likelihood/PSFUtils.h"

#include <algorithm>
#include <fstream>

#include "astro/SkyProj.h"

#include "st_stream/StreamFormatter.h"
#include "st_facilities/Util.h"

#include "Likelihood/BinnedConfig.h"
#include "Likelihood/BinnedExposureBase.h"
#include "Likelihood/CountsMapBase.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/ConvolveHealpix.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SpatialFunction.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/WcsMap2.h"
#include "Likelihood/HealpixProjMap.h"

namespace Likelihood {

  namespace PSFUtils {

    void get_index_and_remainder(double val, size_t n, size_t& i, double& rem) {      
      // This truncates the values between 0 and n-1
      double trunc = val < 0 ? 0. : val > n - 1 ? n - 1 : val;
      // This returns the integer less than or equal to the truncated value
      // Unless the truncated value is n-1, in which case it returns n-2
      i = std::min( size_t(std::floor(trunc)),n-2);
     
      // This sets the remainder to reproduce the behavior of st_facilities::Utils::bilinear
      // This will extrapolate values outside of the array
      rem = val - i;

    }

    double bilinear_on_grid(double x, double y,
			    const std::vector< std::vector<double> > &z) {
      // Get the grid size from the vector dimensions
      size_t nx = z.size();
      size_t ny = z[0].size();
      
      // Storage for the index and remainder in each dimension
      size_t ix(0);
      size_t iy(0);
      double rx(0.);
      double ry(0.);
      get_index_and_remainder(x,nx,ix,rx);
      get_index_and_remainder(y,ny,iy,ry);
  
      // The complements of the remainders (i.e., 0 -> 1 and 1 -> 0)
      double sx(1-rx);
      double sy(1-ry);

      // The values at the nearest points
      double z1 = z[ix][iy];
      double z2 = z[ix][iy+1];
      double z3 = z[ix+1][iy];
      double z4 = z[ix+1][iy+1];

      // The interpolated value
      double val = (sx*sy*z1) + (sx*ry*z2) + (rx*sy*z3) + (rx*ry*z4);
      return val;      
   }    
   

    double test_bilinear(double x, double y,
			 const std::vector< std::vector<double> > &z) {
      size_t nx = z.size();
      size_t ny = z[0].size();
      std::vector<double> xx(nx);
      std::vector<double> yy(ny);
      for ( size_t ix(0); ix < nx; ix++ ) {
	xx[ix] = ix+1;
      }
      for ( size_t iy(0); iy < ny; iy++ ) {
	yy[iy] = iy+1;
      }
      return st_facilities::Util::bilinear(xx,x,yy,y,z);
    } 


    double maxRadius(const std::vector<Likelihood::Pixel> & pixels,
		     const astro::SkyDir & dir) {
      std::vector<Likelihood::Pixel>::const_iterator pixel = pixels.begin();
      double maxValue(0);
      for ( ; pixel != pixels.end(); ++pixel) {
         double dist = pixel->dir().difference(dir);
         if (dist > maxValue) {
            maxValue = dist;
         }
      }
      return maxValue*180./M_PI;
    }
    
    double maxPsfRadius(const PointSource& src, const CountsMap& dataMap) {
      std::vector<astro::SkyDir> pixelDirs;
      dataMap.getBoundaryPixelDirs(pixelDirs);
      
      const astro::SkyDir & srcDir = src.getDir();
      double radius = srcDir.difference(pixelDirs.at(0));
      for (unsigned int i = 1; i < pixelDirs.size(); i++) {
	double new_rad = srcDir.difference(pixelDirs.at(i));
	if (new_rad < radius ) {
	  radius = new_rad;
	}
      }
      return radius*180./M_PI;
    }


    void fillHealpixFromWcsMap(const WcsMap2& inputMap,
			       const astro::HealpixProj& proj,
			       const std::vector<int>& pixList,
			       Healpix_Map<float>& hpm) {
      // EAC, here we loop over the pixel in the partial-sky HEALPix map
      // and look up the value from the interpolated map
      WcsMap2 nc_map = const_cast<WcsMap2&>(inputMap);
      nc_map.setInterpolation(true);
      int index(0);
      int energyLayer(0);
      hpm.fill(0.);
      for ( std::vector<int>::const_iterator itr = pixList.begin(); itr != pixList.end(); 
	    itr++, index++ ) {
	astro::SkyDir dir( double(*itr), 0., proj );
	double value = inputMap(dir,energyLayer);
	hpm[*itr] = value;
      }  
    }


    void fillHealpixFromWcsMap(const WcsMap2& inputMap,
			       const astro::HealpixProj& proj,
			       Healpix_Map<float>& hpm) {
      // EAC, here we are assuming that the resolution of the WCS-based map is 
      // signficantly higher than the resolution of the HEALPix-based map that we are filling.
      // So we just loop over the WCS map pixels and fill the correspond HEALPix map pixels,
      // taking into account the difference in pixel solid angles.  
      hpm.fill(0.);
      const std::vector< std::vector<float> >& solidAngles_in = inputMap.solidAngles();
      const double solidAngle_out = astro::radToDeg( astro::radToDeg( ASTRO_4PI / float(hpm.Npix() ) ) );
      double xval(1.0);   
      for ( int i(0); i < inputMap.nxpix(); i++, xval += 1. ) {
	double yval(1.0); 
	for ( int j(0); j < inputMap.nypix(); j++, yval += 1. ) {
	  std::pair<double,double> crds = inputMap.getProj()->pix2sph(xval,yval);
	  std::pair<double,double> converted = proj.sph2pix(crds.first,crds.second);
	  int fillPixel = int(converted.first);
	  double solidAngleRatio = solidAngles_in[i][j] / solidAngle_out;
	  hpm[fillPixel] += (inputMap.pixelValue(xval,yval));
	}
      }
    }

    MeanPsf* build_psf(const Source& src, 
		       const CountsMapBase& dataMap,
		       const Observation& obs) {
      const PointSource& pointSrc = dynamic_cast<const PointSource&>(src);
      std::vector<double> energies;
      dataMap.getEnergies(energies);

      const astro::SkyDir & dir(pointSrc.getDir());
      MeanPsf* meanPsf = new MeanPsf(dir.ra(), dir.dec(), energies, obs);
      return meanPsf; 
    } 

    bool haveMapCubeFunction(DiffuseSource& src) {
      Source::FuncMap & srcFuncs = src.getSrcFuncs();
      return srcFuncs["SpatialDist"]->genericName() == "MapCubeFunction";
    }

   
    double computeResampFactor(const DiffuseSource & src,
			       const CountsMapBase & dataMap) {
      double data_pixel_size = dataMap.pixelSize();
      double model_pixel_size = data_pixel_size;
      try {
	model_pixel_size = src.mapBaseObject()->projmap().pixelSize();
      } catch (MapBaseException &) {
	// do nothing
      }
      double resamp_factor = 
	std::max(2, static_cast<int>(data_pixel_size/model_pixel_size));
      return resamp_factor;
    }
    
    void createOffsetMap(const Source& src,
			 const CountsMap& dataMap,
			 std::vector< std::vector< double > >& pixelOffsets) {
      
      const PointSource & pointSrc = dynamic_cast<const PointSource&>(src);
      const astro::SkyDir & dir = pointSrc.getDir();
      const std::vector<Pixel> & pixels = dataMap.pixels();
      const double pixel_size = dataMap.pixelSize();

      pixelOffsets.resize(dataMap.naxis2(),
			  std::vector<double>(dataMap.naxis1()));

      bool galactic(dataMap.isGalactic());
      /// Get the pixel center in pixel coordinates
      double src_lon = galactic ? dir.l() : dir.ra();
      double src_lat = galactic ? dir.b() : dir.dec();

      std::pair<double, double> src_coords(dataMap.projection().sph2pix(src_lon, 
									src_lat));

      std::vector<Pixel>::const_iterator pixel(pixels.begin());
      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
	
	int ix = j/dataMap.naxis1();
	int iy = j%dataMap.naxis1();
	
	double pix_sep = 
	  pixel_size*std::sqrt(std::pow(src_coords.first-(iy+1),2) +  
			       std::pow(src_coords.second-(ix+1),2));
	double ang_sep = dir.difference(pixel->dir())*180./M_PI;

	if(pix_sep > 1E-6) {
	  pixelOffsets[ix][iy] = ang_sep/pix_sep-1.0;
	} else {
	  pixelOffsets[ix][iy] = 0.0;
	}
      }
    }

    int makeDiffuseMap(const Source& src, 
		       const CountsMapBase& dataMap,
		       const MeanPsf& meanpsf,
		       const BinnedExposureBase & bexpmap,
		       const PsfIntegConfig& config,
		       st_stream::StreamFormatter& formatter,
		       std::vector<float>& modelmap) {

      switch ( dataMap.projection().method() ) {
      case astro::ProjBase::WCS:
	return makeDiffuseMap_wcs(src,static_cast<const CountsMap&>(dataMap),meanpsf,bexpmap,config,formatter,modelmap);
      case astro::ProjBase::HEALPIX:
	if ( dataMap.mapRadius() > 45 ) {
	  // The ROI takes up a large fraction of the sky, 
	  // let's do the convolution using HEALPix
	  return makeDiffuseMap_healpix(src,static_cast<const CountsMapHealpix&>(dataMap),meanpsf,bexpmap,config,formatter,modelmap);
	} else {
	  // The ROI only take up a relatively small fraction of the sky (< 15%)
	  // we will do the convlution in the native projection, then convert to healpix
	  return makeDiffuseMap_native(src,static_cast<const CountsMapHealpix&>(dataMap),meanpsf,bexpmap,config,formatter,modelmap);
	}
      default:
	break;
      }
      std::string errMsg("Unrecognized projection type for map object for source: ");
      errMsg += src.getName();
      throw std::runtime_error(errMsg);
      return 1; 
    }
    
    int makeDiffuseMap_wcs(const Source& src, 
			   const CountsMap& dataMap,
			   const MeanPsf& meanpsf,
			   const BinnedExposureBase & bexpmap,
			   const PsfIntegConfig& config,
			   st_stream::StreamFormatter& formatter,
			   std::vector<float>& modelmap) {

      const DiffuseSource& diffuseSrc = dynamic_cast<const DiffuseSource&>(src);
   
      bool haveSpatialFunction = dynamic_cast<const SpatialFunction *>(diffuseSrc.spatialDist()) != 0;

      // If the diffuse source is represented by an underlying map, then
      // rebin according to the minimum bin size.
      try {
	const MapBase & tmp = *diffuseSrc.mapBaseObject();
	double pixSize = tmp.projmap().pixelSize();
	if ( pixSize < config.minbinsz() ) {
	  unsigned int factor = static_cast<unsigned int>( config.minbinsz()/pixSize);
	  formatter.info(4) << "\nrebinning factor: " 
			    << factor << std::endl;
	  if (factor > 1) {
            tmp.rebin(factor);
	  }
	}
      } catch (MapBaseException &) {
	// do nothing
      }
      
      const std::vector<Pixel>& pixels = dataMap.pixels();
      std::vector<double> energies;
      dataMap.getEnergies(energies);
      
      long npts = energies.size()*pixels.size();
      modelmap.resize(npts, 0);
      
      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      
      const astro::SkyDir & mapRefDir = dataMap.refDir();
      double resamp_fact = 1;
      if ( config.resample()) {
	resamp_fact = std::max(config.resamp_factor(), 
			       computeResampFactor(diffuseSrc, dataMap));
      }
      formatter.info(4) << "\nresampling factor: " 
			<< resamp_fact << std::endl;
      double crpix1, crpix2;
      int naxis1, naxis2;
      double cdelt1 = dataMap.cdelt1()/resamp_fact;
      double cdelt2 = dataMap.cdelt2()/resamp_fact;
      size_t nx_offset(0), ny_offset(0);
      size_t nx_offset_upper(0), ny_offset_upper(0);
      if (haveSpatialFunction) {
	naxis1 = static_cast<int>(dataMap.naxis1()*resamp_fact);
	naxis2 = static_cast<int>(dataMap.naxis2()*resamp_fact);      
	crpix1 = resamp_fact*(dataMap.crpix1() - 0.5)+0.5;
	crpix2 = resamp_fact*(dataMap.crpix2() - 0.5)+0.5;  
      } else if (dataMap.conformingMap()) {
	double radius = std::min(180., maxRadius(pixels, mapRefDir) + 10.);
	// Conforming maps have abs(CDELT1) == abs(CDELT2).  This
	// expression for the mapsize ensures that the number of
	// pixels in each dimension is even.
	int mapsize(2*static_cast<int>(radius/std::fabs(cdelt1)));
	formatter.info(4) << "mapsize: " << mapsize << std::endl;
	naxis1 = mapsize;
	naxis2 = mapsize;
	crpix1 = (naxis1 + 1.)/2.;
	crpix2 = (naxis2 + 1.)/2.;
	nx_offset = 
	  static_cast<size_t>((mapsize - dataMap.naxis1()*resamp_fact)/2);
	ny_offset = 
	  static_cast<size_t>((mapsize - dataMap.naxis2()*resamp_fact)/2);
	nx_offset_upper = 
	  static_cast<size_t>((mapsize - dataMap.naxis1()*resamp_fact)/2);
	ny_offset_upper = 
	  static_cast<size_t>((mapsize - dataMap.naxis2()*resamp_fact)/2);
	/// For cases where the resampling factor is an odd number, 
	/// there may be a row or column of pixels not accounted for
	/// by n[xy]_offset.  Here we add that row or column back in if
	/// it is missing.
	int xtest = static_cast<int>((naxis1 - nx_offset - nx_offset_upper) 
				     - dataMap.naxis1()*resamp_fact);
	if (config.resample() && xtest != 0) {
	  nx_offset += 1;
	}
	int ytest = static_cast<int>((naxis2 - ny_offset - ny_offset_upper) 
				     - dataMap.naxis2()*resamp_fact);
	if (config.resample() && ytest != 0) {
	  ny_offset += 1;
	}
	if (!config.resample()) { 
	  // Use integer or half-integer reference pixel based on
	  // input counts map, even though naxis1 and naxis2 both
	  // must be even.
	  if (dataMap.naxis1() % 2 == 1) {
            crpix1 += 0.5;
            nx_offset += 1;
	  }
	  if (dataMap.naxis2() % 2 == 1) {
            crpix2 += 0.5;
            ny_offset += 1;
	  }
	}
      } else {
	// The counts map was not created by gtbin, so just adopt the
	// map geometry without adding padding for psf leakage since
	// this cannot be done in general without redefining the
	// reference pixel and reference direction.
	naxis1 = static_cast<int>(dataMap.naxis1()*resamp_fact);
	naxis2 = static_cast<int>(dataMap.naxis2()*resamp_fact);
	// Ensure an even number of pixels in each direction.
	if (naxis1 % 2 == 1) {
	  naxis1 += 1;
	}
	if (naxis2 % 2 == 1) {
	  naxis2 += 1;
	}
	crpix1 = dataMap.crpix1()*resamp_fact;
	crpix2 = dataMap.crpix2()*resamp_fact;
	nx_offset_upper += 1;
	ny_offset_upper += 1;
      }
            
      size_t counter(0);
      std::vector<double>::const_iterator energy = energies.begin();
      for (int k(0); energy != energies.end(); ++energy, k++) {
	bool interpolate(true);	
	WcsMap2 * convolvedMap(0);      
	WcsMap2 diffuseMap(diffuseSrc,  
			   (dataMap.projection().isGalactic() ? mapRefDir.l() : mapRefDir.ra()), 
			   (dataMap.projection().isGalactic() ? mapRefDir.b() : mapRefDir.dec()),
			   crpix1, crpix2, cdelt1, cdelt2, naxis1, naxis2,
			   *energy, dataMap.proj_name(), 
			   dataMap.projection().isGalactic(), 
			   interpolate);
	
	if(haveSpatialFunction) {
	  const SpatialFunction* m = 
	    dynamic_cast<const SpatialFunction *>(diffuseSrc.spatialDist());
	  convolvedMap = static_cast<WcsMap2*>(diffuseMap.convolve(*energy, meanpsf, 
								   bexpmap, *m));
	} else {
	  convolvedMap = static_cast<WcsMap2*>(diffuseMap.convolve(*energy, meanpsf, 
								   bexpmap, config.performConvolution() ) );
	}
	
	size_t rfac(static_cast<size_t>(resamp_fact));
	double solid_angle;
	double added(0.0);
	for (size_t j(ny_offset); j < naxis2 - ny_offset_upper; j++) {
	  for (size_t i(nx_offset); i < naxis1 - nx_offset_upper; i++) {
            if ((i % rfac == 0) && (j % rfac == 0)) {
	      counter++;
	      if (config.verbose() && (counter % (npts/20)) == 0) {
		formatter.warn() << ".";
	      }
            }
            size_t pix_index = ((j-ny_offset)/rfac)*dataMap.naxis1() 
	      + ((i-nx_offset)/rfac);
            solid_angle = pixels.at(pix_index).solidAngle();
            size_t indx = k*dataMap.naxis1()*dataMap.naxis2() + pix_index;
            modelmap[indx] += (convolvedMap->image()[0][j][i]
			       /resamp_fact/resamp_fact
			       *solid_angle);
	    added += modelmap[indx];
	  }
	}
	delete convolvedMap;
      }
      //   computeNpredArray();
      // Delete model map for map-based diffuse sources to save memory.  The
      // map will be reloaded dynamically if it is needed again.
      try {
	MapBase * mapBaseObj = 
	  const_cast<MapBase *>(diffuseSrc.mapBaseObject());
	mapBaseObj->deleteMap();
	formatter.info(4) << "SourceMap::makeDiffuseSource: "
			  << "called mapBaseObj->deleteMap()"
			  << std::endl;
      } catch (MapBaseException & eObj) {
	// Not a map-based source, so do nothing.
      }
      return 0;
    }
    
    int makeDiffuseMap_healpix(const Source& src, 
			       const CountsMapHealpix& dataMap,
			       const MeanPsf& meanpsf,
			       const BinnedExposureBase & bexpmap,
			       const PsfIntegConfig& config,
			       st_stream::StreamFormatter& formatter,			       
			       std::vector<float>& modelmap) {

      const DiffuseSource& diffuseSrc = dynamic_cast<const DiffuseSource&>(src);
      // If the diffuse source is represented by an underlying map, then
      // rebin according to the minimum bin size.
      try {
	const MapBase & tmp = *diffuseSrc.mapBaseObject();
	double pixSize = tmp.projmap().pixelSize();
	if ( pixSize < config.minbinsz() ) {   
	  unsigned int factor = static_cast<unsigned int>(config.minbinsz()/pixSize);
	  formatter.info(4) << "\nrebinning factor: " 
			    << factor << std::endl;
	  if (factor > 1) {
	    tmp.rebin(factor);
	  }
	}
      } catch (MapBaseException &) {
	// do nothing
      }
      
      double solidAngle = dataMap.solidAngle();
      std::vector<double> energies;
      dataMap.getEnergies(energies);
  
      long npts = energies.size()*dataMap.nPixels();
      modelmap.resize(npts, 0);
      const astro::SkyDir & mapRefDir = dataMap.refDir();
      double resamp_fact = 1.;
      if ( config.resample() ) {
	resamp_fact = std::max(config.resamp_factor(),
			       computeResampFactor(diffuseSrc, dataMap));
      }
      // We need to move the resampling factor to the nearest power of two.
      int resamp_test(1);
      while ( resamp_test < resamp_fact ) {
	resamp_test *= 2;
      }
      resamp_fact = resamp_test;
      
      Healpix_Ordering_Scheme scheme = dataMap.healpixProj()->healpix().Scheme();
      int nside_orig = dataMap.healpixProj()->healpix().Nside(); 
      int resamp_nside = resamp_fact*nside_orig;
      size_t counter(0);
      static const double ALLSKY_RADIUS(180.);
      bool interpolate(false);
      std::vector<double>::const_iterator energy = energies.begin();
      double scaleFactor = solidAngle/(resamp_fact*resamp_fact);
      // This is the index in the output vector, it does _not_ get reset between energy layers
      int outidx(0);
      for (int k(0); energy != energies.end(); ++energy, k++) {
	formatter.warn() << ".";
	HealpixProjMap diffuseMap(diffuseSrc, resamp_nside,
				  scheme,SET_NSIDE,
				  *energy,dataMap.projection().isGalactic(),
				  ALLSKY_RADIUS,mapRefDir.ra(), mapRefDir.dec(),
				  interpolate);
	ProjMap* cmap = diffuseMap.convolve(*energy,meanpsf,bexpmap,config.performConvolution());
	HealpixProjMap* convolvedMap = static_cast<HealpixProjMap*>(cmap);
	Healpix_Map<float> outmap(nside_orig,scheme,SET_NSIDE);
	outmap.Import_degrade(convolvedMap->image()[0]);
	if ( nside_orig == resamp_nside ) {
	  outmap = convolvedMap->image()[0];
	} else {
	  outmap.Import_degrade(convolvedMap->image()[0]);
	}
	for (size_t i(0); i < dataMap.nPixels(); i++,outidx++) {
	  int glo = dataMap.localToGlobalIndex(i);
	  modelmap[outidx] = outmap[glo];
	}
      }
      try {
	MapBase * mapBaseObj = 
	  const_cast<MapBase *>(diffuseSrc.mapBaseObject());
	mapBaseObj->deleteMap();
	formatter.info(4) << "SourceMap::makeDiffuseSource: "
			  << "called mapBaseObj->deleteMap()"
			  << std::endl;
      } catch (MapBaseException & eObj) {
	// Not a map-based source, so do nothing.
      }
      return 0;
    }
    
   
    int makeDiffuseMap_native(const Source& src, 
			      const CountsMapHealpix& dataMap,
			      const MeanPsf& meanpsf,
			      const BinnedExposureBase & bexpmap,
			      const PsfIntegConfig& config,
			      st_stream::StreamFormatter& formatter,
			      std::vector<float>& modelmap) {

      const DiffuseSource& diffuseSrc = dynamic_cast<const DiffuseSource&>(src);

      bool haveSpatialFunction = dynamic_cast<const SpatialFunction *>(diffuseSrc.spatialDist()) != 0;
      // If the diffuse source is represented by an underlying map, then
      // rebin according to the minimum bin size.
      std::string nativeProj("STG");
      try {
	const MapBase & tmp = *diffuseSrc.mapBaseObject();
	nativeProj = tmp.projmap().getProj()->projType();
	double pixSize = tmp.projmap().pixelSize();
	if ( pixSize < config.minbinsz() ) {
	  unsigned int factor = static_cast<unsigned int>(config.minbinsz()/pixSize);
	  formatter.info(4) << "\nrebinning factor: " 
			    << factor << std::endl;
	  if (factor > 1) {
	    tmp.rebin(factor);
	  }
	}
      } catch (MapBaseException &) {
	// do nothing
      }

      const std::vector<Pixel> & pixels = dataMap.pixels();
      std::vector<double> energies;
      dataMap.getEnergies(energies);

      long npts = energies.size()*pixels.size();
      modelmap.resize(npts, 0);
      
      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      
      const astro::SkyDir & mapRefDir = dataMap.refDir();
      double resamp_fact = 1.;
      if ( config.resample() ) {
	resamp_fact = std::max(config.resamp_factor(), 
			       computeResampFactor(diffuseSrc, dataMap));
      }
      formatter.info(4) << "\nresampling factor: " 
			<< resamp_fact << std::endl;
      
      double crpix1, crpix2;
      int naxis1(0), naxis2(0);
      double cdelt1 = -1.*dataMap.pixelSize()/resamp_fact;
      double cdelt2 = dataMap.pixelSize()/resamp_fact;
      size_t nx_offset(0), ny_offset(0);
      size_t nx_offset_upper(0), ny_offset_upper(0);
      double ref1 = dataMap.isGalactic() ? mapRefDir.l() : mapRefDir.ra();
      double ref2 = dataMap.isGalactic() ? mapRefDir.b() : mapRefDir.dec();

      // For spatial functions the map only needs to be big enough to
      // encompass the ROI
      if (haveSpatialFunction) {
	double radius = std::min(180., maxRadius(pixels, mapRefDir) + std::fabs(cdelt1) );
	int mapsize(2*static_cast<int>(radius/std::fabs(cdelt1)));
	naxis1 = mapsize;
	naxis2 = mapsize;
	crpix1 = (mapsize+1.0)*0.5;
	crpix2 = (mapsize+1.0)*0.5;
      } else {
	
	double radius = std::min(180., maxRadius(pixels, mapRefDir) + 10.);
	// Conforming maps have abs(CDELT1) == abs(CDELT2).  This
	// expression for the mapsize ensures that the number of
	// pixels in each dimension is even.
	int mapsize(2*static_cast<int>(radius/std::fabs(cdelt1)));
	
	formatter.info(4) << "mapsize: " << mapsize << std::endl;
	naxis1 = mapsize;
	naxis2 = mapsize;
	crpix1 = (naxis1 + 1.)/2.;
	crpix2 = (naxis2 + 1.)/2.;
	
	int naxis_data = int( dataMap.mapRadius() / dataMap.pixelSize() );
	
	nx_offset = 
	  static_cast<size_t>((mapsize - naxis_data*resamp_fact)/2);
	ny_offset = 
	  static_cast<size_t>((mapsize - naxis_data*resamp_fact)/2);
	nx_offset_upper = 
	  static_cast<size_t>((mapsize - naxis_data*resamp_fact)/2);
	ny_offset_upper = 
	  static_cast<size_t>((mapsize - naxis_data*resamp_fact)/2);
	/// For cases where the resampling factor is an odd number, 
	/// there may be a row or column of pixels not accounted for
	/// by n[xy]_offset.  Here we add that row or column back in if
	/// it is missing.
	int xtest = static_cast<int>((naxis1 - nx_offset - nx_offset_upper) 
				     - naxis_data*resamp_fact);
	if (config.resample() && xtest != 0) {
	  nx_offset += 1;
	}
	int ytest = static_cast<int>((naxis2 - ny_offset - ny_offset_upper) 
				     - naxis_data*resamp_fact);
	if (config.resample() && ytest != 0) {
	  ny_offset += 1;
	}
	if (!config.resample()) { 
	  // Use integer or half-integer reference pixel based on
	  // input counts map, even though naxis1 and naxis2 both
	  // must be even.
	  if (naxis_data % 2 == 1) {
	    crpix1 += 0.5;
	    crpix2 += 0.5;
	    nx_offset += 1;
	    ny_offset += 1;
	  }
	}
      }
  

      size_t counter(0);
      std::vector<double>::const_iterator energy = energies.begin();
      Healpix_Map<float> hpm(dataMap.healpixProj()->healpix().Nside(),
			     dataMap.healpixProj()->healpix().Scheme(),
			     SET_NSIDE); 
      
      const double solidAngle = dataMap.solidAngle();
      
      static const double d2tosr = ASTRO_D2R * ASTRO_D2R;
      
      for (int k(0); energy != energies.end(); ++energy, k++) {
	bool interpolate(true);
	bool enforceEnergyRange(false);
	bool computeIntegrals(false);
	formatter.warn() << ".";
	
	WcsMap2 * convolvedMap(0);
	WcsMap2 diffuseMap(diffuseSrc, ref1, ref2,
			   crpix1, crpix2, cdelt1, cdelt2, naxis1, naxis2,
			   *energy, nativeProj, 
			   dataMap.projection().isGalactic(), 
			   interpolate,enforceEnergyRange,computeIntegrals);
		
	if(haveSpatialFunction) {
	  const SpatialFunction* m = 
	    dynamic_cast<const SpatialFunction *>(diffuseSrc.spatialDist());
	  convolvedMap = static_cast<WcsMap2*>(diffuseMap.convolve(*energy, meanpsf, 
								   bexpmap, *m));
	} else {
	  convolvedMap = static_cast<WcsMap2*>(diffuseMap.convolve(*energy, meanpsf, 
								   bexpmap, config.performConvolution()));
	}
	
	// This function assumes that the convolved map is at 
	// a higher resolution than the output map
	// In this case that should be true b/c of the resampling done above
	//fillHealpixFromWcsMap(diffuseMap,*dataMap.healpixProj(),dataMap.pixelIndices(),hpm);
	//fillHealpixFromWcsMap(*convolvedMap,*dataMap.healpixProj(),hpm);
	fillHealpixFromWcsMap(*convolvedMap,*dataMap.healpixProj(),dataMap.pixelIndices(),hpm);
	
	double added(0.0);
	for ( int i(0); i < dataMap.nPixels(); i++, counter++ ) {
	  int iglo = dataMap.localToGlobalIndex(i);
	  modelmap[counter] = hpm[iglo] * solidAngle * d2tosr;
	  added += modelmap[counter];
	}
      }
      // Delete model map for map-based diffuse sources to save memory.  The
      // map will be reloaded dynamically if it is needed again.
      try {
	MapBase * mapBaseObj = 
	  const_cast<MapBase *>(diffuseSrc.mapBaseObject());
	mapBaseObj->deleteMap();
	formatter.info(4) << "SourceMap::makeDiffuseSource: "
			  << "called mapBaseObj->deleteMap()"
			  << std::endl;
      } catch (MapBaseException & eObj) {
	// Not a map-based source, so do nothing.
      }
      return 0;
    }
      
    
    int makePointSourceMap(const Source& src, 
			   const CountsMapBase& dataMap,
			   const PsfIntegConfig& config,
			   const MeanPsf& meanpsf,
			   st_stream::StreamFormatter& formatter,
			   std::vector<float>& modelmap) {
      switch ( dataMap.projection().method() ) {
      case astro::ProjBase::WCS:
	return makePointSourceMap_wcs(src,static_cast<const CountsMap&>(dataMap),config,meanpsf,formatter,modelmap);	 
      case astro::ProjBase::HEALPIX:
	return makePointSourceMap_healpix(src,static_cast<const CountsMapHealpix&>(dataMap),config,meanpsf,formatter,modelmap);
      default:
	break;
      }
      std::string errMsg("SourceMap::makePointSourceMap, did not recognize CountsMapBase type at: ");
      errMsg += dataMap.filename();
      throw std::runtime_error(errMsg);
      return -1;
    }
    
    int makePointSourceMap_wcs(const Source& src, 
			       const CountsMap& dataMap,
			       const PsfIntegConfig& config,
			       const MeanPsf& meanpsf,
			       st_stream::StreamFormatter& formatter,
			       std::vector<float>& modelmap) {
      
      const PointSource& pointSrc = dynamic_cast<const PointSource&>(src);
      const std::vector<Pixel> & pixels(dataMap.pixels());
      std::vector<double> energies;
      dataMap.getEnergies(energies);
      
      long npts = energies.size()*pixels.size();
      modelmap.resize(npts, 0);
      
      const astro::SkyDir & dir(pointSrc.getDir());
      
      const std::vector<double> & exposure = meanpsf.exposure();
      double ref_pixel_size = dataMap.pixelSize();
      
      std::vector< std::vector< double > > pixelOffsets;
      createOffsetMap(src,dataMap,pixelOffsets);
      
      bool galactic(dataMap.isGalactic());
      /// Get the pixel center in pixel coordinates
      double src_lon, src_lat;
      if (galactic) {
	src_lon = dir.l();
	src_lat = dir.b();
      } else {
	src_lon = dir.ra();
	src_lat = dir.dec();
      }
      
      std::pair<double, double> srcCoord = 
	std::pair<double, double>(dataMap.projection().sph2pix(src_lon, 
							       src_lat));
      
      std::vector< std::pair<double, double> > pixCoords(pixels.size());
      std::vector<Pixel>::const_iterator pixel(pixels.begin());
      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
	
	double lon, lat;
	if (galactic) {
	  lon = pixel->dir().l();
	  lat = pixel->dir().b();
	} else {
	  lon = pixel->dir().ra();
	  lat = pixel->dir().dec();
	}
	pixCoords[j] = std::pair<double, double>(pixel->proj().sph2pix(lon, lat));
      }
      
      if (config.performConvolution()) {
	long icount(0);
	std::vector<double> mapCorrections(energies.size(), 1.);
	std::vector<double> mapIntegrals(energies.size());
	bool apply_map_corrections = config.applyPsfCorrections() &&
	  dataMap.withinBounds(dir, energies.at(energies.size()/2), 4);
	  
	//std::vector<Pixel>::const_iterator pixel(pixels.begin());
	pixel = pixels.begin();
	for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
	  std::vector<double>::const_iterator energy = energies.begin();
	  for (int k = 0; energy != energies.end(); ++energy, k++) {
            unsigned long indx = k*pixels.size() + j;
            if (config.verbose() && (icount % (npts/20)) == 0) {
	      formatter.warn() << ".";
            }
            double offset(dir.difference(pixel->dir())*180./M_PI);
            double psf_value(psfValueEstimate(meanpsf, energies.at(k), 
                                              dir, *pixel, srcCoord,
					      pixCoords[j],pixelOffsets,ref_pixel_size,config));
	    // This is the value without the exposure, which we need for the map integrals
	    double value = psf_value*pixel->solidAngle();
	    mapIntegrals[k] += value;
	    // Now we factor in the exposure
 	    value *= exposure.at(k);
	    modelmap.at(indx) += value;
            icount++;
	  }
	}
	// We have finished the loop of the pixels, now we apply the map corrections, if requested
	if ( apply_map_corrections ) {
	  double psfRadius = maxPsfRadius(pointSrc,dataMap);	  
	  std::vector<double>::const_iterator energy = energies.begin();
	  for ( int kk = 0; energy != energies.end(); ++energy, kk++) {
	    // This is the correction factor for this energy layer
	    double correction_factor = meanpsf.integral(psfRadius, energies[kk] ) / mapIntegrals[kk];
	    // Loop over the pixels and apply the correction factor
	    std::vector<float>::iterator itr_corr_b =  modelmap.begin() + kk*pixels.size();
	    std::vector<float>::iterator itr_corr_e =  itr_corr_b + pixels.size();
	    for ( std::vector<float>::iterator itr_corr = itr_corr_b; itr_corr != itr_corr_e; itr_corr++) {
	      *itr_corr *= correction_factor;
	    }
	  }
	}
      } else {
	const std::vector<Pixel>::const_iterator targetPixel = 
	  Pixel::find(pixels.begin(), pixels.end(),
		      Pixel(dir.ra(), dir.dec(), 1), 2.);
	if (targetPixel != pixels.end()) {
	  size_t ipix = targetPixel - pixels.begin();
	  std::vector<double>::const_iterator energy = energies.begin();
	  for (int k = 0; energy != energies.end(); ++energy, k++) {
            size_t indx = k*pixels.size() + ipix;
            modelmap.at(indx) = exposure.at(k);
	  }
	}
      }
      return 0;
    }
    
    int makePointSourceMap_healpix(const Source& src, 
				   const CountsMapHealpix& dataMap,
				   const PsfIntegConfig& config,
				   const MeanPsf& meanpsf,
				   st_stream::StreamFormatter& formatter,
				   std::vector<float>& modelmap) {

      const PointSource& pointSrc = dynamic_cast<const PointSource &>(src);
      const std::vector<Pixel> & pixels(dataMap.pixels());
      std::vector<double> energies;
      dataMap.getEnergies(energies);
      
      long npts = energies.size()*pixels.size();
      modelmap.resize(npts, 0);
      
      const astro::SkyDir & dir(pointSrc.getDir());
      const std::vector<double> & exposure = meanpsf.exposure();
      
      int nPix = dataMap.nPixels();
      
      if (config.performConvolution()) {    
	size_t indx(0);
	const Healpix_Base& hp = dataMap.healpixProj()->healpix();
	Healpix_Map<float> hpmap(hp.Nside(),hp.Scheme(),SET_NSIDE);
	for (int k(0); k < energies.size(); k++ ) {
	  formatter.warn() << ".";
	  hpmap.fill(0.);
	  ConvolveHealpix::fillMapWithPSF_refDir(meanpsf,energies[k],dir,dataMap.isGalactic(),hpmap);
	  for (int i(0); i < nPix; i++, indx++ ) {
	    modelmap.at(indx) = exposure.at(k) * hpmap[ dataMap.localToGlobalIndex(i) ];
	  }
	}      
      } else {
	dataMap.projection();
	double ipix(0);
	double jpix(0);
	st_facilities::Util::skyDir2pixel(dataMap.projection(),dir,ipix,jpix);
	int iglo = dataMap.globalToLocalIndex(ipix);
	// Pixel is not-defined, return with error code
	if ( iglo < 0 ) return -1;
	if ( ipix >= 0 ) {
	  std::vector<double>::const_iterator energy = energies.begin();
	  size_t indx = size_t(ipix);
	  for (int k = 0; energy != energies.end(); ++energy, ++k ) {
	    size_t indx = k*pixels.size() + iglo;
	    modelmap.at(indx) = dataMap.solidAngle() * exposure.at(k);
	  }      
	}
      }   
      return 0;
    }

    

    double psfValueEstimate(const MeanPsf & meanPsf, double energy, 
			    const astro::SkyDir & srcDir, 
			    const Pixel & pixel,
			    const std::pair<double, double>& srcCoord,
			    const std::pair<double, double>& pixCoord,
			    const std::vector< std::vector< double > >& pixelOffsets,
			    double ref_pixel_size,
			    const PsfIntegConfig& config) {
   
      switch ( config.integ_type() ) {
      case PsfIntegConfig::adaptive:
	return integrate_psf_adaptive(meanPsf, energy, srcDir, pixel,
				      srcCoord, pixCoord, pixelOffsets, ref_pixel_size, config);
      case PsfIntegConfig::pixel_center:
	return psfValue_pixelCenter(meanPsf, energy, srcDir, pixel);
      case PsfIntegConfig::annular:
	return psfValue_annular(meanPsf, energy, srcDir, pixel);
      default:
	break;
      }
      throw std::runtime_error(" PSFUtils::psfValueEstimate: Unknown PSF estimator: ");
    }

    double psfValue_pixelCenter(const MeanPsf & meanPsf, double energy, 
				const astro::SkyDir & srcDir, 
				const Pixel & pixel) {
      double offset(srcDir.difference(pixel.dir())*180./M_PI);
      double pixelSolidAngle(pixel.solidAngle());
      // FIXME, should there be a solid angle correction here?
      return meanPsf(energy, offset);
    }

    double psfValue_annular(const MeanPsf & meanPsf, double energy, 
			    const astro::SkyDir & srcDir, 
			    const Pixel & pixel) {
      
      double offset(srcDir.difference(pixel.dir())*180./M_PI);
      double pixelSolidAngle(pixel.solidAngle());
      
      /// To estimate the psf value averaged over a pixel, average the psf
      /// over an annulus centered on the source position with approximately
      /// the same extent in theta as the pixel in question.
      static double sqrt2(std::sqrt(2.));
      double pixel_value(0);
      double pixel_size(std::sqrt(pixelSolidAngle)*180./M_PI);
      if (pixel_size/2. >= offset) {
	/// Average over an acceptance cone with the same solid angle
	/// as the central pixel.
	double radius(std::acos(1. - pixelSolidAngle/2./M_PI)*180./M_PI);
	pixel_value = meanPsf.integral(radius, energy)/pixelSolidAngle;
	// if (offset < pixel_size/sqrt2) {
	//    pixel_value = integrate_psf(meanPsf, energy, srcDir, pixel);
      } else {
	// Use integral over annulus with pixel_size width centered on
	// the offset angle to estimate average psf value within a pixel
	// at the offset.
	double theta1(offset - pixel_size/2.);
	double theta2(offset + pixel_size/2.);
	pixel_value = ( (meanPsf.integral(theta2, energy)
			 - meanPsf.integral(theta1, energy))
			/(2.*M_PI*(std::cos(theta1*M_PI/180.)
				   - std::cos(theta2*M_PI/180.))) );
      }
      return pixel_value;
    }

    double integrate_psf_adaptive(const MeanPsf & meanPsf, 
				  double energy, 
				  const astro::SkyDir & srcDir, 
				  const Pixel & pixel, 
				  const std::pair<double, double>& srcCoord,
				  const std::pair<double, double>& pixCoord,
				  const std::vector< std::vector< double > >& pixelOffsets,
				  double ref_pixel_size,
				  const PsfIntegConfig& config) {
      
      double offset(srcDir.difference(pixel.dir())*180./M_PI);
      const int max_npts = 64;
      double ftol_threshold = config.psfEstimatorFtol();
      double peak_threshold = config.psfEstimatorPeakTh();
      double peak_val = meanPsf.peakValue(energy);
      double v0 = meanPsf(energy, offset);      
      double v1 = 0;
      double ferr = 0;
      double peak_ratio = v0/peak_val;
      
      if(peak_ratio < peak_threshold)
	return v0;
      
      int npts = 1;
      
      while(1) {
	npts *= 2;	
	v1 = integrate_psf_fast(meanPsf, energy, srcDir, pixel, srcCoord, pixCoord, pixelOffsets, ref_pixel_size, npts);
	
	if(v1-v0 == 0 || v1 == 0)
	  break;
	
	ferr = std::fabs((v1-v0)/v1);
	if(ferr < ftol_threshold || npts >= max_npts)
	  break;
	
	v0=v1;
      }
      
      return v1;
    }

    double integrate_psf_fast(const MeanPsf & meanPsf, double energy, 
			      const astro::SkyDir & srcDir, 
			      const Pixel & pixel, 
			      const std::pair<double, double>& srcCoord,
			      const std::pair<double, double>& pixCoord, 
			      const std::vector< std::vector< double > >& pixelOffsets,
			      double ref_pixel_size,
			      size_t npts)  {
      
      double psf_value(0);
      double dstep(1./static_cast<double>(npts));
      for (size_t i(0); i < npts; i++) {
	double x(pixCoord.first + i*dstep - 0.5 + 0.5*dstep);
	for (size_t j(0); j < npts; j++) {
	  double y(pixCoord.second + j*dstep - 0.5 + 0.5*dstep);
	  double pix_offset = 
	    ref_pixel_size*std::sqrt(std::pow(srcCoord.first-x,2) +  
				 std::pow(srcCoord.second-y,2));
	  // We have to subtract off 1 to get from the FITS convention to the 
	  // array indices
	  double scale = bilinear_on_grid(y-1,x-1,pixelOffsets);	  
	  psf_value += meanPsf(energy, pix_offset*(1+scale));
	}
      }
      return psf_value/static_cast<double>(npts*npts);
    }

    double integrate_psf(const MeanPsf & meanPsf, double energy,
			 const astro::SkyDir & srcDir, const Pixel & pixel, size_t npts) {

      bool galactic(pixel.proj().isGalactic());
      
      /// Get the pixel center in pixel coordinates
      double lon, lat;
      if (galactic) {
	lon = pixel.dir().l();
	lat = pixel.dir().b();
      } else {
	lon = pixel.dir().ra();
	lat = pixel.dir().dec();
      }
      std::pair<double, double> pix_coords(pixel.proj().sph2pix(lon, lat));
      
      // loop over subpixels in longitudinal and latitudinal directions
      double psf_value(0);
      double dstep(1./static_cast<double>(npts));
      for (size_t i(0); i < npts; i++) {
	double x(pix_coords.first + i*dstep - 0.5 + 0.5*dstep);
	for (size_t j(0); j < npts; j++) {
	  double y(pix_coords.second + j*dstep - 0.5 + 0.5*dstep);
	  std::pair<double, double> pix_dir(pixel.proj().pix2sph(x, y));
	  double offset(0);
	  if (galactic) {
            astro::SkyDir my_dir(pix_dir.first, pix_dir.second,
                                 astro::SkyDir::GALACTIC);
            offset = my_dir.difference(srcDir)*180./M_PI;
	  } else {
            astro::SkyDir my_dir(pix_dir.first, pix_dir.second,
                                 astro::SkyDir::EQUATORIAL);
            offset = my_dir.difference(srcDir)*180./M_PI;
	  }
	  psf_value += meanPsf(energy, offset);
	}
      }
      return psf_value/static_cast<double>(npts*npts);
    }

  } // namespace PSFUtils
 
} // namespace Likelihood
