/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/SourceMap.cxx,v 1.114 2016/03/03 21:51:38 mdwood Exp $
 */

#include <cmath>

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "astro/SkyProj.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMapBase.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/CountsMapHealpix.h"
#include "Likelihood/ConvolveHealpix.h"
#include "Likelihood/AppHelpers.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/MapBase.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/Source.h"
#define ST_DLL_EXPORTS
#include "Likelihood/SourceMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/SpatialFunction.h"

#include "Likelihood/WcsMap2.h"
#include "Likelihood/HealpixProjMap.h"

namespace {
   double my_acos(double mu) {
      if (mu > 1) {
         return 0;
      } else if (mu < -1) {
         return M_PI;
      } else {
         return acos(mu);
      }
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
}

namespace Likelihood {

  
void SourceMap::fillHealpixFromWcsMap(const WcsMap2& inputMap,
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


void SourceMap::fillHealpixFromWcsMap(const WcsMap2& inputMap,
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


SourceMap::SourceMap(Source * src, const CountsMapBase * dataMap,
                     const Observation & observation, 
                     bool applyPsfCorrections,
                     bool performConvolution,
                     bool resample,
                     double resamp_factor,
                     double minbinsz,
                     bool verbose,
		     const SourceMap* weights)
   : m_name(src->getName()),
     m_srcType(src->getType()),
     m_dataMap(dataMap),
     m_observation(observation),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_weights(weights),
     m_deleteDataMap(false),
     m_pixelOffset(),
     m_pixelCoordX(),
     m_pixelCoordY(),
     m_psfEstimatorMethod("adaptive"),
     m_psfEstimatorFtol(1E-3),
     m_psfEstimatorPeakTh(1E-6) {
   if (verbose) {
      m_formatter->warn() << "Generating SourceMap for " << m_name;
   }

   bool havePointSource = dynamic_cast<PointSource *>(src) != 0;
   bool haveDiffuseSource = dynamic_cast<DiffuseSource *>(src) != 0;

   if (haveDiffuseSource) {
      makeDiffuseMap(src, dataMap, applyPsfCorrections,
                     performConvolution, resample, resamp_factor,
                     minbinsz, verbose);
   } else if (havePointSource) {
      makePointSourceMap(src, dataMap, applyPsfCorrections,
                         performConvolution, verbose);
   }
   if (verbose) {
      m_formatter->warn() << "!" << std::endl;
   }
   applyPhasedExposureMap();
   computeNpredArray();
}

SourceMap::SourceMap(const ProjMap& weight_map,
		     const CountsMapBase * dataMap,
		     const Observation & observation,
		     bool verbose)
   : m_name("__weights__"),
     m_srcType("Diffuse"),
     m_dataMap(dataMap),
     m_observation(observation),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_weights(0),
     m_deleteDataMap(false),
     m_pixelOffset(),
     m_pixelCoordX(),
     m_pixelCoordY(),
     m_psfEstimatorMethod("adaptive"),
     m_psfEstimatorFtol(1E-3),
     m_psfEstimatorPeakTh(1E-6) {
   if (verbose) {
      m_formatter->warn() << "Generating SourceMap for " << m_name;
   }

   makeProjectedMap(weight_map);
   // This is just to make sure that they are there.  
   // For this type of map they shouldn't be used for anything
   computeNpredArray();
   if (verbose) {
      m_formatter->warn() << "!" << std::endl;
   }
 }


SourceMap::SourceMap(const std::string & sourceMapsFile,
                     const std::string & srcName,
                     const Observation & observation,
		     const SourceMap* weights) 
   : m_name(srcName),
     m_dataMap(AppHelpers::readCountsMap(sourceMapsFile)),
     m_observation(observation),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_weights(weights),
     m_deleteDataMap(true),
     m_pixelOffset(),
     m_pixelCoordX(),
     m_pixelCoordY(),
     m_psfEstimatorMethod("adaptive"),
     m_psfEstimatorFtol(1E-3),
     m_psfEstimatorPeakTh(1E-6) {

   m_model.clear();
   bool ok(false);
   switch ( m_dataMap->projection().method()  ) {
   case astro::ProjBase::WCS:
     readImage(sourceMapsFile);
     ok = true;
     break;
   case astro::ProjBase::HEALPIX:
     readTable_healpix(sourceMapsFile);
     ok = true;
     break;
   default:
     break;
   }

   if ( !ok ) {
     std::string errMsg("SourceMap did not recognize CountsMapBase type at: ");
     errMsg += sourceMapsFile;
     throw std::runtime_error(errMsg);

   }
   applyPhasedExposureMap();
   computeNpredArray();
}


void SourceMap::makeDiffuseMap(Source * src,
			       const CountsMapBase * dataMap,
			       bool applyPsfCorrections,
			       bool performConvolution,
			       bool resample,
			       double resamp_factor,
			       double minbinsz,
			       bool verbose) {
  const CountsMapHealpix* hpxmap(0);
  switch ( dataMap->projection().method() ) {
  case astro::ProjBase::WCS:
    makeDiffuseMap_wcs(src,static_cast<const CountsMap*>(dataMap),applyPsfCorrections,
		       performConvolution,resample,resamp_factor,
		       minbinsz,verbose);
    return;
  case astro::ProjBase::HEALPIX:
    hpxmap = static_cast<const CountsMapHealpix*>(dataMap);
    if ( hpxmap->mapRadius() > 45 ) {
      // The ROI takes up a large fraction of the sky, 
      // let's do the convolution using HEALPix
      makeDiffuseMap_healpix(src,hpxmap,applyPsfCorrections,
			     performConvolution,resample,resamp_factor,
			     minbinsz,verbose);
    } else {
      // The ROI only take up a relatively small fraction of the sky (< 15%)
      // we will do the convlution in the native projection, then convert to healpix
      makeDiffuseMap_native(src,hpxmap,applyPsfCorrections,
			    performConvolution,resample,resamp_factor,
			    minbinsz,verbose);      
    }
    return;
  default:
    break;
  }
  std::string errMsg("Unrecognized projection type for map object for source: ");
  errMsg += src->getName();
  throw std::runtime_error(errMsg);
  return;
}



void SourceMap::makeDiffuseMap_native(Source * src,
				      const CountsMapHealpix * dataMap,
				      bool applyPsfCorrections,
				      bool performConvolution,
				      bool resample,
				      double resamp_factor,
				      double minbinsz,
				      bool verbose) {
  DiffuseSource * diffuseSrc = dynamic_cast<DiffuseSource *>(src);

  bool haveSpatialFunction = dynamic_cast<const SpatialFunction *>(diffuseSrc->spatialDist()) != 0;
  // If the diffuse source is represented by an underlying map, then
  // rebin according to the minimum bin size.
  std::string nativeProj("STG");
  try {
    MapBase & tmp(*diffuseSrc->mapBaseObject());
    nativeProj = tmp.projmap().getProj()->projType();
    double pixSize = tmp.projmap().pixelSize();
    if ( pixSize < minbinsz ) {
      unsigned int factor = static_cast<unsigned int>(minbinsz/pixSize);
      m_formatter->info(4) << "\nrebinning factor: " 
			   << factor << std::endl;
      if (factor > 1) {
	tmp.rebin(factor);
      }
    }
  } catch (MapBaseException &) {
    // do nothing
  }

  const std::vector<Pixel> & pixels(dataMap->pixels());
  std::vector<double> energies;
  dataMap->getEnergies(energies);

  long npts = energies.size()*pixels.size();
  m_model.resize(npts, 0);

  std::vector<Pixel>::const_iterator pixel = pixels.begin();

  const astro::SkyDir & mapRefDir = dataMap->refDir();
  if (!resample) {
    resamp_factor = 1;
  } else {
    resamp_factor = std::max(resamp_factor, 
			     computeResampFactor(*diffuseSrc, *dataMap));
  }
  m_formatter->info(4) << "\nresampling factor: " 
		       << resamp_factor << std::endl;

  double crpix1, crpix2;
  int naxis1(0), naxis2(0);
  double cdelt1 = -1.*dataMap->pixelSize()/resamp_factor;
  double cdelt2 = dataMap->pixelSize()/resamp_factor;
  size_t nx_offset(0), ny_offset(0);
  size_t nx_offset_upper(0), ny_offset_upper(0);
  double ref1 = dataMap->isGalactic() ? mapRefDir.l() : mapRefDir.ra();
  double ref2 = dataMap->isGalactic() ? mapRefDir.b() : mapRefDir.dec();


  // For spatial functions the map only needs to be big enough to
  // encompass the ROI
  if (haveSpatialFunction) {
    double radius = std::min(180., ::maxRadius(pixels, mapRefDir) + std::fabs(cdelt1) );
    int mapsize(2*static_cast<int>(radius/std::fabs(cdelt1)));
    naxis1 = mapsize;
    naxis2 = mapsize;
    crpix1 = (mapsize+1.0)*0.5;
    crpix2 = (mapsize+1.0)*0.5;
  } else {

    double radius = std::min(180., ::maxRadius(pixels, mapRefDir) + 10.);
    // Conforming maps have abs(CDELT1) == abs(CDELT2).  This
    // expression for the mapsize ensures that the number of
    // pixels in each dimension is even.
    int mapsize(2*static_cast<int>(radius/std::fabs(cdelt1)));

    m_formatter->info(4) << "mapsize: " << mapsize << std::endl;
    naxis1 = mapsize;
    naxis2 = mapsize;
    crpix1 = (naxis1 + 1.)/2.;
    crpix2 = (naxis2 + 1.)/2.;
  
    int naxis_data = int( dataMap->mapRadius() / dataMap->pixelSize() );

    nx_offset = 
      static_cast<size_t>((mapsize - naxis_data*resamp_factor)/2);
    ny_offset = 
      static_cast<size_t>((mapsize - naxis_data*resamp_factor)/2);
    nx_offset_upper = 
      static_cast<size_t>((mapsize - naxis_data*resamp_factor)/2);
    ny_offset_upper = 
      static_cast<size_t>((mapsize - naxis_data*resamp_factor)/2);
    /// For cases where the resampling factor is an odd number, 
    /// there may be a row or column of pixels not accounted for
    /// by n[xy]_offset.  Here we add that row or column back in if
    /// it is missing.
    int xtest = static_cast<int>((naxis1 - nx_offset - nx_offset_upper) 
				 - naxis_data*resamp_factor);
    if (resample && xtest != 0) {
      nx_offset += 1;
    }
    int ytest = static_cast<int>((naxis2 - ny_offset - ny_offset_upper) 
				 - naxis_data*resamp_factor);
    if (resample && ytest != 0) {
      ny_offset += 1;
    }
    if (!resample) { 
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
  Healpix_Map<float> hpm(dataMap->healpixProj()->healpix().Nside(),
			 dataMap->healpixProj()->healpix().Scheme(),
			 SET_NSIDE); 

  const double solidAngle = dataMap->solidAngle();
  
  static const double d2tosr = ASTRO_D2R * ASTRO_D2R;

  for (int k(0); energy != energies.end(); ++energy, k++) {
    bool interpolate(true);
    bool enforceEnergyRange(false);
    bool computeIntegrals(false);
    m_formatter->warn() << ".";

    WcsMap2 * convolvedMap(0);
    WcsMap2 diffuseMap(*diffuseSrc, ref1, ref2,
		       crpix1, crpix2, cdelt1, cdelt2, naxis1, naxis2,
		       *energy, nativeProj, 
		       dataMap->projection().isGalactic(), 
		       interpolate,enforceEnergyRange,computeIntegrals);

    const MeanPsf & meanpsf(m_observation.meanpsf());
    const BinnedExposureBase & bexpmap = m_observation.bexpmap();

    if(haveSpatialFunction) {
      const SpatialFunction* m = 
	dynamic_cast<const SpatialFunction *>(diffuseSrc->spatialDist());
      convolvedMap = static_cast<WcsMap2*>(diffuseMap.convolve(*energy, meanpsf, 
								 bexpmap, *m));
    } else {
      convolvedMap = static_cast<WcsMap2*>(diffuseMap.convolve(*energy, meanpsf, 
							       bexpmap, performConvolution));
    }

    // This function assumes that the convolved map is at 
    // a higher resolution than the output map
    // In this case that should be true b/c of the resampling done above
    //fillHealpixFromWcsMap(diffuseMap,*dataMap->healpixProj(),dataMap->pixelIndices(),hpm);
    //fillHealpixFromWcsMap(*convolvedMap,*dataMap->healpixProj(),hpm);
    fillHealpixFromWcsMap(*convolvedMap,*dataMap->healpixProj(),dataMap->pixelIndices(),hpm);

    double added(0.0);
    for ( int i(0); i < dataMap->nPixels(); i++, counter++ ) {
      int iglo = dataMap->localToGlobalIndex(i);
      m_model[counter] = hpm[iglo] * solidAngle * d2tosr;
      added += m_model[counter];
    }
  }
  // Delete model map for map-based diffuse sources to save memory.  The
  // map will be reloaded dynamically if it is needed again.
  try {
    MapBase * mapBaseObj = 
      const_cast<MapBase *>(diffuseSrc->mapBaseObject());
    mapBaseObj->deleteMap();
    m_formatter->info(4) << "SourceMap::makeDiffuseSource: "
			 << "called mapBaseObj->deleteMap()"
			 << std::endl;
  } catch (MapBaseException & eObj) {
    // Not a map-based source, so do nothing.
  }
}



void SourceMap::makeDiffuseMap_healpix(Source * src,
				       const CountsMapHealpix * dataMap,
				       bool applyPsfCorrections,
				       bool performConvolution,
				       bool resample,
				       double resamp_factor,
				       double minbinsz,
				       bool verbose) {

  DiffuseSource * diffuseSrc = dynamic_cast<DiffuseSource *>(src);
  // If the diffuse source is represented by an underlying map, then
  // rebin according to the minimum bin size.
  try {
    MapBase & tmp(*diffuseSrc->mapBaseObject());
    double pixSize = tmp.projmap().pixelSize();
    if ( pixSize < minbinsz ) {   
      unsigned int factor = static_cast<unsigned int>(minbinsz/pixSize);
      m_formatter->info(4) << "\nrebinning factor: " 
			   << factor << std::endl;
      if (factor > 1) {
	tmp.rebin(factor);
      }
    }
  } catch (MapBaseException &) {
    // do nothing
  }
  
  double solidAngle = dataMap->solidAngle();
  std::vector<double> energies;
  dataMap->getEnergies(energies);
  
  long npts = energies.size()*dataMap->nPixels();
  m_model.resize(npts, 0);
  const astro::SkyDir & mapRefDir = dataMap->refDir();
  if (!resample) {
    resamp_factor = 1;
  } else {
    resamp_factor = std::max(resamp_factor, 
			     computeResampFactor(*diffuseSrc, *dataMap));
  }
  // We need to move the resampling factor to the nearest power of two.
  int resamp_test(1);
  while ( resamp_test < resamp_factor ) {
    resamp_test *= 2;
  }
  resamp_factor = resamp_test;

  Healpix_Ordering_Scheme scheme = dataMap->healpixProj()->healpix().Scheme();
  int nside_orig = dataMap->healpixProj()->healpix().Nside(); 
  int resamp_nside = resamp_factor*nside_orig;
  size_t counter(0);
  static const double ALLSKY_RADIUS(180.);
  bool interpolate(false);
  std::vector<double>::const_iterator energy = energies.begin();
  double scaleFactor = solidAngle/(resamp_factor*resamp_factor);
  // This is the index in the output vector, it does _not_ get reset between energy layers
  int outidx(0);
  for (int k(0); energy != energies.end(); ++energy, k++) {
    m_formatter->warn() << ".";
    HealpixProjMap diffuseMap(*diffuseSrc, resamp_nside,
			      scheme,SET_NSIDE,
			      *energy,dataMap->projection().isGalactic(),
			      ALLSKY_RADIUS,mapRefDir.ra(), mapRefDir.dec(),
			      interpolate);
    const MeanPsf & meanpsf = m_observation.meanpsf();
    const BinnedExposureBase & bexpmap = m_observation.bexpmap();
    ProjMap* cmap = diffuseMap.convolve(*energy,meanpsf,bexpmap,performConvolution);
    HealpixProjMap* convolvedMap = static_cast<HealpixProjMap*>(cmap);
    Healpix_Map<float> outmap(nside_orig,scheme,SET_NSIDE);
    outmap.Import_degrade(convolvedMap->image()[0]);
    if ( nside_orig == resamp_nside ) {
      outmap = convolvedMap->image()[0];
    } else {
      outmap.Import_degrade(convolvedMap->image()[0]);
    }
    for (size_t i(0); i < dataMap->nPixels(); i++,outidx++) {
      int glo = dataMap->localToGlobalIndex(i);
      m_model[outidx] = outmap[glo];
    }
  }
  try {
    MapBase * mapBaseObj = 
      const_cast<MapBase *>(diffuseSrc->mapBaseObject());
    mapBaseObj->deleteMap();
    m_formatter->info(4) << "SourceMap::makeDiffuseSource: "
			 << "called mapBaseObj->deleteMap()"
			 << std::endl;
  } catch (MapBaseException & eObj) {
    // Not a map-based source, so do nothing.
  }
}

void SourceMap::makeDiffuseMap_wcs(Source * src, 
				   const CountsMap* dataMap,
				   bool applyPsfCorrections,
				   bool performConvolution,
				   bool resample,
				   double resamp_factor,
				   double minbinsz,
				   bool verbose) {
   DiffuseSource * diffuseSrc = dynamic_cast<DiffuseSource *>(src);
   
   bool haveSpatialFunction = dynamic_cast<const SpatialFunction *>(diffuseSrc->spatialDist()) != 0;

// If the diffuse source is represented by an underlying map, then
// rebin according to the minimum bin size.
   try {
      MapBase & tmp(*diffuseSrc->mapBaseObject());
      double pixSize = tmp.projmap().pixelSize();
      if ( pixSize < minbinsz ) {
	 unsigned int factor = static_cast<unsigned int>(minbinsz/pixSize);
         m_formatter->info(4) << "\nrebinning factor: " 
                              << factor << std::endl;
         if (factor > 1) {
            tmp.rebin(factor);
         }
      }
   } catch (MapBaseException &) {
      // do nothing
   }

   const std::vector<Pixel> & pixels(dataMap->pixels());
   std::vector<double> energies;
   dataMap->getEnergies(energies);

   long npts = energies.size()*pixels.size();
   m_model.resize(npts, 0);

   std::vector<Pixel>::const_iterator pixel = pixels.begin();

   const astro::SkyDir & mapRefDir = dataMap->refDir();
   if (!resample) {
      resamp_factor = 1;
   } else {
      resamp_factor = std::max(resamp_factor, 
                               computeResampFactor(*diffuseSrc, *dataMap));
   }
   m_formatter->info(4) << "\nresampling factor: " 
                        << resamp_factor << std::endl;
   double crpix1, crpix2;
   int naxis1, naxis2;
   double cdelt1 = dataMap->cdelt1()/resamp_factor;
   double cdelt2 = dataMap->cdelt2()/resamp_factor;
   size_t nx_offset(0), ny_offset(0);
   size_t nx_offset_upper(0), ny_offset_upper(0);
   if (haveSpatialFunction) {
     naxis1 = static_cast<int>(dataMap->naxis1()*resamp_factor);
     naxis2 = static_cast<int>(dataMap->naxis2()*resamp_factor);      
     crpix1 = resamp_factor*(dataMap->crpix1() - 0.5)+0.5;
     crpix2 = resamp_factor*(dataMap->crpix2() - 0.5)+0.5;  
   } else if (dataMap->conformingMap()) {
      double radius = std::min(180., ::maxRadius(pixels, mapRefDir) + 10.);
      // Conforming maps have abs(CDELT1) == abs(CDELT2).  This
      // expression for the mapsize ensures that the number of
      // pixels in each dimension is even.
      int mapsize(2*static_cast<int>(radius/std::fabs(cdelt1)));
      m_formatter->info(4) << "mapsize: " << mapsize << std::endl;
      naxis1 = mapsize;
      naxis2 = mapsize;
      crpix1 = (naxis1 + 1.)/2.;
      crpix2 = (naxis2 + 1.)/2.;
      nx_offset = 
         static_cast<size_t>((mapsize - dataMap->naxis1()*resamp_factor)/2);
      ny_offset = 
         static_cast<size_t>((mapsize - dataMap->naxis2()*resamp_factor)/2);
      nx_offset_upper = 
         static_cast<size_t>((mapsize - dataMap->naxis1()*resamp_factor)/2);
      ny_offset_upper = 
         static_cast<size_t>((mapsize - dataMap->naxis2()*resamp_factor)/2);
      /// For cases where the resampling factor is an odd number, 
      /// there may be a row or column of pixels not accounted for
      /// by n[xy]_offset.  Here we add that row or column back in if
      /// it is missing.
      int xtest = static_cast<int>((naxis1 - nx_offset - nx_offset_upper) 
                                   - dataMap->naxis1()*resamp_factor);
      if (resample && xtest != 0) {
         nx_offset += 1;
      }
      int ytest = static_cast<int>((naxis2 - ny_offset - ny_offset_upper) 
                                   - dataMap->naxis2()*resamp_factor);
      if (resample && ytest != 0) {
         ny_offset += 1;
      }
      if (!resample) { 
         // Use integer or half-integer reference pixel based on
         // input counts map, even though naxis1 and naxis2 both
         // must be even.
         if (dataMap->naxis1() % 2 == 1) {
            crpix1 += 0.5;
            nx_offset += 1;
         }
         if (dataMap->naxis2() % 2 == 1) {
            crpix2 += 0.5;
            ny_offset += 1;
         }
      }
   } else {
      // The counts map was not created by gtbin, so just adopt the
      // map geometry without adding padding for psf leakage since
      // this cannot be done in general without redefining the
      // reference pixel and reference direction.
      naxis1 = static_cast<int>(dataMap->naxis1()*resamp_factor);
      naxis2 = static_cast<int>(dataMap->naxis2()*resamp_factor);
      // Ensure an even number of pixels in each direction.
      if (naxis1 % 2 == 1) {
         naxis1 += 1;
      }
      if (naxis2 % 2 == 1) {
         naxis2 += 1;
      }
      crpix1 = dataMap->crpix1()*resamp_factor;
      crpix2 = dataMap->crpix2()*resamp_factor;
      nx_offset_upper += 1;
      ny_offset_upper += 1;
   }

   const MeanPsf & meanpsf(m_observation.meanpsf());   
   const BinnedExposureBase & bexpmap(m_observation.bexpmap());

   size_t counter(0);
   std::vector<double>::const_iterator energy = energies.begin();
   for (int k(0); energy != energies.end(); ++energy, k++) {
      bool interpolate(true);

      WcsMap2 * convolvedMap(0);      
      WcsMap2 diffuseMap(*diffuseSrc,  
			 (dataMap->projection().isGalactic() ? mapRefDir.l() : mapRefDir.ra()), 
			 (dataMap->projection().isGalactic() ? mapRefDir.b() : mapRefDir.dec()),
                         crpix1, crpix2, cdelt1, cdelt2, naxis1, naxis2,
                         *energy, dataMap->proj_name(), 
                         dataMap->projection().isGalactic(), 
                         interpolate);

      if(haveSpatialFunction) {
	const SpatialFunction* m = 
	  dynamic_cast<const SpatialFunction *>(diffuseSrc->spatialDist());
	convolvedMap = static_cast<WcsMap2*>(diffuseMap.convolve(*energy, meanpsf, 
								 bexpmap, *m));
      } else {
	convolvedMap = static_cast<WcsMap2*>(diffuseMap.convolve(*energy, meanpsf, 
								 bexpmap, performConvolution));
      }

      size_t rfac(static_cast<size_t>(resamp_factor));
      double solid_angle;
      double added(0.0);
      for (size_t j(ny_offset); j < naxis2 - ny_offset_upper; j++) {
         for (size_t i(nx_offset); i < naxis1 - nx_offset_upper; i++) {
            if ((i % rfac == 0) && (j % rfac == 0)) {
               counter++;
               if (verbose && (counter % (npts/20)) == 0) {
                  m_formatter->warn() << ".";
               }
            }
            size_t pix_index = ((j-ny_offset)/rfac)*dataMap->naxis1() 
               + ((i-nx_offset)/rfac);
            solid_angle = pixels.at(pix_index).solidAngle();
            size_t indx = k*dataMap->naxis1()*dataMap->naxis2() + pix_index;
            m_model[indx] += (convolvedMap->image()[0][j][i]
                              /resamp_factor/resamp_factor
                              *solid_angle);
	    added += m_model[indx];
         }
      }
      delete convolvedMap;
   }
//   computeNpredArray();
// Delete model map for map-based diffuse sources to save memory.  The
// map will be reloaded dynamically if it is needed again.
   try {
      MapBase * mapBaseObj = 
         const_cast<MapBase *>(diffuseSrc->mapBaseObject());
      mapBaseObj->deleteMap();
      m_formatter->info(4) << "SourceMap::makeDiffuseSource: "
                           << "called mapBaseObj->deleteMap()"
                           << std::endl;
   } catch (MapBaseException & eObj) {
      // Not a map-based source, so do nothing.
   }
}

void SourceMap::makePointSourceMap(Source * src,
				   const CountsMapBase * dataMap,
				   bool applyPsfCorrections,
				   bool performConvolution,
				   bool verbose) {
   bool ok = false;
   switch ( dataMap->projection().method() ) {
   case astro::ProjBase::WCS:
     makePointSourceMap_wcs(src,static_cast<const CountsMap*>(dataMap),
			    applyPsfCorrections,performConvolution,verbose);
     ok = true;
     break;
   case astro::ProjBase::HEALPIX:
     makePointSourceMap_healpix(src,static_cast<const CountsMapHealpix*>(dataMap),
				applyPsfCorrections,performConvolution,verbose);
     ok = true;
     break;
   default:
     break;
   }
   if ( !ok ) {
     std::string errMsg("SourceMap::makePointSourceMap, did not recognize CountsMapBase type at: ");
     errMsg += dataMap->filename();
     throw std::runtime_error(errMsg);
   }
   return;
}

void SourceMap::makePointSourceMap_wcs(Source * src,
				       const CountsMap * dataMap,
				       bool applyPsfCorrections,
				       bool performConvolution,
				       bool verbose) {

   PointSource * pointSrc = dynamic_cast<PointSource *>(src);
   const std::vector<Pixel> & pixels(dataMap->pixels());
   std::vector<double> energies;
   dataMap->getEnergies(energies);

   long npts = energies.size()*pixels.size();
   m_model.resize(npts, 0);

   const astro::SkyDir & dir(pointSrc->getDir());
   MeanPsf meanPsf(dir.ra(), dir.dec(), energies, m_observation);
   
   const std::vector<double> & exposure = meanPsf.exposure();
   double pixel_size = m_dataMap->pixelSize();
   
   bool galactic(dataMap->isGalactic());
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
     std::pair<double, double>(dataMap->projection().sph2pix(src_lon, 
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

   createOffsetMap(src,dataMap);

   if(::getenv("USE_ADAPTIVE_PSF_ESTIMATOR")) {
     m_psfEstimatorMethod = "adaptive";
   } else if(::getenv("USE_ANNULAR_PSF_ESTIMATOR")) {
     m_psfEstimatorMethod = "annular";
   } else if(::getenv("USE_OLD_PSF_ESTIMATOR")) {
     m_psfEstimatorMethod = "pixel_center";
   } else {
     m_psfEstimatorMethod = "adaptive";
   }

   if (::getenv("PSF_ADAPTIVE_ESTIMATOR_FTOL"))
     m_psfEstimatorFtol = atof(::getenv("PSF_ADAPTIVE_ESTIMATOR_FTOL"));

   if (::getenv("PSF_ADAPTIVE_ESTIMATOR_PEAK_TH"))
     m_psfEstimatorPeakTh = atof(::getenv("PSF_ADAPTIVE_ESTIMATOR_PEAK_TH"));

   if (performConvolution) {
      long icount(0);
      std::vector<double> mapCorrections(energies.size(), 1.);
      if (applyPsfCorrections &&
          dataMap->withinBounds(dir, energies.at(energies.size()/2), 4)) {
	getMapCorrections(pointSrc, meanPsf, pixels, srcCoord, 
			  pixCoords, energies,
			  mapCorrections);
      }

      //std::vector<Pixel>::const_iterator pixel(pixels.begin());
      pixel = pixels.begin();
      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
         std::vector<double>::const_iterator energy = energies.begin();
         for (int k = 0; energy != energies.end(); ++energy, k++) {
            unsigned long indx = k*pixels.size() + j;
            if (verbose && (icount % (npts/20)) == 0) {
               m_formatter->warn() << ".";
            }
            double offset(dir.difference(pixel->dir())*180./M_PI);
            double psf_value(psfValueEstimate(meanPsf, energies.at(k), 
                                              dir, *pixel, srcCoord,
					      pixCoords[j]));
            double value(psf_value*exposure.at(k));
            value *= pixel->solidAngle()*mapCorrections.at(k);
            m_model.at(indx) += value;
            icount++;
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
            m_model.at(indx) = exposure.at(k);
         }
      }
   }
}


void SourceMap::makePointSourceMap_healpix(Source * src,
					   const CountsMapHealpix * dataMap,
					   bool applyPsfCorrections,
					   bool performConvolution,
					   bool verbose) {
  PointSource * pointSrc = dynamic_cast<PointSource *>(src);
  const std::vector<Pixel> & pixels(dataMap->pixels());
  std::vector<double> energies;
  dataMap->getEnergies(energies);
  
  long npts = energies.size()*pixels.size();
  m_model.resize(npts, 0);

  const astro::SkyDir & dir(pointSrc->getDir());
  MeanPsf meanPsf(dir.ra(), dir.dec(), energies, m_observation);  
  const std::vector<double> & exposure = meanPsf.exposure();

  int nPix = dataMap->nPixels();

  if (performConvolution) {    
    size_t indx(0);
    const Healpix_Base& hp = dataMap->healpixProj()->healpix();
    Healpix_Map<float> hpmap(hp.Nside(),hp.Scheme(),SET_NSIDE);
    for (int k(0); k < energies.size(); k++ ) {
      m_formatter->warn() << ".";
      hpmap.fill(0.);
      ConvolveHealpix::fillMapWithPSF_refDir(meanPsf,energies[k],dir,dataMap->isGalactic(),hpmap);
      for (int i(0); i < nPix; i++, indx++ ) {
	m_model.at(indx) = exposure.at(k) * hpmap[ dataMap->localToGlobalIndex(i) ];
      }
    }      
  } else {
    dataMap->projection();
    double ipix(0);
    double jpix(0);
    st_facilities::Util::skyDir2pixel(dataMap->projection(),dir,ipix,jpix);
    int iglo = dataMap->globalToLocalIndex(ipix);
    if ( iglo < 0 ) return;
    if ( ipix >= 0 ) {
      std::vector<double>::const_iterator energy = energies.begin();
      size_t indx = size_t(ipix);
      for (int k = 0; energy != energies.end(); ++energy, ++k ) {
	size_t indx = k*pixels.size() + iglo;
	m_model.at(indx) = dataMap->solidAngle() * exposure.at(k);
      }      
    }
  }
}
  
SourceMap::~SourceMap() {
   if (m_deleteDataMap) {
      delete m_dataMap;
   }
   delete m_formatter;
}

void SourceMap::createOffsetMap(Source * src, const CountsMap * dataMap) {

  PointSource * pointSrc = dynamic_cast<PointSource *>(src);
  const astro::SkyDir & dir(pointSrc->getDir());
  const std::vector<Pixel> & pixels(dataMap->pixels());
  const double pixel_size = m_dataMap->pixelSize();

  m_pixelOffset.resize(dataMap->naxis2(),
		      std::vector<double>(dataMap->naxis1()));
  m_pixelCoordX.resize(dataMap->naxis1());
  m_pixelCoordY.resize(dataMap->naxis2());

  for( int j = 0; j < dataMap->naxis1(); j++)
    m_pixelCoordX[j] = double(j+1);

  for( int j = 0; j < dataMap->naxis2(); j++)
    m_pixelCoordY[j] = double(j+1);

   bool galactic(dataMap->isGalactic());
   /// Get the pixel center in pixel coordinates
   double src_lon, src_lat;
   if (galactic) {
      src_lon = dir.l();
      src_lat = dir.b();
   } else {
      src_lon = dir.ra();
      src_lat = dir.dec();
   }

  std::pair<double, double> src_coords(dataMap->projection().sph2pix(src_lon, 
								     src_lat));

  std::vector<Pixel>::const_iterator pixel(pixels.begin());
  for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
	
    int ix = j/dataMap->naxis1();
    int iy = j%dataMap->naxis1();
	
    double pix_sep = 
      pixel_size*std::sqrt(std::pow(src_coords.first-(iy+1),2) +  
			   std::pow(src_coords.second-(ix+1),2));
    double ang_sep = dir.difference(pixel->dir())*180./M_PI;
    m_pixelOffset[ix][iy] = ang_sep/pix_sep-1.0;

  }
}

void SourceMap::getMapCorrections(PointSource * src, const MeanPsf & meanPsf,
                                  const std::vector<Pixel> & pixels,
				  const std::pair<double,double>& srcCoord,
				  const std::vector< std::pair<double,double> > & pixCoords,
                                  const std::vector<double> & energies,
                                  std::vector<double> & mapCorrections) const {
   std::ofstream * output(0);
   char * log_file(::getenv("PSF_INTEGRANDS_LOG"));
   if (log_file) {
      std::cout << "log_file: " << log_file << std::endl;
      output = new std::ofstream(log_file);
   }
   const astro::SkyDir & srcDir = src->getDir();
   double psfRadius(maxPsfRadius(src));
   std::vector<unsigned int> containedPixels;
   for (unsigned int j = 0; j < pixels.size(); j++) {
      if (srcDir.difference(pixels.at(j).dir())*180./M_PI <= psfRadius) {
         containedPixels.push_back(j);
      }
   }
   mapCorrections.clear();
   mapCorrections.reserve(energies.size());
   for (unsigned int k = 0; k < energies.size()-1; k++) {
      double map_integral(0);
      std::vector<unsigned int>::const_iterator j = containedPixels.begin();
      for ( ; j != containedPixels.end(); ++j) {
         const Pixel & pix = pixels.at(*j);
         double solid_angle(pix.solidAngle());
         double psf_value(psfValueEstimate(meanPsf, energies.at(k), srcDir, 
					   pix, srcCoord, pixCoords[*j]));
         map_integral += solid_angle*psf_value;
         if (output) {
            double offset(srcDir.difference(pix.dir())*180./M_PI);
            try {
               *output << energies[k] << "  "
                       << offset << "  "
                       << solid_angle << "  "
                       << psf_value << "  "
                       << meanPsf.integral(offset, energies[k]) << std::endl;
            } catch (std::out_of_range &) {
            }
         }
      }
      if (map_integral == 0) {
         /// source effectively lies on map boundary, so apply no
         /// correction
         mapCorrections.push_back(1);
      } else {
         /// Correct for undersampling of the PSF at high eneriges.
         double value(meanPsf.integral(psfRadius, energies.at(k))/map_integral);
         mapCorrections.push_back(value);
      }
   }
   if (output) {
      output->close();
   }
   mapCorrections.push_back(mapCorrections.back());
}

double SourceMap::maxPsfRadius(PointSource * src) const {
   std::vector<astro::SkyDir> pixelDirs;
   m_dataMap->getBoundaryPixelDirs(pixelDirs);

   const astro::SkyDir & srcDir = src->getDir();
   double radius = srcDir.difference(pixelDirs.at(0));
   for (unsigned int i = 1; i < pixelDirs.size(); i++) {
      double new_rad = srcDir.difference(pixelDirs.at(i));
      if (new_rad < radius ) {
         radius = new_rad;
      }
   }
   return radius*180./M_PI;
}

void SourceMap::computeNpredArray() {
   const std::vector<Pixel> & pixels(m_dataMap->pixels());
   
   std::vector<double> energies;
   switch ( m_dataMap->projection().method() ) {
   case astro::ProjBase::WCS:
     m_dataMap->getAxisVector(2, energies);
     break;
   case astro::ProjBase::HEALPIX:
     m_dataMap->getAxisVector(1, energies);
     break;
   }
   
   m_npreds.clear();
   m_npreds.resize(energies.size(), 0);
   m_weightedNpreds.clear();
   m_weightedNpreds.resize(energies.size(), 0);
   for (size_t k(0); k < energies.size(); k++) {
      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
         size_t indx(k*pixels.size() + j);
	 double addend = m_model.at(indx);
         m_npreds[k] += addend;
	 addend *= m_weights != 0 ? m_weights->model()[indx] : 1.0;
	 m_weightedNpreds[k] += addend;
      }
   }
}

bool SourceMap::haveMapCubeFunction(DiffuseSource * src) const {
   Source::FuncMap & srcFuncs = src->getSrcFuncs();
   return srcFuncs["SpatialDist"]->genericName() == "MapCubeFunction";
}

double SourceMap::computeResampFactor(const DiffuseSource & src,
                                      const CountsMapBase & dataMap) const {
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

void SourceMap::makeProjectedMap(const ProjMap& weight_map) {
   const std::vector<Pixel> & pixels(m_dataMap->pixels());
   std::vector<double> energies;
   m_dataMap->getEnergies(energies);
   m_model.resize(energies.size()*pixels.size());
   for (size_t k(0); k < energies.size(); k++) {
      std::vector<Pixel>::const_iterator pixel(pixels.begin());
      for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
         size_t indx(k*pixels.size() + j);
	 try {
	   m_model.at(indx) = weight_map.operator()(pixel->dir(),
						    energies[k]);
	 } catch (...) {
	   m_model.at(indx) = 1.0;
	 }	 
      }
   }
}    

void SourceMap::applyPhasedExposureMap() {
   if (!m_observation.have_phased_expmap()) {
      return;
   }
   const ProjMap * phased_expmap = &(m_observation.phased_expmap());
   const std::vector<Pixel> & pixels(m_dataMap->pixels());
   std::vector<double> energies;
   m_dataMap->getEnergies(energies);
   for (size_t k(0); k < energies.size(); k++) {
      std::vector<Pixel>::const_iterator pixel(pixels.begin());
      for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
         size_t indx(k*pixels.size() + j);
         m_model.at(indx) *= phased_expmap->operator()(pixel->dir(),
                                                       energies[k]);
      }
   }
}

double SourceMap::
psfValueEstimate(const MeanPsf & meanPsf, double energy, 
                 const astro::SkyDir & srcDir, 
                 const Pixel & pixel,
		 const std::pair<double, double>& srcCoord,
		 const std::pair<double, double>& pixCoord) const {
   
   if (m_psfEstimatorMethod == "adaptive") {
     // Use the adaptive pixel integrator
     return integrate_psf_adaptive(meanPsf, energy, srcDir, pixel,
				   srcCoord, pixCoord);
   } else if(m_psfEstimatorMethod == "pixel_center") {

     double offset(srcDir.difference(pixel.dir())*180./M_PI);
     double pixelSolidAngle(pixel.solidAngle());
   
     // Use the central pixel value as in the previous implementation 
     // (ST 09-33-00)
     return meanPsf(energy, offset);
   } else if(m_psfEstimatorMethod == "annular") {

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
   } else {
     std::string errMsg(" SourceMap::psfValueEstimate: Unknown PSF estimator: ");
     errMsg += m_psfEstimatorMethod;
     throw std::runtime_error(errMsg);
   }
}

double SourceMap::integrate_psf_adaptive(const MeanPsf & meanPsf, 
					 double energy, 
					 const astro::SkyDir & srcDir, 
					 const Pixel & pixel, 
					 const std::pair<double, double>& srcCoord,
					 const std::pair<double, double>& pixCoord) const {
    
  double offset(srcDir.difference(pixel.dir())*180./M_PI);
  double pixel_size = m_dataMap->pixelSize();
  const int max_npts = 64;
  double ftol_threshold = m_psfEstimatorFtol;
  double peak_threshold = m_psfEstimatorPeakTh;

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
      v1 = integrate_psf_fast(meanPsf, energy, srcDir, pixel, srcCoord, pixCoord, npts);

      if(v1-v0 == 0 || v1 == 0)
	break;

      ferr = std::fabs((v1-v0)/v1);
      if(ferr < ftol_threshold || npts >= max_npts)
	break;

      v0=v1;
  }

  return v1;
}

double SourceMap::integrate_psf_fast(const MeanPsf & meanPsf, double energy, 
				     const astro::SkyDir & srcDir, 
				     const Pixel & pixel, 
				     const std::pair<double, double>& srcCoord,
				     const std::pair<double, double>& pixCoord, 
				     size_t npts) const {

   double pixel_size = m_dataMap->pixelSize();
   double psf_value(0);
   double dstep(1./static_cast<double>(npts));
   for (size_t i(0); i < npts; i++) {
      double x(pixCoord.first + i*dstep - 0.5 + 0.5*dstep);
      for (size_t j(0); j < npts; j++) {
         double y(pixCoord.second + j*dstep - 0.5 + 0.5*dstep);
	 double pix_offset = 
	   pixel_size*std::sqrt(std::pow(srcCoord.first-x,2) +  
				std::pow(srcCoord.second-y,2));
	 double scale = st_facilities::Util::bilinear(m_pixelCoordY, y, 
						      m_pixelCoordX, x, 
						      m_pixelOffset);

	 psf_value += meanPsf(energy, pix_offset*(1+scale));
      }
   }
   return psf_value/static_cast<double>(npts*npts);
}

double SourceMap::
integrate_psf(const MeanPsf & meanPsf, double energy,
              const astro::SkyDir & srcDir, const Pixel & pixel, size_t npts) const {
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

void SourceMap::setImage(const std::vector<float>& model) {
  if(model.size() != m_model.size())
    throw std::runtime_error("Wrong size for input model map.");

  m_model = model;
  applyPhasedExposureMap();
  computeNpredArray();
}

void SourceMap::setWeights(const SourceMap* weights) {
  m_weights = weights;
  computeNpredArray();
}

void SourceMap::readImage(const std::string& sourceMapsFile) {
  std::auto_ptr<const tip::Image>
  image(tip::IFileSvc::instance().readImage(sourceMapsFile, m_name));
  m_model.clear();
  image->get(m_model);
}

void SourceMap::readTable_healpix(const std::string& sourceMapsFile) {
  std::auto_ptr<const tip::Table> 
    table(tip::IFileSvc::instance().readTable(sourceMapsFile,m_name));

  // This is a bit tricky, basically all the data we care about
  // are in the columns called "CHANNELx"
  // Note also that tip work in lowercase
  std::vector<tip::FieldIndex_t> dataColumns;
  const tip::Table::FieldCont& colNames = table->getValidFields();
  for ( tip::Table::FieldCont::const_iterator itr = colNames.begin(); 
	itr != colNames.end(); itr++ ) {
    if ( itr->find("channel") == 0 ) { 
      dataColumns.push_back( table->getFieldIndex(*itr) );     
    } else {
      continue;
    }
  }
  
  int ncol = dataColumns.size();
  tip::Index_t nrow = table->getNumRecords();
  m_model.clear();
  m_model.reserve(ncol*nrow);

  double readVal(0.);
  for ( std::vector<tip::FieldIndex_t>::const_iterator itrData = dataColumns.begin();
	itrData != dataColumns.end(); itrData++ ) {    
    const tip::IColumn* col = table->getColumn(*itrData);
    for ( tip::Index_t irow(0); irow < nrow; irow++ ) {
      col->get(irow,readVal);
      m_model.push_back(readVal);
    }
  }
}


} // Likelihood
