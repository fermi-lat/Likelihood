/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.104 2014/05/20 21:53:28 jchiang Exp $
 */

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
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

#include "Likelihood/WcsMap2.h"

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

SourceMap::SourceMap(Source * src, const CountsMap * dataMap,
                     const Observation & observation, 
                     bool applyPsfCorrections,
                     bool performConvolution,
                     bool resample,
                     double resamp_factor,
                     double minbinsz,
                     bool verbose)
   : m_name(src->getName()),
     m_srcType(src->getType()),
     m_dataMap(dataMap),
     m_observation(observation),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_deleteDataMap(false) {
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

SourceMap::SourceMap(const std::string & sourceMapsFile,
                     const std::string & srcName,
                     const Observation & observation) 
   : m_name(srcName),
     m_dataMap(new CountsMap(sourceMapsFile)),
     m_observation(observation),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_deleteDataMap(true) {
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(sourceMapsFile, srcName));
   m_model.clear();
   image->get(m_model);
   applyPhasedExposureMap();
   computeNpredArray();
}

void SourceMap::makeDiffuseMap(Source * src, 
                               const CountsMap * dataMap,
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
      double cdelt1 = std::fabs(tmp.wcsmap().cdelt1());
      double cdelt2 = std::fabs(tmp.wcsmap().cdelt2());
      if (cdelt1 < minbinsz || cdelt2 < minbinsz) {
         unsigned int factor = 
            std::max(static_cast<unsigned int>(minbinsz/cdelt1),
                     static_cast<unsigned int>(minbinsz/cdelt2));
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
   dataMap->getAxisVector(2, energies);

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
   if (dataMap->conformingMap()) {
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
   size_t counter(0);
   std::vector<double>::const_iterator energy = energies.begin();
   for (int k(0); energy != energies.end(); ++energy, k++) {
      bool interpolate;
      WcsMap2 diffuseMap(*diffuseSrc, mapRefDir.ra(), mapRefDir.dec(),
                         crpix1, crpix2, cdelt1, cdelt2, naxis1, naxis2,
                         *energy, dataMap->proj_name(), 
                         dataMap->projection().isGalactic(), 
                         interpolate=true);
      const MeanPsf & meanpsf(m_observation.meanpsf());
      const BinnedExposure & bexpmap(m_observation.bexpmap());
      WcsMap2 convolvedMap(diffuseMap.convolve(*energy, meanpsf, 
                                               bexpmap, performConvolution));
      size_t rfac(static_cast<size_t>(resamp_factor));
      double solid_angle;
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
            m_model[indx] += (convolvedMap.image()[0][j][i]
                              /resamp_factor/resamp_factor
                              *solid_angle);
         }
      }
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
                                   const CountsMap * dataMap,
                                   bool applyPsfCorrections,
                                   bool performConvolution,
                                   bool verbose) {
   PointSource * pointSrc = dynamic_cast<PointSource *>(src);
   const std::vector<Pixel> & pixels(dataMap->pixels());
   std::vector<double> energies;
   dataMap->getAxisVector(2, energies);

   long npts = energies.size()*pixels.size();
   m_model.resize(npts, 0);

   std::vector<Pixel>::const_iterator pixel = pixels.begin();

   const astro::SkyDir & dir(pointSrc->getDir());
   MeanPsf meanPsf(dir.ra(), dir.dec(), energies, m_observation);
   
   const std::vector<double> & exposure = meanPsf.exposure();
   double pixel_size(std::min(std::fabs(m_dataMap->cdelt1()), 
                              std::fabs(m_dataMap->cdelt2())));
   
   if (performConvolution) {
      long icount(0);
      std::vector<double> mapCorrections(energies.size(), 1.);
      if (applyPsfCorrections &&
          dataMap->withinBounds(dir, energies.at(energies.size()/2))) {
            getMapCorrections(pointSrc, meanPsf, pixels, energies,
                              mapCorrections);
      }

      std::vector<Pixel>::const_iterator pixel(pixels.begin());
      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
         std::vector<double>::const_iterator energy = energies.begin();
         for (int k = 0; energy != energies.end(); ++energy, k++) {
            unsigned long indx = k*pixels.size() + j;
            if (verbose && (icount % (npts/20)) == 0) {
               m_formatter->warn() << ".";
            }
            double offset(dir.difference(pixel->dir())*180./M_PI);
            double psf_value(psfValueEstimate(meanPsf, energies.at(k),
                                              offset, pixel->solidAngle()));
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

SourceMap::~SourceMap() {
   if (m_deleteDataMap) {
      delete m_dataMap;
   }
   delete m_formatter;
}

void SourceMap::getMapCorrections(PointSource * src, const MeanPsf & meanPsf,
                                  const std::vector<Pixel> & pixels,
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
         double offset(srcDir.difference(pix.dir())*180./M_PI);
         double psf_value(psfValueEstimate(meanPsf, energies.at(k), offset,
                                           pix.solidAngle()));
         map_integral += solid_angle*psf_value;
         if (output) {
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
   m_dataMap->getAxisVector(2, energies);

   m_npreds.resize(energies.size(), 0);
   for (size_t k(0); k < energies.size(); k++) {
      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      for (size_t j(0); pixel != pixels.end(); ++pixel, j++) {
         size_t indx(k*pixels.size() + j);
         m_npreds.at(k) += m_model.at(indx);
      }
   }
}

bool SourceMap::haveMapCubeFunction(DiffuseSource * src) const {
   Source::FuncMap & srcFuncs = src->getSrcFuncs();
   return srcFuncs["SpatialDist"]->genericName() == "MapCubeFunction";
}

double SourceMap::computeResampFactor(const DiffuseSource & src,
                                      const CountsMap & dataMap) const {
   double data_pixel_size = std::min(std::fabs(dataMap.cdelt1()), 
                                     std::fabs(dataMap.cdelt2()));
   double model_pixel_size = data_pixel_size;
   try {
      model_pixel_size = 
         std::min(std::fabs(src.mapBaseObject()->wcsmap().cdelt1()),
                  std::fabs(src.mapBaseObject()->wcsmap().cdelt2()));
   } catch (MapBaseException &) {
      // do nothing
   }
   double resamp_factor = 
      std::max(2, static_cast<int>(data_pixel_size/model_pixel_size));
   return resamp_factor;
}

void SourceMap::applyPhasedExposureMap() {
   const WcsMap2 * phased_expmap(&(m_observation.phased_expmap()));
   if (phased_expmap == 0) {
      return;
   }
   const std::vector<Pixel> & pixels(m_dataMap->pixels());
   std::vector<double> energies;
   m_dataMap->getAxisVector(2, energies);
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
                 double offset, double pixelSolidAngle) const {
/// To estimate the psf value averaged over a pixel, average the psf
/// over an annulus centered on the source position with approximately
/// the same extent in theta as the pixel in question.
   if (::getenv("USE_OLD_PSF_ESTIMATOR")) {
      // Use the central pixel value as in the previous implementation 
      // (ST 09-33-00)
      return meanPsf(energy, offset);
   }
   double pixel_value(0);
   double pixel_size(std::sqrt(pixelSolidAngle)*180./M_PI);
   if (pixel_size/2. >= offset) {
      /// Average over an acceptance cone with the same solid angle
      /// as the central pixel.
      double radius(std::acos(1. - pixelSolidAngle/2./M_PI)*180./M_PI);
      pixel_value = meanPsf.integral(radius, energy)/pixelSolidAngle;
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

} // namespace Likelihood
