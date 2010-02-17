/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.74 2010/02/17 07:12:16 jchiang Exp $
 */

#include <algorithm>
#include <deque>
#include <iostream>
#include <memory>

#include "st_stream/StreamFormatter.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/Source.h"
#define ST_DLL_EXPORTS
#include "Likelihood/SourceMap.h"
#undef ST_DLL_EXPORTS
#include "Likelihood/TrapQuad.h"

#include "Likelihood/WcsMap.h"

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

/// @bug Use of this function (by CountsSpectra.cxx) is necessitated by
/// annoying linkage problems on Windows.
void getNpreds(const SourceMap & srcMap, std::vector<double> & npreds) {
   npreds = srcMap.npreds();
}

std::string SourceMap::s_expMapFileName;
MeanPsf * SourceMap::s_meanPsf(0);
BinnedExposure * SourceMap::s_binnedExposure(0);
unsigned int SourceMap::s_refCount(0);

std::vector<double> SourceMap::s_phi;
std::vector<double> SourceMap::s_mu;
std::vector<double> SourceMap::s_theta;

SourceMap::SourceMap(Source * src, const CountsMap * dataMap,
                     const Observation & observation, 
                     bool applyPsfCorrections,
                     bool performConvolution,
                     bool resample,
                     double resamp_factor)
   : m_name(src->getName()), m_srcType(src->getType()),
     m_dataMap(dataMap), 
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_deleteDataMap(false) {
   s_refCount++;
   if (s_mu.size() == 0 || s_phi.size() == 0 || s_theta.size() == 0) {
      prepareAngleArrays();
   }

   const std::vector<Pixel> & pixels(dataMap->pixels());
   
   std::vector<double> energies;
   dataMap->getAxisVector(2, energies);

   m_formatter->warn() << "Generating SourceMap for " << m_name;
   long npts = energies.size()*pixels.size();
   m_model.resize(npts, 0);
   long icount(0);

   m_npreds.resize(energies.size(), 0);

   std::vector<Pixel>::const_iterator pixel = pixels.begin();

   bool havePointSource = dynamic_cast<PointSource *>(src) != 0;
   bool haveDiffuseSource = dynamic_cast<DiffuseSource *>(src) != 0;

   if (haveDiffuseSource) {
      computeExposureAndPsf(observation);
      DiffuseSource * diffuseSrc = dynamic_cast<DiffuseSource *>(src);
      const astro::SkyDir & mapRefDir = dataMap->refDir();
      if (!resample) {
         resamp_factor = 1;
      }
      double crpix1, crpix2;
      int naxis1, naxis2;
      double cdelt1 = dataMap->cdelt1()/resamp_factor;
      double cdelt2 = dataMap->cdelt2()/resamp_factor;
      if (dataMap->conformingMap()) {
         double radius = std::min(180., ::maxRadius(pixels, mapRefDir) + 10.);
         // Conforming maps have abs(CDELT1) == abs(CDELT2).  This
         // expression for the mapsize ensures that the number of
         // pixels in each dimension is even.
         int mapsize(2*static_cast<int>(radius/std::fabs(cdelt1)));
         naxis1 = mapsize;
         naxis2 = mapsize;
         crpix1 = (naxis1 + 1.)/2.;
         crpix2 = (naxis2 + 1.)/2.;
         if (!resample) { 
            // Use integer or half-integer reference pixel based on
            // input counts map, even though naxis1 and naxis2 both
            // must be even.
            if (dataMap->naxis1() % 2 == 1) {
               crpix1 += 0.5;
            }
            if (dataMap->naxis2() % 2 == 1) {
               crpix2 += 0.5;
            }
         }
      } else {
         // The counts map was not created by gtbin, so just adopt the
         // map geometry without adding padding for psf leakage since
         // this cannot be done in general without redefining the
         // reference pixel and reference direction.
         naxis1 = dataMap->naxis1()*resamp_factor;
         naxis2 = dataMap->naxis2()*resamp_factor;
         // Ensure an even number of pixels in each direction.
         if (naxis1 % 2 == 1) {
            naxis1 += 1;
         }
         if (naxis2 % 2 == 1) {
            naxis2 += 1;
         }
         crpix1 = dataMap->crpix1()*resamp_factor;
         crpix2 = dataMap->crpix2()*resamp_factor;
      }
      unsigned int indx(0);
      std::vector<double>::const_iterator energy = energies.begin();
      for (int k(0); energy != energies.end(); ++energy, k++) {
         WcsMap diffuseMap(*diffuseSrc, mapRefDir.ra(), mapRefDir.dec(),
                           crpix1, crpix2, cdelt1, cdelt2, naxis1, naxis2,
                           *energy, dataMap->proj_name(), 
                           dataMap->projection().isGalactic(), true);
         WcsMap convolvedMap(diffuseMap.convolve(*energy, *s_meanPsf, 
                                                 *s_binnedExposure,
                                                 performConvolution));
         for (pixel = pixels.begin(); pixel != pixels.end();
              ++pixel, indx++) {
            if ((indx % (npts/20)) == 0) {
               m_formatter->warn() << ".";
            }
            if (pixel->solidAngle() > 0) {
               m_model.at(indx) = (convolvedMap(pixel->dir())
                                   *pixel->solidAngle());
               m_npreds.at(k) += m_model.at(indx);
            }
         }
      }
   } else if (havePointSource) {
      PointSource * pointSrc = dynamic_cast<PointSource *>(src);

      const astro::SkyDir & dir(pointSrc->getDir());
      MeanPsf meanPsf(dir.ra(), dir.dec(), energies, observation);

      const std::vector<double> & exposure = meanPsf.exposure();

      if (performConvolution) {
         std::vector<double> mapCorrections(energies.size(), 1.);
         if (applyPsfCorrections &&
             dataMap->withinBounds(dir, energies.at(energies.size()/2))) {
            getMapCorrections(pointSrc, meanPsf, pixels, energies,
                              mapCorrections);
         }

         for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
            std::vector<double>::const_iterator energy = energies.begin();
            for (int k = 0; energy != energies.end(); ++energy, k++) {
               unsigned long indx = k*pixels.size() + j;
               if ((icount % (npts/20)) == 0) {
                  m_formatter->warn() << ".";
               }
               double value = (meanPsf(energies.at(k), 
                                       dir.difference(pixel->dir())*180./M_PI)
                               *exposure.at(k));
               value *= pixel->solidAngle()*mapCorrections.at(k);
               m_model.at(indx) += value;
               m_npreds.at(k) += value;
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
               m_npreds.at(k) = m_model.at(indx);
            }
         }
      }
   }
   m_formatter->warn() << "!" << std::endl;
}

SourceMap::~SourceMap() {
   s_refCount--;
   if (s_refCount == 0) {
      delete s_meanPsf;
      s_meanPsf = 0;
      delete s_binnedExposure;
      s_binnedExposure = 0;
   }
   if (m_deleteDataMap) {
      delete m_dataMap;
   }
   delete m_formatter;
}

void SourceMap::setBinnedExposure(const std::string & filename) {
   if (s_binnedExposure != 0) {
      delete s_binnedExposure;
      s_binnedExposure = 0;
   }
   s_binnedExposure = new BinnedExposure(filename);
   s_expMapFileName = filename;
}

void SourceMap::getMapCorrections(PointSource * src, const MeanPsf & meanPsf,
                                  const std::vector<Pixel> & pixels,
                                  const std::vector<double> & energies,
                                  std::vector<double> & mapCorrections) const {
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
         map_integral += pix.solidAngle()*
            meanPsf(energies.at(k), srcDir.difference(pix.dir())*180./M_PI);
      }
      if (map_integral == 0) {
// source effectively lies on map boundary, so apply no correction
         mapCorrections.push_back(1);
      } else {
         mapCorrections.push_back(meanPsf.integral(psfRadius, energies.at(k))
                                  /map_integral);
      }
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

SourceMap::SourceMap(const std::string & sourceMapsFile,
                     const std::string & srcName) 
   : m_name(srcName), m_dataMap(new CountsMap(sourceMapsFile)),
     m_formatter(new st_stream::StreamFormatter("SourceMap", "", 2)),
     m_deleteDataMap(true) {
   s_refCount++;
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(sourceMapsFile, srcName));
   m_model.clear();
   image->get(m_model);

   const std::vector<Pixel> & pixels(m_dataMap->pixels());
   
   std::vector<double> energies;
   m_dataMap->getAxisVector(2, energies);

   m_npreds.resize(energies.size(), 0);
   for (unsigned int k = 0; k < energies.size(); k++) {
      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
         unsigned long indx = k*pixels.size() + j;
         m_npreds.at(k) += m_model.at(indx);
      }
   }
   if (s_mu.size() == 0 || s_phi.size() == 0 || s_theta.size() == 0) {
      prepareAngleArrays();
   }
}

bool SourceMap::haveMapCubeFunction(DiffuseSource * src) const {
   Source::FuncMap & srcFuncs = src->getSrcFuncs();
   return srcFuncs["SpatialDist"]->genericName() == "MapCubeFunction";
}

void SourceMap::computeExposureAndPsf(const Observation & observation) {
   std::vector<double> energies;
   m_dataMap->getAxisVector(2, energies);
   if (s_meanPsf == 0) {
      double ra = m_dataMap->refDir().ra();
      double dec = m_dataMap->refDir().dec();
      s_meanPsf = new MeanPsf(ra, dec, energies, observation);
   }
   if (s_binnedExposure == 0) {
      st_stream::StreamFormatter formatter("SourceMap", "computeMap", 2);
      if (s_expMapFileName == "" || s_expMapFileName == "none") {
         formatter.info() << "\nBinned exposure map not specified. "
                          << "\nEnter filename of existing map or "
                          << "the name to be used for a new map: ";
         std::cin >> s_expMapFileName;
      }
      try {
         s_binnedExposure = new BinnedExposure(s_expMapFileName);
      } catch (tip::TipException &) {
         s_binnedExposure = new BinnedExposure(energies, observation);
         s_binnedExposure->writeOutput(s_expMapFileName);
      }
   }
}

void SourceMap::prepareAngleArrays(int nmu, int nphi) {
   double radius = 30.;
   double mumin = cos(radius*M_PI/180);

// Sample more densely near theta = 0:
   std::deque<double> my_mu;
   double nscale = static_cast<double>((nmu-1)*(nmu-1));
   for (int i = 0; i < nmu; i++) {
      my_mu.push_front(1. - i*i/nscale*(1. - mumin));
   }
   s_mu.resize(my_mu.size());
   std::copy(my_mu.begin(), my_mu.end(), s_mu.begin());

   s_theta.resize(s_mu.size());
   for (unsigned int i = 0; i < s_mu.size(); i++) {
      s_theta.at(i) = acos(s_mu.at(i))*180./M_PI;
   }

   s_phi.clear();
   double phistep = 2.*M_PI/(nphi - 1.);
   for (int i = 0; i < nphi; i++) {
      s_phi.push_back(phistep*i);
   }
}
   void SourceMap::setBinnedExpMapName(const std::string & filename) {
      s_expMapFileName = filename;
   }

    const std::string & SourceMap::binnedExpMap() {
      return s_expMapFileName;
   }
} // namespace Likelihood
