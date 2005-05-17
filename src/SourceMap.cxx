/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.32 2005/05/17 00:26:40 jchiang Exp $
 */

#include <algorithm>
#include <deque>
#include <memory>

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

#include "st_facilities/Util.h"

#include "Likelihood/BinnedExposure.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/FitsImage.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/Observation.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/TrapQuad.h"

#include "Verbosity.h"

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
}

namespace Likelihood {

#include "fitsio.h"

MeanPsf * SourceMap::s_meanPsf(0);
BinnedExposure * SourceMap::s_binnedExposure(0);
unsigned int SourceMap::s_refCount(0);

std::vector<double> SourceMap::s_phi;
std::vector<double> SourceMap::s_mu;
std::vector<double> SourceMap::s_theta;

SourceMap::SourceMap(Source * src, const CountsMap * dataMap,
                     const Observation & observation) 
   : m_name(src->getName()), m_dataMap(dataMap), m_deleteDataMap(false) {
   s_refCount++;
   if (s_mu.size() == 0 || s_phi.size() == 0 || s_theta.size() == 0) {
      prepareAngleArrays();
   }

   std::vector<Pixel> pixels;
   dataMap->getPixels(pixels);
   
   std::vector<double> energies;
   dataMap->getAxisVector(2, energies);

   if (print_output()) std::cerr << "Generating SourceMap for " << m_name;
   long npts = energies.size()*pixels.size();
   m_model.resize(npts, 0);
   long icount(0);

   m_npreds.resize(energies.size(), 0);

   std::vector<Pixel>::const_iterator pixel = pixels.begin();

   bool havePointSource = dynamic_cast<PointSource *>(src) != 0;
   bool haveDiffuseSource = dynamic_cast<DiffuseSource *>(src) != 0;

   if (haveDiffuseSource) {
      DiffuseSource * diffuseSrc = dynamic_cast<DiffuseSource *>(src);
      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
         computeSrcDirs(*pixel, src);
         std::vector<double>::const_iterator energy = energies.begin();
         for (int k = 0; energy != energies.end(); ++energy, k++) {
            unsigned long indx = k*pixels.size() + j;
            if (print_output() && (icount % (npts/20)) == 0) std::cerr << ".";
            double value(0);
            if (haveMapCubeFunction(diffuseSrc)) {
               recomputeSrcStrengths(diffuseSrc, *energy);
            }
            value = sourceRegionIntegral(*energy, observation);
            value *= pixel->solidAngle();
            m_model.at(indx) += value;
            m_npreds.at(k) += value;
            icount++;
         }
      }
   } else if (havePointSource) {
      PointSource * pointSrc = dynamic_cast<PointSource *>(src);

      const astro::SkyDir & dir(pointSrc->getDir());
      MeanPsf meanPsf(dir.ra(), dir.dec(), energies, observation);

      std::vector<double> exposure;
      for (unsigned int k = 0; k < energies.size(); k++) {
         exposure.push_back(meanPsf.exposure(energies.at(k)));
      }

      std::vector<double> mapCorrections(energies.size(), 1.);
      if (dataMap->withinBounds(dir, energies.at(energies.size()/2))) {
         getMapCorrections(pointSrc, meanPsf, pixels, energies,
                           mapCorrections);
      }

      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
         std::vector<double>::const_iterator energy = energies.begin();
         for (int k = 0; energy != energies.end(); ++energy, k++) {
            unsigned long indx = k*pixels.size() + j;
            if (print_output() && (icount % (npts/20)) == 0) {
               std::cerr << ".";
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
   }
   if (print_output()) {
      std::cerr << "!" << std::endl;
   }
}

SourceMap::~SourceMap() {
   s_refCount--;
   if (s_refCount == 0) {
      delete s_meanPsf;
      s_meanPsf = 0;
      delete s_binnedExposure;
      s_binnedExposure = 0;
   }
   if (m_deleteDataMap) delete m_dataMap;
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
   mapCorrections.reserve(energies.size());
   for (unsigned int k = 0; k < energies.size()-1; k++) {
      double map_integral(0);
      std::vector<unsigned int>::const_iterator j = containedPixels.begin();
      for ( ; j != containedPixels.end(); ++j) {
         const Pixel & pix = pixels.at(*j);
         map_integral += pix.solidAngle()*
            meanPsf(energies.at(k), srcDir.difference(pix.dir())*180./M_PI);
      }
      mapCorrections.push_back(meanPsf.integral(psfRadius, energies.at(k))
                               /map_integral);
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
     m_deleteDataMap(true) {
   s_refCount++;
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(sourceMapsFile, srcName));
   std::vector<float> image_data;
   image->get(image_data);
   m_model.resize(image_data.size());
   std::copy(image_data.begin(), image_data.end(), m_model.begin());

   std::vector<Pixel> pixels;
   m_dataMap->getPixels(pixels);
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

void SourceMap::save(const std::string & filename) const {
   if (st_facilities::Util::fileExists(filename)) {
      throw std::runtime_error("SourceMap::save: " + filename 
                               + " already exists.");
   }

   fitsfile * fptr;
   int status(0);

   fits_create_file(&fptr, filename.c_str(), &status);
   fitsReportError(stderr, status);

   long naxes[] = {m_dataMap->imageDimension(0),
                   m_dataMap->imageDimension(1),
                   m_dataMap->imageDimension(2) + 1};
   fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &status);
   fitsReportError(stderr, status);
   
   long group(0);
   const std::vector<double> & data = model();
   fits_write_3d_dbl(fptr, group, naxes[0], naxes[1], naxes[0], naxes[1],
                     naxes[2], const_cast<double *>(&data[0]), &status);
   fitsReportError(stderr, status);

   fits_update_key(fptr, TSTRING, "EXTNAME", 
                   const_cast<char *>(m_name.c_str()), "SourceMap name",
                   &status);
   fitsReportError(stderr, status);

   fits_close_file(fptr, &status);
   fitsReportError(stderr, status);
}

void SourceMap::fitsReportError(FILE *stream, int status) const {
   if (status != 0) {
      fits_report_error(stream, status);
      throw std::runtime_error("SourceMap::save: cfitsio error.");
   }
}

double SourceMap::Aeff::operator()(double costheta) const {
   double inclination = acos(costheta)*180./M_PI;
   static double phi(0);
   return m_observation.respFuncs().totalResponse(inclination, phi, m_energy,
                                                  m_energy, m_separation,
                                                  m_type);
}

bool SourceMap::haveMapCubeFunction(DiffuseSource * src) const {
   Source::FuncMap & srcFuncs = src->getSrcFuncs();
   return srcFuncs["SpatialDist"]->genericName() == "MapCubeFunction";
}

double SourceMap::sourceRegionIntegral(double energy,
                                       const Observation & observation) const {
   std::vector<double> energies;
   m_dataMap->getAxisVector(2, energies);
   if (s_meanPsf == 0) {
      double ra = m_dataMap->mapCenter().ra();
      double dec = m_dataMap->mapCenter().dec();
      s_meanPsf = new MeanPsf(ra, dec, energies, observation);
//      s_meanPsf->write("mean_psf.dat");
   }
   if (s_binnedExposure == 0) {
      s_binnedExposure = new BinnedExposure(energies, observation);
      s_binnedExposure->writeOutput("binned_exposure.fits");
   }
   MeanPsf & psf = *s_meanPsf;
   BinnedExposure & exposure = *s_binnedExposure;

// Loop over source region locations.
   unsigned int indx(0);
   std::vector<double> mu_integrand;
   for (unsigned int i = 0; i < s_mu.size(); i++) {
      double psf_val = psf(energy, s_theta.at(i));
      std::vector<double> phi_integrand;
      for (unsigned int j = 0; j < s_phi.size(); j++) {
         const astro::SkyDir & srcDir = m_srcDirs.at(indx);
         double exposure_val = exposure(energy, srcDir.ra(), srcDir.dec());
         phi_integrand.push_back(exposure_val*psf_val
                                 *m_srcStrengths.at(indx++));
      }
      TrapQuad phiQuad(s_phi, phi_integrand);
      mu_integrand.push_back(phiQuad.integral());
   }
   TrapQuad muQuad(s_mu, mu_integrand);
   double value = muQuad.integral();

   return value;
}

void SourceMap::computeSrcDirs(const Pixel & pixel, Source * src) {
   DiffuseSource * diffuseSrc = dynamic_cast<DiffuseSource *>(src);

// Rotation matrix from Equatorial coords to local coord system
   FitsImage::EquinoxRotation eqRot(pixel.dir().ra(), pixel.dir().dec());

   m_srcDirs.clear();
   m_srcStrengths.clear();
// Loop over source region locations.
   for (unsigned int i = 0; i < s_mu.size(); i++) {
      for (unsigned int j = 0; j < s_phi.size(); j++) {
         astro::SkyDir srcDir;
         getCelestialDir(s_phi[j], s_mu[i], eqRot, srcDir);
         m_srcDirs.push_back(srcDir);
         if (!haveMapCubeFunction(diffuseSrc)) {
            m_srcStrengths.push_back(diffuseSrc->spatialDist(srcDir));
         }
      }
   }
}

void SourceMap::recomputeSrcStrengths(DiffuseSource * src, double energy) {
   m_srcStrengths.clear();
   m_srcStrengths.reserve(m_srcDirs.size());
   std::vector<astro::SkyDir>::const_iterator dir;
   for (dir = m_srcDirs.begin(); dir != m_srcDirs.end(); ++dir) {
      m_srcStrengths.push_back(src->spatialDist(SkyDirArg(*dir, energy)));
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

void SourceMap::getCelestialDir(double phi, double mu, 
                                FitsImage::EquinoxRotation & eqRot,
                                astro::SkyDir & dir) const {
   double sp = sin(phi);
   double arg = mu/sqrt(1 - (1 - mu*mu)*sp*sp);
   double alpha;
   if (cos(phi) < 0) {
      alpha = 2*M_PI - my_acos(arg);
   } else {
      alpha = my_acos(arg);
   }
   double delta = asin(sqrt(1 - mu*mu)*sp);

// The direction in "Equinox rotated" coordinates
   astro::SkyDir indir(alpha*180/M_PI, delta*180/M_PI);

// Convert to the unrotated coordinate system (should probably use 
// Hep3Vector methods here instead).
   eqRot.do_rotation(indir, dir);
}

} // namespace Likelihood
