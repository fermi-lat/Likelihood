/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.16 2004/11/04 01:21:12 jchiang Exp $
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
#include "Likelihood/ExposureCube.h"
#include "Likelihood/FitsImage.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/PointSource.h"
#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/TrapQuad.h"

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

SourceMap::SourceMap(Source * src, const CountsMap * dataMap) 
   : m_name(src->getName()), m_dataMap(dataMap), m_deleteDataMap(false) {
   s_refCount++;
   if (s_mu.size() == 0 || s_phi.size() == 0) {
      prepareAngleArrays(100, 50);
   }

   std::vector<Pixel> pixels;
   dataMap->getPixels(pixels);
   
   std::vector<double> energies;
   dataMap->getAxisVector(2, energies);

   std::cerr << "Generating SourceMap for " << m_name;
   long npts = energies.size()*pixels.size();
   m_model.resize(npts, 0);
   long icount(0);

   m_npreds.resize(energies.size(), 0);

   bool havePointSource = dynamic_cast<PointSource *>(src) != 0;
   bool haveDiffuseSource = dynamic_cast<DiffuseSource *>(src) != 0;

   std::vector<Pixel>::const_iterator pixel = pixels.begin();
   for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
      if (haveDiffuseSource) {
         computeSrcDirs(*pixel);
      }
      std::vector<double>::const_iterator energy = energies.begin();
      for (int k = 0; energy != energies.end(); ++energy, k++) {
         unsigned long indx = k*pixels.size() + j;
         if ((icount % (npts/20)) == 0) std::cerr << ".";
         double value(0);
         if (havePointSource) {
/// @todo Ensure the desired event types are correctly included in this
/// calculation.
            for (int evtType = 0; evtType < 2; evtType++) {
               Aeff aeff(src, pixel->dir(), *energy, evtType);
               value += ExposureCube::instance()->value(pixel->dir(), aeff);
            }
         } else if (haveDiffuseSource) {
            value = sourceRegionIntegral(src, *pixel, *energy);
         }
         value *= pixel->solidAngle();
         m_model.at(indx) += value;
         m_npreds.at(k) += value;
         icount++;
      }
   }
   std::cerr << "!" << std::endl;
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
   if (s_mu.size() == 0 || s_phi.size() == 0) {
      prepareAngleArrays(100, 50);
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
   return ResponseFunctions::totalResponse(inclination, phi, m_energy,
                                           m_energy, m_separation, m_type);

}

double SourceMap::sourceRegionIntegral(Source * src, const Pixel & pixel,
                                       double energy) const {
   DiffuseSource * diffuseSrc = dynamic_cast<DiffuseSource *>(src);
   std::vector<double> energies;
   m_dataMap->getAxisVector(2, energies);
   if (s_meanPsf == 0) {
      double ra, dec;
      RoiCuts::instance()->getRaDec(ra, dec);
      s_meanPsf = new MeanPsf(ra, dec, energies);
   }
   if (s_binnedExposure == 0) {
      s_binnedExposure = new BinnedExposure(energies);
   }
   MeanPsf & psf = *s_meanPsf;
   BinnedExposure & exposure = *s_binnedExposure;

// Loop over source region locations.
   unsigned int indx(0);
   std::vector<double> mu_integrand;
   for (unsigned int i = 0; i < s_mu.size(); i++) {
      std::vector<double> phi_integrand;
      for (unsigned int j = 0; j < s_phi.size(); j++) {
         const astro::SkyDir & srcDir = m_srcDirs.at(indx++);
         double separation = srcDir.difference(pixel.dir())*180./M_PI;
         double psf_val = psf(energy, separation);
         double exposure_val = exposure(energy, srcDir.ra(), srcDir.dec());
         double source_strength = diffuseSrc->spatialDist(srcDir);
         phi_integrand.push_back(exposure_val*psf_val*source_strength);
      }
      TrapQuad phiQuad(s_phi, phi_integrand);
      mu_integrand.push_back(phiQuad.integral());
   }
   TrapQuad muQuad(s_mu, mu_integrand);
   double value = muQuad.integral();

   return value;
}

void SourceMap::computeSrcDirs(const Pixel & pixel) {
// Rotation matrix from Equatorial coords to local coord system
   FitsImage::EquinoxRotation eqRot(pixel.dir().ra(), pixel.dir().dec());

   m_srcDirs.clear();
// Loop over source region locations.
   for (unsigned int i = 0; i < s_mu.size(); i++) {
      for (unsigned int j = 0; j < s_phi.size(); j++) {
         astro::SkyDir srcDir;
         getCelestialDir(s_phi[j], s_mu[i], eqRot, srcDir);
         m_srcDirs.push_back(srcDir);
      }
   }
}

void SourceMap::prepareAngleArrays(int nmu, int nphi) {
   double radius = RoiCuts::instance()->extractionRegion().radius()
      *sqrt(2.) + 10.;

   double mumin = cos(radius*M_PI/180);

// Sample more densely near theta = 0:
   std::deque<double> my_mu;
   double nscale = static_cast<double>((nmu-1)*(nmu-1));
   for (int i = 0; i < nmu; i++) {
      my_mu.push_front(1. - i*i/nscale*(1. - mumin));
   }
   s_mu.resize(my_mu.size());
   std::copy(my_mu.begin(), my_mu.end(), s_mu.begin());

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
