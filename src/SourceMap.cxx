/**
 * @file SourceMap.cxx
 * @brief Spatial distribution of a source folded through the instrument
 *        response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/SourceMap.cxx,v 1.10 2004/10/07 00:01:05 jchiang Exp $
 */

#include <algorithm>
#include <deque>
#include <memory>

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

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

MeanPsf * SourceMap::s_meanPsf(0);
BinnedExposure * SourceMap::s_binnedExposure(0);
std::vector<double> SourceMap::s_phi;
std::vector<double> SourceMap::s_mu;

SourceMap::SourceMap(Source * src, const CountsMap & dataMap) 
   : m_name(src->getName()), m_dataMap(dataMap) {
   if (s_mu.size() == 0 || s_phi.size() == 0) {
      prepareAngleArrays(100, 50);
   }

   std::vector<Pixel> pixels;
   dataMap.getPixels(pixels);
   
   std::vector<double> energies;
   dataMap.getAxisVector(2, energies);

   std::cerr << "Generating SourceMap for " << m_name;
   long npts = energies.size()*pixels.size();
   m_model.resize(npts, 0);
   long icount(0);

   m_npreds.resize(energies.size(), 0);

   std::vector<double>::const_iterator energy = energies.begin();
   for (int k = 0; energy != energies.end(); ++energy, k++) {
      std::vector<Pixel>::const_iterator pixel = pixels.begin();
      for (int j = 0; pixel != pixels.end(); ++pixel, j++) {
         unsigned long indx = k*pixels.size() + j;
         if ((icount % (npts/10)) == 0) std::cerr << ".";
         double value(0);
         if (dynamic_cast<PointSource *>(src) != 0) {
/// @todo Ensure the desired event types are correctly included in this
/// calculation.
            for (int evtType = 0; evtType < 2; evtType++) {
               Aeff aeff(src, pixel->dir(), *energy, evtType);
               value += ExposureCube::instance()->value(pixel->dir(), aeff)
                  *pixel->solidAngle();
            }
         } else if ( dynamic_cast<DiffuseSource *>(src) != 0 ) {
            value = sourceRegionIntegral(src, *pixel, *energy);
         }
         m_model.at(indx) += value;
         m_npreds.at(k) += value;
         icount++;
      }
   }
   std::cerr << "!" << std::endl;
}

SourceMap::SourceMap(const std::string & sourceMapsFile,
                     const std::string & srcName) 
   : m_name(srcName), m_dataMap(CountsMap(sourceMapsFile)) {
   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(sourceMapsFile, srcName));
   std::vector<float> image_data;
   image->get(image_data);
   m_model.resize(image_data.size());
   std::copy(image_data.begin(), image_data.end(), m_model.begin());

   std::vector<Pixel> pixels;
   m_dataMap.getPixels(pixels);
   std::vector<double> energies;
   m_dataMap.getAxisVector(2, energies);

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
   m_dataMap.getAxisVector(2, energies);
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

// Rotation matrix from Equatorial coords to local coord system
   FitsImage::EquinoxRotation eqRot(pixel.dir().ra(), pixel.dir().dec());

// Loop over source region locations.
   std::vector<double> mu_integrand;
   for (unsigned int i = 0; i < s_mu.size(); i++) {
      std::vector<double> phi_integrand;
      for (unsigned int j = 0; j < s_phi.size(); j++) {
         astro::SkyDir srcDir;
         getCelestialDir(s_phi[j], s_mu[i], eqRot, srcDir);
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
