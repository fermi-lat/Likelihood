/** 
 * @file Event.cxx
 * @brief Event class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Event.cxx,v 1.11 2003/05/06 23:47:55 jchiang Exp $
 */

#include "Likelihood/Event.h"
#include "Likelihood/Psf.h"
#include "Likelihood/Aeff.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/TrapQuad.h"
#include <assert.h>

extern double my_acos(double mu);

double my_acos(double mu) {
   if (mu > 1) {
      return 0;
   } else if (mu < -1) {
      return M_PI;
   } else {
      return acos(mu);
   }
}   

namespace Likelihood {

Event::Event(double ra, double dec, double energy, 
             double time, double sc_ra, double sc_dec, double muZenith) {
   
   m_appDir = astro::SkyDir(ra, dec);
   m_energy = energy;
   m_arrTime = time;
   m_scDir = astro::SkyDir(sc_ra, sc_dec);
   m_muZenith = muZenith;
}

Event::Event(const Event &event) {
   m_appDir = event.m_appDir;
   m_energy = event.m_energy;
   m_arrTime = event.m_arrTime;
   m_scDir = event.m_scDir;
   m_respEg = event.m_respEg;
   m_respGal = event.m_respGal;
   m_respDiffuseSrcs = event.m_respDiffuseSrcs;
}

double Event::diffuseResponse(double energy, 
                              const std::string &diffuseComponent) const {
// Since the energy resolution is presently assumed to be infinite,
// simply return the (second member of the pair of the) first (and
// only) element of the diffuse_response vector.

   std::map<std::string, diffuse_response>::const_iterator it;
   if ((it = m_respDiffuseSrcs.find(diffuseComponent))
       != m_respDiffuseSrcs.end()) {
      return it->second[0].second;
   } else {
      std::cerr << "Event::diffuseResponse: Diffuse component " 
                << diffuseComponent 
                << " does not have an associated diffuse response."
                << std::endl;
      assert(false);
      return 0;
   }
   return 0;
}

void Event::computeResponse(DiffuseSource &src, double sr_radius) {
   Psf *psf = Psf::instance();
   Aeff *aeff = Aeff::instance();

// Create the EquinoxRotation object for this Event.
//   FitsImage::EquinoxRotation eqRot(m_appDir.ra(), m_appDir.dec());
// The following instantiates the eqRot object to do the same rotation
// as performed by my evt_exposr FTOOL
   RoiCuts *roi_cuts = RoiCuts::instance();
   std::pair<astro::SkyDir, double> roi = roi_cuts->getExtractionRegion();
   astro::SkyDir roiCenter = roi.first;
   FitsImage::EquinoxRotation eqRot(roiCenter.ra(), roiCenter.dec());

// Do this calculation assuming infinite energy resolution and thus
// compute a single value.

   int nmu = 70;
   double mumin = cos(sr_radius*M_PI/180);
   double mustep = (1. - mumin)/(nmu - 1.);
   std::vector<double> mu;
   for (int i = 0; i < nmu; i++) mu.push_back(mustep*i + mumin);

   int nphi = 40;
   double phistep = 2.*M_PI/(nphi - 1.);
   std::vector<double> phi;
   for (int i = 0; i < nphi; i++) phi.push_back(phistep*i);

   std::vector<double> mu_integrand;
   for (int i = 0; i < nmu; i++) {
      std::vector<double> phi_integrand;
      for (int j = 0; j < nphi; j++) {
         astro::SkyDir srcDir;
         getCelestialDir(phi[j], mu[i], eqRot, srcDir);
// !!!need to canonicalize angular units for inputs to Psf!!!
         double inc = m_scDir.SkyDir::difference(srcDir)*180/M_PI;
         if (inc < Response::incMax()) {
            double separation = m_appDir.SkyDir::difference(srcDir);
            double psf_val = (*psf)(separation, m_energy, inc);
            double aeff_val = (*aeff)(m_energy, inc);
            double srcDist_val = src.spatialDist(srcDir);
            phi_integrand.push_back(psf_val*aeff_val*srcDist_val);
         } else {
            phi_integrand.push_back(0);
         }
      }
      TrapQuad phiQuad(phi, phi_integrand);
      mu_integrand.push_back(phiQuad.integral());
   }
   TrapQuad muQuad(mu, mu_integrand);
   diffuse_response diff_resp;
   diff_resp.push_back(std::make_pair(m_energy, muQuad.integral()));

//     std::cout << diff_resp[0].first << "  "
//               << diff_resp[0].second << std::endl;

   m_respDiffuseSrcs[src.getName()] = diff_resp;
}

void Event::computeResponse(std::vector<DiffuseSource> &srcs, 
                            double sr_radius) {
   Psf *psf = Psf::instance();
   Aeff *aeff = Aeff::instance();

// Create the EquinoxRotation object for this Event.
//   FitsImage::EquinoxRotation eqRot(m_appDir.ra(), m_appDir.dec());
// The following instantiates the eqRot object to do the same rotation
// as performed by my evt_exposr FTOOL
   RoiCuts *roi_cuts = RoiCuts::instance();
   std::pair<astro::SkyDir, double> roi = roi_cuts->getExtractionRegion();
   astro::SkyDir roiCenter = roi.first;
   FitsImage::EquinoxRotation eqRot(roiCenter.ra(), roiCenter.dec());

// Do this calculation assuming infinite energy resolution and thus
// compute a single value.

   int nmu = 70;
   double mumin = cos(sr_radius*M_PI/180);
   double mustep = (1. - mumin)/(nmu - 1.);
   std::vector<double> mu;
   for (int i = 0; i < nmu; i++) mu.push_back(mustep*i + mumin);

   int nphi = 40;
   double phistep = 2.*M_PI/(nphi - 1.);
   std::vector<double> phi;
   for (int i = 0; i < nphi; i++) phi.push_back(phistep*i);

   std::vector< std::vector<double> > mu_integrands;
   mu_integrands.resize(srcs.size());
   for (int i = 0; i < nmu; i++) {
      std::vector< std::vector<double> > phi_integrands;
      phi_integrands.resize(srcs.size());
      for (int j = 0; j < nphi; j++) {
         astro::SkyDir srcDir;
         getCelestialDir(phi[j], mu[i], eqRot, srcDir);
// !!!need to canonicalize angular units for inputs to Psf!!!
         double inc = m_scDir.SkyDir::difference(srcDir)*180/M_PI;
         if (inc < Response::incMax()) {
            double separation = m_appDir.SkyDir::difference(srcDir);
            double psf_val = (*psf)(separation, m_energy, inc);
            double aeff_val = (*aeff)(m_energy, inc);
            for (unsigned int k = 0; k < srcs.size(); k++) {
               double srcDist_val = srcs[k].spatialDist(srcDir);
               phi_integrands[k].push_back(psf_val*aeff_val*srcDist_val);
            }
         } else {
            for (unsigned int k = 0; k < srcs.size(); k++)
               phi_integrands[k].push_back(0);
         }
      }
      for (unsigned int k = 0; k < srcs.size(); k++) {
         TrapQuad phiQuad(phi, phi_integrands[k]);
         mu_integrands[k].push_back(phiQuad.integral());
      }
   }
   for (unsigned int k = 0; k < srcs.size(); k++) {
      TrapQuad muQuad(mu, mu_integrands[k]);
      diffuse_response diff_resp;
      diff_resp.push_back(std::make_pair(m_energy, muQuad.integral()));

      m_respDiffuseSrcs[srcs[k].getName()] = diff_resp;
   }
}

void Event::getCelestialDir(double phi, double mu, 
                            FitsImage::EquinoxRotation &eqRot,
                            astro::SkyDir &dir) {
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
