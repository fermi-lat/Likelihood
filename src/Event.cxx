/** 
 * @file Event.cxx
 * @brief Event class implementation
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/Event.cxx,v 1.21 2003/12/04 04:15:29 jchiang Exp $
 */

#include <cassert>
#include <iostream>
#include <sstream>
#include <utility>

#include "latResponse/IPsf.h"
#include "latResponse/IAeff.h"
#include "latResponse/Irfs.h"
#include "latResponse/../src/Glast25.h"

#include "Likelihood/ResponseFunctions.h"
#include "Likelihood/Event.h"
#include "Likelihood/RoiCuts.h"
#include "Likelihood/ScData.h"
#include "Likelihood/DiffuseSource.h"
#include "Likelihood/TrapQuad.h"
#include "Likelihood/Exception.h"

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
   double totalResponse(double energy, double time, 
                        const astro::SkyDir &srcDir,
                        const astro::SkyDir &appDir, int type) {
// This implementation neglects energy dispersion.
      Likelihood::ResponseFunctions * respFuncs 
         = Likelihood::ResponseFunctions::instance();
      Likelihood::ScData * scData = Likelihood::ScData::instance();
   
      astro::SkyDir zAxis = scData->zAxis(time);
      astro::SkyDir xAxis = scData->xAxis(time);

      double myResponse = 0;
      std::map<unsigned int, latResponse::Irfs *>::iterator respIt
         = respFuncs->begin();
      for ( ; respIt != respFuncs->end(); respIt++) {
         if (respIt->second->irfID() == type) {
            latResponse::IPsf *psf = respIt->second->psf();
            latResponse::IAeff *aeff = respIt->second->aeff();
            double psf_val = psf->value(appDir, energy, srcDir, zAxis, xAxis);
            double aeff_val = aeff->value(energy, srcDir, zAxis, xAxis);
            myResponse += psf_val*aeff_val;
         }
      }
      return myResponse;
   }

}

namespace Likelihood {

Event::Event(double ra, double dec, double energy, 
             double time, double sc_ra, double sc_dec, 
             double muZenith, int type) {
   m_appDir = astro::SkyDir(ra, dec);
   m_energy = energy;
   m_arrTime = time;
   m_scDir = astro::SkyDir(sc_ra, sc_dec);
   m_muZenith = muZenith;
   m_type = type;
}

Event::Event(const Event &event) {
   m_appDir = event.m_appDir;
   m_energy = event.m_energy;
   m_arrTime = event.m_arrTime;
   m_scDir = event.m_scDir;
   m_muZenith = event.m_muZenith;
   m_type = event.m_type;
   m_respEg = event.m_respEg;
   m_respGal = event.m_respGal;
   m_respDiffuseSrcs = event.m_respDiffuseSrcs;
}

// double Event::diffuseResponse(double energy, 
//                               const std::string &diffuseComponent) const 
double Event::diffuseResponse(double,
                              const std::string &diffuseComponent) const 
   throw(Exception) {
// Since the energy resolution is presently assumed to be infinite,
// simply return the (second member of the pair of the) first (and
// only) element of the diffuse_response vector.

   std::map<std::string, diffuse_response>::const_iterator it;
   if ((it = m_respDiffuseSrcs.find(diffuseComponent))
       != m_respDiffuseSrcs.end()) {
      return it->second[0].second;
   } else {
      std::ostringstream errorMessage;
      errorMessage << "Event::diffuseResponse: \nDiffuse component " 
                   << diffuseComponent 
                   << " does not have an associated diffuse response.\n";
      throw Exception(errorMessage.str());
   }
   return 0;
}

void Event::computeResponse(std::vector<DiffuseSource *> &srcs, 
                            double sr_radius) {

// Create the EquinoxRotation object for this Event.
//   FitsImage::EquinoxRotation eqRot(m_appDir.ra(), m_appDir.dec());
// The following instantiates the eqRot object to do the same rotation
// as performed by my evt_exposr FTOOL
   RoiCuts *roi_cuts = RoiCuts::instance();
   astro::SkyDir roiCenter = roi_cuts->extractionRegion().center();
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
         double inc = m_scDir.SkyDir::difference(srcDir)*180/M_PI;
         if (inc < latResponse::Glast25::incMax()) {
            double totalResp = ::totalResponse(m_energy, m_arrTime,
                                               srcDir, m_appDir, m_type);
            for (unsigned int k = 0; k < srcs.size(); k++) {
               double srcDist_val = srcs[k]->spatialDist(srcDir);
               phi_integrands[k].push_back(totalResp*srcDist_val);
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

      m_respDiffuseSrcs[srcs[k]->getName()] = diff_resp;
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
