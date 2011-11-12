/**
 * @file DiffRespIntegrand2.cxx
 * @brief Implementation for functor classes that are integrands for
 * the diffuse response calculation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/DiffRespIntegrand.cxx,v 1.5 2011/03/18 18:39:11 jchiang Exp $
 */

#include <cmath>

#include "st_facilities/GaussianQuadrature.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/EquinoxRotation.h"
#include "Likelihood/Event.h"
#include "Likelihood/ResponseFunctions.h"

#include "Likelihood/DiffRespIntegrand2.h"

namespace {
   double my_acos(double mu) {
      if (mu > 1) {
         return 0;
      } else if (mu < -1) {
         return M_PI;
      } else {
         return std::acos(mu);
      }
   }
}

namespace Likelihood {

double DiffRespIntegrand2::
do2DIntegration(const Event & event,
		const ResponseFunctions & respFuncs,
		const DiffuseSource & src,
		const EquinoxRotation & eqRot,
		double mumin, double mumax, double phimin, double phimax,
		double muerr, double phierr)
{
  DiffRespIntegrand2 phiIntegrand(event, respFuncs, src, eqRot,
				  mumin, mumax, muerr);
  int ierr;
  return st_facilities::GaussianQuadrature::dgaus8(phiIntegrand, phimin, phimax,
						   phierr, ierr);
}

DiffRespIntegrand2::
DiffRespIntegrand2(const Event & event,
                  const ResponseFunctions & respFuncs,
                  const DiffuseSource & src,
                  const EquinoxRotation & eqRot, 
                  double mumin, double mumax, double err)
   : m_event(event), m_respFuncs(respFuncs), m_src(src), m_eqRot(eqRot),
     m_mumin(mumin), m_mumax(mumax), m_err(err) {}

double DiffRespIntegrand2::operator()(double phi) const {
   DiffRespMuIntegrand muIntegrand(phi, *this);

   double err(m_err);
   int ierr;

   return st_facilities::GaussianQuadrature::dgaus8(muIntegrand, m_mumin,
						    m_mumax, err, ierr);
}

void DiffRespIntegrand2::
getSrcDir(double mu, double phi, const EquinoxRotation & eqRot,
          astro::SkyDir & srcDir) {
   double sp(std::sin(phi));
   double arg(mu/std::sqrt(1. - (1. - mu*mu)*sp*sp));
   double alpha(::my_acos(arg));
   if (std::cos(phi) < 0) {
      alpha = 2*M_PI - alpha;
   }
   double delta(std::asin(std::sqrt(1. - mu*mu)*sp));

   astro::SkyDir indir(alpha*180./M_PI, delta*180./M_PI);

   eqRot.do_rotation(indir, srcDir);
}

double DiffRespIntegrand2::
phiValue(double mu, const astro::SkyDir & srcDir,
         const EquinoxRotation & eqRot) {
   astro::SkyDir refDir;
   getSrcDir(mu, 0, eqRot, refDir);

   double cos_psi = refDir().dot(srcDir());
   double arg = (cos_psi - mu*mu)/(1. - mu*mu);
   double phi0;
   if (arg <= -1) {
      phi0 = M_PI;
   } else if (arg >= 1) {
      phi0 = 0;
   } else {
      phi0 = std::acos(arg);
   }
   double phi1 = -phi0;

   astro::SkyDir srcDir0, srcDir1;
   getSrcDir(mu, phi0, eqRot, srcDir0);
   getSrcDir(mu, phi1, eqRot, srcDir1);
   if (srcDir.difference(srcDir0) < 1e-8) {
      return phi0;
   }
   return phi1;
}

DiffRespIntegrand2::DiffRespMuIntegrand::
DiffRespMuIntegrand(double phi, const DiffRespIntegrand2 & phiIntegrand) 
   : m_phi(phi), m_phiIntegrand(phiIntegrand) {}

double DiffRespIntegrand2::DiffRespMuIntegrand::
operator()(double mu) const {
   astro::SkyDir srcDir;
   DiffRespIntegrand2::DiffRespMuIntegrand::getSrcDir(mu, srcDir);

   const Event & event(m_phiIntegrand.m_event);
   const ResponseFunctions & respFuncs(m_phiIntegrand.m_respFuncs);
   const DiffuseSource & src(m_phiIntegrand.m_src);

   double inc(event.getScDir().difference(srcDir)*180./M_PI);
   if (inc > 90) {
      return 0;
   }

   // Neglect energy dispersion.
   double trueEnergy(event.getEnergy());
   double totalResp = 
      respFuncs.totalResponse(trueEnergy, event.getEnergy(), 
                              event.zAxis(), event.xAxis(), 
                              srcDir, event.getDir(), event.getType());
   double srcDist_val(src.spatialDist(SkyDirArg(srcDir, trueEnergy)));
   
   return totalResp*srcDist_val;
}

void DiffRespIntegrand2::DiffRespMuIntegrand::
getSrcDir(double mu, astro::SkyDir & srcDir) const {
   DiffRespIntegrand2::getSrcDir(mu, m_phi, m_phiIntegrand.m_eqRot, srcDir);
}

} // namespace Likelihood 
