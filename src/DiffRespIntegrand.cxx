/**
 * @file DiffRespIntegrand.cxx
 *
 * @brief Functor classes that are integrands for the diffuse response
 * calculation. In this implementation the outer integral is over "phi"
 * and the inner is over "mu" for fixed "phi".
 *
 * @author S. Fegan <sfegan@llr.in2p3.fr>, 
 *         J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/DiffRespIntegrand.cxx,v 1.7 2013/01/09 00:44:41 jchiang Exp $
 */

#include <cmath>

#include "st_facilities/GaussianQuadrature.h"

#include "Likelihood/DiffuseSource.h"
#include "Likelihood/EquinoxRotation.h"
#include "Likelihood/Event.h"
#include "Likelihood/ResponseFunctions.h"

#include "Likelihood/DiffRespIntegrand.h"

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

double DiffRespIntegrand::
do2DIntegration(const Event & event,
		const ResponseFunctions & respFuncs,
		const DiffuseSource & src,
		const EquinoxRotation & eqRot,
		double trueEnergy, bool use_edisp,
		double mumin, double mumax, double phimin, double phimax,
		double muerr, double phierr)
{
  DiffRespIntegrand muIntegrand(event, respFuncs, src, eqRot,
				phimin, phimax, phierr);
  int ierr;
  return st_facilities::GaussianQuadrature::dgaus8(muIntegrand, mumin, mumax,
						   muerr, ierr);
}

DiffRespIntegrand::
DiffRespIntegrand(const Event & event,
                  const ResponseFunctions & respFuncs,
                  const DiffuseSource & src,
                  const EquinoxRotation & eqRot, 
		  double trueEnergy, bool use_edisp,
                  double phimin, double phimax, double err)
   : m_event(event), m_respFuncs(respFuncs), m_src(src), m_eqRot(eqRot),
     m_trueEnergy(trueEnergy), m_useEdisp(useEdisp),
     m_phimin(phimin), m_phimax(phimax), m_err(err) {}

double DiffRespIntegrand::operator()(double mu) const {
   DiffRespPhiIntegrand phiIntegrand(mu, *this);

   double err(m_err);
   int ierr;

//    std::cout << mu << "  "
//              << ::my_acos(mu) << "  "
//              << m_event.getEnergy() << std::endl;

   return st_facilities::GaussianQuadrature::dgaus8(phiIntegrand, m_phimin,
                                                    m_phimax, err, ierr);
}

void DiffRespIntegrand::
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

double DiffRespIntegrand::
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

DiffRespIntegrand::DiffRespPhiIntegrand::
DiffRespPhiIntegrand(double mu, const DiffRespIntegrand & muIntegrand) 
   : m_mu(mu), m_muIntegrand(muIntegrand) {}

double DiffRespIntegrand::DiffRespPhiIntegrand::
operator()(double phi) const {
   astro::SkyDir srcDir;
   DiffRespIntegrand::DiffRespPhiIntegrand::getSrcDir(phi, srcDir);

   const Event & event(m_muIntegrand.m_event);
   const ResponseFunctions & respFuncs(m_muIntegrand.m_respFuncs);
   const DiffuseSource & src(m_muIntegrand.m_src);

   double inc(event.getScDir().difference(srcDir)*180./M_PI);
   if (inc > 90) {
      return 0;
   }

   // Neglect energy dispersion.
   double trueEnergy = m_muIntegrand.m_useEdisp ? 
      m_muIntegrand.m_trueEnergy : event.getEnergy();
   double totalResp = 
      respFuncs.totalResponse(trueEnergy, event.getEnergy(), 
                              event.zAxis(), event.xAxis(), 
                              srcDir, event.getDir(), event.getType(),
                              event.getArrTime(),
			      m_phiIntegrand.m_useEdisp);
   double srcDist_val(src.spatialDist(SkyDirArg(srcDir, trueEnergy)));
   
   return totalResp*srcDist_val;
}

void DiffRespIntegrand::DiffRespPhiIntegrand::
getSrcDir(double phi, astro::SkyDir & srcDir) const {
   DiffRespIntegrand::getSrcDir(m_mu, phi, m_muIntegrand.m_eqRot, srcDir);
}

} // namespace Likelihood 
