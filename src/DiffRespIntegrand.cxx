/**
 * @file DiffRespIntegrand.cxx
 * @brief Implementation for functor classes that are integrands for
 * the diffuse response calculation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/DiffRespIntegrand.cxx,v 1.2 2007/12/11 07:14:49 jchiang Exp $
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

DiffRespIntegrand::
DiffRespIntegrand(const Event & event,
                  const ResponseFunctions & respFuncs,
                  const DiffuseSource & src,
                  const EquinoxRotation & eqRot)
   : m_event(event), m_respFuncs(respFuncs), m_src(src), m_eqRot(eqRot) {}

double DiffRespIntegrand::operator()(double mu) const {
   DiffRespPhiIntegrand phiIntegrand(mu, *this);

   double phimin(0);
   double phimax(2*M_PI);
   double err(1e-1);
   int ierr;

//    std::cout << mu << "  "
//              << ::my_acos(mu) << "  "
//              << m_event.getEnergy() << std::endl;

   return st_facilities::GaussianQuadrature::dgaus8(phiIntegrand, phimin,
                                                    phimax, err, ierr);
}

DiffRespIntegrand::DiffRespPhiIntegrand::
DiffRespPhiIntegrand(double mu, const DiffRespIntegrand & muIntegrand) 
   : m_mu(mu), m_muIntegrand(muIntegrand) {}

double DiffRespIntegrand::DiffRespPhiIntegrand::
operator()(double phi) const {
   astro::SkyDir srcDir;
   getSrcDir(phi, srcDir);

   const Event & event(m_muIntegrand.m_event);
   const ResponseFunctions & respFuncs(m_muIntegrand.m_respFuncs);
   const DiffuseSource & src(m_muIntegrand.m_src);

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

void DiffRespIntegrand::DiffRespPhiIntegrand::
getSrcDir(double phi, astro::SkyDir & srcDir) const {
   double sp(std::sin(phi));
   double arg(m_mu/std::sqrt(1. - (1. - m_mu*m_mu)*sp*sp));
   double alpha(::my_acos(arg));
   if (std::cos(phi) < 0) {
      alpha = 2*M_PI - alpha;
   }
   double delta(std::asin(std::sqrt(1. - m_mu*m_mu)*sp));

   astro::SkyDir indir(alpha*180./M_PI, delta*180./M_PI);

   m_muIntegrand.m_eqRot.do_rotation(indir, srcDir);
}

} // namespace Likelihood 
