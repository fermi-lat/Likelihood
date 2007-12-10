/**
 * @file DiffRespIntegrand.h
 * @brief Functor classes that are integrands for the diffuse response
 * calculation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef Likelihood_DiffRespIntegrand_h
#define Likelihood_DiffRespIntegrand_h

namespace astro {
   class SkyDir;
}

#include "Likelihood/EquinoxRotation.h"

namespace Likelihood {

class DiffuseSource;
class Event;
class ResponseFunctions;

/**
 * @class DiffRespIntegrand
 */

class DiffRespIntegrand {

public:

   DiffRespIntegrand(const Event & event,
                     const ResponseFunctions & respFuncs,
                     const DiffuseSource & src);

   double operator()(double mu) const;

private:

   /**
    * @class DiffRespPhiIntegrand
    *
    * @brief Nested class for phi-integrand
    *
    */

   class DiffRespPhiIntegrand {
      
   public:

      DiffRespPhiIntegrand(double mu, const DiffRespIntegrand & muIntegrand);

      double operator()(double phi) const;

   private:

      double m_mu;

      const DiffRespIntegrand & m_muIntegrand;

      void getSrcDir(double phi, astro::SkyDir & srcDir) const;

   };

   const Event & m_event;
   const ResponseFunctions & m_respFuncs;
   const DiffuseSource & m_src;

   EquinoxRotation m_eqRot;

};

} // namespace Likelihood

#endif // Likelihood_DiffRespIntegrand_h
