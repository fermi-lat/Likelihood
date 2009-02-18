/**
 * @file DiffRespIntegrand.h
 * @brief Functor classes that are integrands for the diffuse response
 * calculation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DiffRespIntegrand.h,v 1.2 2007/12/11 07:14:48 jchiang Exp $
 */

#ifndef Likelihood_DiffRespIntegrand_h
#define Likelihood_DiffRespIntegrand_h

namespace astro {
   class SkyDir;
}

namespace Likelihood {

class DiffuseSource;
class EquinoxRotation;
class Event;
class ResponseFunctions;

/**
 * @class DiffRespIntegrand
 */

class DiffRespIntegrand {

public:

   DiffRespIntegrand(const Event & event,
                     const ResponseFunctions & respFuncs,
                     const DiffuseSource & src,
                     const EquinoxRotation & eqRot,
                     double phimin=0,
                     double phimax=2*M_PI);

   double operator()(double mu) const;

   static void getSrcDir(double mu, double phi, const EquinoxRotation & eqRot,
                         astro::SkyDir & srcDir);

   static double phiValue(double mu, const astro::SkyDir & srcDir,
                          const EquinoxRotation & eqRot);

#ifndef SWIG
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
#endif

private:

   const Event & m_event;
   const ResponseFunctions & m_respFuncs;
   const DiffuseSource & m_src;
   const EquinoxRotation & m_eqRot;

   double m_phimin;
   double m_phimax;

};

} // namespace Likelihood

#endif // Likelihood_DiffRespIntegrand_h
