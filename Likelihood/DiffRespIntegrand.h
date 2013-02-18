/**
 * @file DiffRespIntegrand.h
 * @brief Functor classes that are integrands for the diffuse response
 * calculation.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/DiffRespIntegrand.h,v 1.4 2011/11/12 19:02:04 sfegan Exp $
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

   static double do2DIntegration(const Event & event,
				 const ResponseFunctions & respFuncs,
				 const DiffuseSource & src,
				 const EquinoxRotation & eqRot,
				 double trueEnergy, bool useEdisp,
				 double mumin=-1.0,
				 double mumax=1.0,
				 double phimin=0,
				 double phimax=2*M_PI,
				 double muerr = 0.01,
				 double phierr = 0.1);

   DiffRespIntegrand(const Event & event,
                     const ResponseFunctions & respFuncs,
                     const DiffuseSource & src,
                     const EquinoxRotation & eqRot,
		     double trueEnergy, bool useEdisp,
                     double phimin=0,
                     double phimax=2*M_PI,
		     double err = 0.1);


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

   double m_trueEnergy;
   bool   m_useEdisp;
   double m_phimin;
   double m_phimax;
   double m_err;

};

} // namespace Likelihood

#endif // Likelihood_DiffRespIntegrand_h
