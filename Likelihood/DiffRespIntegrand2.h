/**
 * @file DiffRespIntegrand2.h
 * @brief Functor classes that are integrands for the diffuse response
 * calculation. In this implementation the outer integral is over "phi"
 * and the inner is over "mu" for fixed "phi".
 *
 * @author S. Fegan <sfegan@llr.in2p3.fr>, 
 *         J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/DiffRespIntegrand2.h,v 1.1 2011/11/12 19:02:04 sfegan Exp $
 */

#ifndef Likelihood_DiffRespIntegrand2_h
#define Likelihood_DiffRespIntegrand2_h

namespace astro {
   class SkyDir;
}

namespace Likelihood {

class DiffuseSource;
class EquinoxRotation;
class Event;
class ResponseFunctions;

/**
 * @class DiffRespIntegrand2
 */

class DiffRespIntegrand2 {

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

   DiffRespIntegrand2(const Event & event,
		      const ResponseFunctions & respFuncs,
		      const DiffuseSource & src,
		      const EquinoxRotation & eqRot,
		      double trueEnergy, bool useEdisp,
		      double mumin=-1.0,
		      double mumax=1.0,
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
   class DiffRespMuIntegrand {
      
   public:

      DiffRespMuIntegrand(double phi, const DiffRespIntegrand2 & muIntegrand);

      double operator()(double mu) const;

   private:

      double m_phi;

      const DiffRespIntegrand2 & m_phiIntegrand;

      void getSrcDir(double mu, astro::SkyDir & srcDir) const;

   };
#endif

private:

   const Event & m_event;
   const ResponseFunctions & m_respFuncs;
   const DiffuseSource & m_src;
   const EquinoxRotation & m_eqRot;

   double m_trueEnergy;
   bool   m_useEdisp;
   double m_mumin;
   double m_mumax;
   double m_err;

};

} // namespace Likelihood

#endif // Likelihood_DiffRespIntegrand2_h
