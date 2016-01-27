/** 
 * @file SpatialGaussian.h
 * @brief Declaration for the SpatialGaussian Function class
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/SpatialGaussian.h,v 1.2 2015/12/03 23:47:44 mdwood Exp $
 *
 */

#ifndef Likelihood_SpatialGaussian_h
#define Likelihood_SpatialGaussian_h

#include <utility>

#include "optimizers/Function.h"
#include "optimizers/Parameter.h"

#include "Likelihood/MapBase.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/SpatialFunction.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

class Event;
class ResponseFunctions;

/** 
 * @class SpatialGaussian
 *
 * @brief A Function object that is an analytic representation of a 2D
 * spatial gaussian at a SkyDir location with functional form:
 *
 * f(r, sigma) = 1/(2 * pi * sigma^2) * exp(-r^2/(2 * sigma^2))
 *
 * where r is the angular distance from the SkyDir location.
 *
 */
    
class SpatialGaussian : public SpatialFunction {

public:

   SpatialGaussian();

   SpatialGaussian(double ra, double dec, double sigma);
                         
   SpatialGaussian(const SpatialGaussian &);

   SpatialGaussian & operator=(const SpatialGaussian &);

   virtual ~SpatialGaussian();

   double value(const astro::SkyDir &) const;
   double value(double delta, double sigma) const;

   double spatialResponse(const astro::SkyDir &, double energy, const MeanPsf& psf) const;
   double spatialResponse(double delta, double energy, const MeanPsf& psf) const;

   virtual double diffuseResponse(const ResponseFunctor& fn, double energy,
				  double dtheta) const;

   virtual double getDiffRespLimits(const astro::SkyDir &, 
				    double & mumin, double & mumax,
				    double & phimin, double & phimax) const;

   virtual SpatialGaussian * clone() const {
      return new SpatialGaussian(*this);
   }

   virtual void update();

   static double convolve(const ResponseFunctor& fn, double energy, double separation,
			  double sigma, double tol = 0.0001);

#ifndef SWIG
   /**
    * @class RadialIntegrand
    *
    * @brief Integrand for radial part of convolution.
    *
    */
   class RadialIntegrand {
      
   public:

   RadialIntegrand(const ResponseFunctor& fn, double energy, double x, double sigma): 
     m_fn(fn), m_energy(energy), m_x(x), m_sigma(sigma) { }

     double operator()(double x) const;

     static int s_ncall;
   private:
     const Likelihood::ResponseFunctor& m_fn;
     double m_energy;
     double m_x;
     double m_sigma;
   };
#endif

protected:

   double value(const optimizers::Arg &) const;

   double derivByParamImp(const optimizers::Arg &, const std::string &) const;

private:

   double         m_sigma;
};

} // namespace Likelihood

#endif // Likelihood_SpatialGaussian_h
