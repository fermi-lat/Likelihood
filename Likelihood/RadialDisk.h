/** 
 * @file RadialDisk.h
 * @brief Declaration for the RadialDisk Function class
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/RadialDisk.h,v 1.3 2016/01/27 02:48:07 mdwood Exp $
 *
 */

#ifndef Likelihood_RadialDisk_h
#define Likelihood_RadialDisk_h

#include <utility>

#include "optimizers/Function.h"

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
 * @class RadialDisk
 *
 * @brief A Function object that is an analytic representation of a 2D
 * spatial disk at a SkyDir location with functional form:
 *
 * f(r,R) = 1/(pi * R^2) * H(1-r/R)
 *
 * where r is the angular distance from the SkyDir location and H(x)
 * is the Heaviside step function.
 *
 */
    
class RadialDisk : public SpatialFunction {

public:

   RadialDisk();

   RadialDisk(double ra, double dec, double radius);
                         
   RadialDisk(const RadialDisk &);

   RadialDisk & operator=(const RadialDisk &);

   virtual ~RadialDisk();

   double value(const astro::SkyDir &) const;
   double value(double delta, double radius) const;

   double spatialResponse(const astro::SkyDir &, double energy, const MeanPsf& psf) const;
   double spatialResponse(double delta, double energy, const MeanPsf& psf) const;

   virtual double diffuseResponse(const ResponseFunctor& fn, double energy,
				  double separation) const;

   virtual double getDiffRespLimits(const astro::SkyDir &, 
				    double & mumin, double & mumax,
				    double & phimin, double & phimax) const;

   virtual RadialDisk * clone() const {
      return new RadialDisk(*this);
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

   double         m_radius;
};

} // namespace Likelihood

#endif // Likelihood_RadialDisk_h
