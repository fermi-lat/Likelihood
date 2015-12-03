/** 
 * @file SpatialDisk.h
 * @brief Declaration for the SpatialDisk Function class
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialDisk.h,v 1.1 2015/10/17 17:19:14 mdwood Exp $
 *
 */

#ifndef Likelihood_SpatialDisk_h
#define Likelihood_SpatialDisk_h

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
 * @class SpatialDisk
 *
 * @brief A Function object that is an analytic representation of a 2D
 * spatial disk at a SkyDir location.
 *
 */
    
class SpatialDisk : public SpatialFunction {

public:

   SpatialDisk();

   SpatialDisk(double ra, double dec, double width);
                         
   SpatialDisk(const SpatialDisk &);

   SpatialDisk & operator=(const SpatialDisk &);

   virtual ~SpatialDisk();

   double value(const astro::SkyDir &) const;
   double value(double delta, double width) const;

   double spatialResponse(const astro::SkyDir &, double energy, const MeanPsf& psf) const;
   double spatialResponse(double delta, double energy, const MeanPsf& psf) const;

   virtual double diffuseResponse(const ResponseFunctor& fn, double energy,
				  double separation) const;

   virtual double getDiffRespLimits(const astro::SkyDir &, 
				    double & mumin, double & mumax,
				    double & phimin, double & phimax) const;

   virtual SpatialDisk * clone() const {
      return new SpatialDisk(*this);
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

   double         m_width;
};

} // namespace Likelihood

#endif // Likelihood_SpatialDisk_h
