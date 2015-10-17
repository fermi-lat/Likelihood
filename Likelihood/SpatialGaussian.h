/** 
 * @file SpatialGaussian.h
 * @brief Declaration for the SpatialGaussian Function class
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialGaussian.h,v 1.27 2015/03/21 05:38:03 jchiang Exp $
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
 * spatial gaussian at a SkyDir location.
 *
 */
    
class SpatialGaussian : public SpatialFunction {

public:

   SpatialGaussian();

   SpatialGaussian(double ra, double dec, double width);
                         
   SpatialGaussian(const SpatialGaussian &);

   SpatialGaussian & operator=(const SpatialGaussian &);

   virtual ~SpatialGaussian();

   double value(const astro::SkyDir &) const;
   double value(double delta, double width) const;
   double value(const astro::SkyDir &, double energy, const MeanPsf& psf) const;
   double value(double delta, double energy, const MeanPsf& psf) const;

   virtual SpatialGaussian * clone() const {
      return new SpatialGaussian(*this);
   }

   virtual void update();

   static double convolve(const MeanPsf& psf, double energy, double x,
			  double sigma, int n = 100);

protected:

   double value(const optimizers::Arg &) const;

   double derivByParamImp(const optimizers::Arg &, const std::string &) const;

private:

   double         m_width;
};

} // namespace Likelihood

#endif // Likelihood_SpatialGaussian_h
