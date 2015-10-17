/** 
 * @file SpatialDisk.h
 * @brief Declaration for the SpatialDisk Function class
 * @author M. Wood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialDisk.h,v 1.27 2015/03/21 05:38:03 jchiang Exp $
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
   double value(const astro::SkyDir &, double energy, const MeanPsf& psf) const;
   double value(double delta, double energy, const MeanPsf& psf) const;

   virtual SpatialDisk * clone() const {
      return new SpatialDisk(*this);
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

#endif // Likelihood_SpatialDisk_h
