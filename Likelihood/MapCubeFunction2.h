/**
 * @file MapCubeFunction2.h
 * @brief Encapsulation of a 3D FITS image.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MapCubeFunction2.h,v 1.4 2012/01/06 07:11:58 jchiang Exp $
 */

#ifndef Likelihood_MapCubeFunction2_h
#define Likelihood_MapCubeFunction2_h

#include "optimizers/Function.h"

#include "Likelihood/MapBase.h"

namespace astro {
   class SkyProj;
}

namespace Likelihood {

   class WcsMap2;

/**
 * @class MapCubeFunction2
 *
 * @brief Encapsulation of a 3D FITS image, e.g., with RA, Dec, energy
 * dimensions, for modeling diffuse sources with spectral variation as
 * a function of position on the sky.
 *
 */

class MapCubeFunction2 : public optimizers::Function, public MapBase {

public:
   
   MapCubeFunction2();

   MapCubeFunction2(const std::string & fitsFile);

   MapCubeFunction2(const MapCubeFunction2 &);

   virtual ~MapCubeFunction2();

   virtual MapCubeFunction2 & operator=(const MapCubeFunction2 &);

   virtual MapCubeFunction2 * clone() const {
      return new MapCubeFunction2(*this);
   }

   double mapIntegral() const;

   /// @return The angular integral of the differential flux as a
   /// function of energy (implementation for
   /// MapBase::mapIntegral(double) const)
   virtual double mapIntegral(double energy) const;

   virtual void integrateSpatialDist(const std::vector<double> & energies,
                                     const ExposureMap & expmap,
                                     std::vector<double> & exposure) const;
protected:

   virtual double value(optimizers::Arg &) const;

   virtual double derivByParamImp(optimizers::Arg & dir,
                                  const std::string & paramName) const {
      // There is only the normalization, so the derivative is easy:
      return value(dir)/getParamValue(paramName);
   }

};

} // namespace Likelihood

#endif // Likelihood_MapCubeFunction2_h
