/**
 * @file MapCubeFunction.h
 * @brief Encapsulation of a 3D FITS image.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef Likelihood_MapCubeFunction_h
#define Likelihood_MapCubeFunction_h

#include "optimizers/Function.h"

namespace Likelihood {

   class SkyDirArg;

/**
 * @class MapCubeFunction
 * @brief Encapsulation of a 3D FITS image, e.g., with RA, Dec, energy
 * dimensions, for modeling diffuse sources with spectral variation as
 * a function of position on the sky.
 * @author J. Chiang
 *
 * $Header$
 */

class MapCubeFunction : public optimizers::Function {

public:
   
   MapCubeFunction(const std::string & fitsFile);

   virtual ~MapCubeFunction() {}

   virtual double value(optimizers::Arg &) const;

   virtual double derivByParam(optimizers::Arg & dir,
                               const std::string & paramName) const {
// There is only the normalization, so the derivative is easy:
      return value(dir)/getParamValue(paramName);
   }

private:

   std::string m_coordSys;
   std::vector<double> m_image;
   std::vector<double> m_lon;
   double m_lonMin, m_lonMax;
   std::vector<double> m_lat;
   double m_latMin, m_latMax;
   std::vector<double> m_energies;

   void readEnergyVector(const std::string & fitsFile);

   static int findIndex(std::vector<double> xx, double x);

};

} // namespace Likelihood

#endif // Likelihood_MapCubeFunction_h
