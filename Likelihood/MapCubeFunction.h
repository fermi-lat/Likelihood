/**
 * @file MapCubeFunction.h
 * @brief Encapsulation of a 3D FITS image.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MapCubeFunction.h,v 1.8 2009/02/17 06:48:47 jchiang Exp $
 */

#ifndef Likelihood_MapCubeFunction_h
#define Likelihood_MapCubeFunction_h

#include "optimizers/Function.h"

#include "Likelihood/MapBase.h"

namespace astro {
   class SkyProj;
}

namespace Likelihood {

/**
 * @class MapCubeFunction
 * @brief Encapsulation of a 3D FITS image, e.g., with RA, Dec, energy
 * dimensions, for modeling diffuse sources with spectral variation as
 * a function of position on the sky.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MapCubeFunction.h,v 1.8 2009/02/17 06:48:47 jchiang Exp $
 */

class MapCubeFunction : public optimizers::Function, public MapBase {

public:
   
   MapCubeFunction();

   MapCubeFunction(const std::string & fitsFile);

   MapCubeFunction(const MapCubeFunction &);

   virtual ~MapCubeFunction();

   virtual MapCubeFunction & operator=(const MapCubeFunction &);

   virtual double value(optimizers::Arg &) const;

   virtual double derivByParam(optimizers::Arg & dir,
                               const std::string & paramName) const {
// There is only the normalization, so the derivative is easy:
      return value(dir)/getParamValue(paramName);
   }

   virtual MapCubeFunction * clone() const {
      return new MapCubeFunction(*this);
   }

   virtual void readFitsFile(const std::string & fitsFile,
                             const std::string & extension="");

   const std::string & fitsFile() const {
      return m_fitsFile;
   }

   /// @todo Consider moving this method to FitsImage class.
   double mapIntegral() const;

private:

   std::string m_fitsFile;

   astro::SkyProj * m_proj;

   int m_nlon;
   int m_nlat;

   bool m_isPeriodic;

   std::vector<double> m_energies;

   std::vector<float> m_image;

   void init();

   int findIndex(const std::vector<double> & xx, double x) const;

   double powerLawIntegral(double x1, double x2, double y1, double y2) const;
};

} // namespace Likelihood

#endif // Likelihood_MapCubeFunction_h
