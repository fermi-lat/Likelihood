/**
 * @file MapCubeFunction.h
 * @brief Encapsulation of a 3D FITS image.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MapCubeFunction.h,v 1.3 2005/02/15 07:04:41 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/MapCubeFunction.h,v 1.3 2005/02/15 07:04:41 jchiang Exp $
 */

class MapCubeFunction : public optimizers::Function {

public:
   
   MapCubeFunction() {
      init();
   }

   MapCubeFunction(const std::string & fitsFile) {
      init();
      readFitsFile(fitsFile);
   }

   virtual ~MapCubeFunction() {}

   virtual double value(optimizers::Arg &) const;

   virtual double derivByParam(optimizers::Arg & dir,
                               const std::string & paramName) const {
// There is only the normalization, so the derivative is easy:
      return value(dir)/getParamValue(paramName);
   }

   virtual MapCubeFunction * clone() const {
      return new MapCubeFunction(*this);
   }

   void readFitsFile(const std::string & fitsFile);

   const std::string & fitsFile() const {
      return m_fitsFile;
   }

   /// @todo Consider moving this method to FitsImage class.
   double mapIntegral() const;

private:

   std::string m_fitsFile;

   std::string m_coordSys;
   std::vector<double> m_image;
   std::vector<double> m_lon;
   double m_lonMin, m_lonMax;
   std::vector<double> m_lat;
   double m_latMin, m_latMax;
   std::vector<double> m_energies;

   void init();

   void readEnergyVector(const std::string & fitsFile);

   static int findIndex(std::vector<double> xx, double x);

   double powerLawIntegral(double x1, double x2, double y1, double y2) const;
};

} // namespace Likelihood

#endif // Likelihood_MapCubeFunction_h
