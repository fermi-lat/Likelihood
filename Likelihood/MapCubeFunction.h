/**
 * @file MapCubeFunction.h
 * @brief Encapsulation of a 3D FITS image.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/MapCubeFunction.h,v 1.14 2011/03/16 00:19:37 jchiang Exp $
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
 *
 * @brief Encapsulation of a 3D FITS image, e.g., with RA, Dec, energy
 * dimensions, for modeling diffuse sources with spectral variation as
 * a function of position on the sky.
 *
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
                             const std::string & extension="",
                             bool loadMap=true);

   virtual void readFitsFile();

   virtual void deleteMap();

   const std::string & fitsFile() const {
      return m_fitsFile;
   }

   /// @todo Consider moving this method to FitsImage class.
   double mapIntegral() const;

   /// @return The angular integral of the differential flux as a
   /// function of energy (implementation for
   /// MapBase::mapIntegral(double) const)
   virtual double mapIntegral(double energy) const;

   virtual void integrateSpatialDist(const std::vector<double> & energies,
                                     const ExposureMap & expmap,
                                     std::vector<double> & exposure) const {
      throw std::runtime_error("MapCubeFunction::integrateSpatialDist: "
                               "not implemented.");
   }

private:

   astro::SkyProj * m_proj;

   int m_nlon;
   int m_nlat;

   bool m_isPeriodic;

   std::vector<double> m_energies;

   std::vector<float> m_image;

   std::vector<double> m_mapIntegrals;

   void init();

   int findIndex(const std::vector<double> & xx, double x) const;

   double powerLawIntegral(double x1, double x2, double y1, double y2) const;

   void computeMapIntegrals();

   double bilinear_interpolation(size_t k, int ix, int iy,
                                 double x, double y) const;
};

} // namespace Likelihood

#endif // Likelihood_MapCubeFunction_h
