/** 
 * @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/Likelihood/SpatialMap.h,v 1.6 2015/11/25 18:52:41 echarles Exp $
 *
 */

#ifndef Likelihood_SpatialMap_h
#define Likelihood_SpatialMap_h

#include <utility>

#include "optimizers/Function.h"

#include "Likelihood/MapBase.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

class Event;
class ResponseFunctions;
// EAC, switch to projection-specific methods
class WcsMap2;
class HealpixProjMap;

/** 
 * @class SpatialMap
 *
 * @brief A Function object that returns the value of a FITS image
 * at a SkyDir location.
 *
 */
    
class SpatialMap : public optimizers::Function, public MapBase {

public:

   SpatialMap();

   SpatialMap(const std::string & fitsFile, const std::string & extension="");

   SpatialMap(const SpatialMap &);

   SpatialMap & operator=(const SpatialMap &);

   virtual ~SpatialMap();

   double value(const optimizers::Arg &) const;

   double value(const astro::SkyDir &) const;

   double derivByParamImp(const optimizers::Arg &, const std::string &) const {
      return 0;
   }

   virtual optimizers::Function * clone() const {
      return new SpatialMap(*this);
   }

   const std::string & fitsFile() const {
      return m_fitsFile;
   }

   virtual double mapIntegral(double energy) const {
      (void)(energy);
      return projmap().mapIntegral();
   }

   virtual void integrateSpatialDist(const std::vector<double> & energies,
                                     const ExposureMap & expmap,
                                     std::vector<double> & exposure) const;

   void integrateSpatialDist_wcs(const std::vector<double> & energies,
                                 const ExposureMap & expmap,
                                 const WcsMap2& wcsmap,
                                 std::vector<double> & exposure) const;
   
   void integrateSpatialDist_healpix(const std::vector<double> & energies,
				     const ExposureMap & expmap,
				     const HealpixProjMap& healmap,
				     std::vector<double> & exposure) const;


private:

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // Likelihood_SpatialMap_h
