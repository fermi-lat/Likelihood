/** 
 * @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.25 2012/04/17 20:28:11 jchiang Exp $
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


   double value(const astro::SkyDir &) const;

   virtual optimizers::Function * clone() const {
      return new SpatialMap(*this);
   }

   const std::string & fitsFile() const {
      return m_fitsFile;
   }

   virtual double mapIntegral(double energy) const {
      (void)(energy);
      return wcsmap().mapIntegral();
   }

   virtual void integrateSpatialDist(const std::vector<double> & energies,
                                     const ExposureMap & expmap,
                                     std::vector<double> & exposure) const;

protected:

   double value(optimizers::Arg &) const;

   double derivByParamImp(optimizers::Arg &, const std::string &) const {
      return 0;
   }

};

} // namespace Likelihood

#endif // Likelihood_SpatialMap_h
