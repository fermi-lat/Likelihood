/** 
 * @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.16 2005/10/04 05:38:08 jchiang Exp $
 *
 */

#ifndef Likelihood_SpatialMap_h
#define Likelihood_SpatialMap_h

#include "optimizers/Function.h"

namespace astro {
   class SkyDir;
}

namespace Likelihood {

class WcsMap;

/** 
 * @class SpatialMap
 *
 * @brief A Function object that returns the value of a FITS image
 * at a SkyDir location.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.16 2005/10/04 05:38:08 jchiang Exp $
 *
 */
    
class SpatialMap : public optimizers::Function {

public:

   SpatialMap();

   SpatialMap(const std::string & fitsFile, const std::string & extension="");

   SpatialMap(const SpatialMap &);

   SpatialMap & operator=(const SpatialMap &);

   virtual ~SpatialMap();

   double value(optimizers::Arg &) const;

   double value(const astro::SkyDir &) const;

   void readFitsFile(const std::string & fitsFile,
                     const std::string & extension="");

   double derivByParam(optimizers::Arg &, const std::string &) const {
      return 0;
   }

   virtual optimizers::Function * clone() const {
      return new SpatialMap(*this);
   }

   const std::string & fitsFile() const {
      return m_fitsFile;
   }

private:

   WcsMap * m_wcsmap;

   std::string m_fitsFile;
   std::string m_extension;

   void init();

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // Likelihood_SpatialMap_h
