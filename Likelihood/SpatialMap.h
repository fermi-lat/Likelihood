/** 
 * @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.15 2005/10/03 15:02:37 jchiang Exp $
 *
 */

#ifndef Likelihood_SpatialMap_h
#define Likelihood_SpatialMap_h

#include "optimizers/Function.h"

namespace astro {
   class SkyProj;
}

namespace Likelihood {

/** 
 * @class SpatialMap
 *
 * @brief A Function object that returns the value of a FITS image
 * at a SkyDir location.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.15 2005/10/03 15:02:37 jchiang Exp $
 *
 */
    
class SpatialMap : public optimizers::Function {

public:

   SpatialMap() : m_proj(0) {
      init();
   }

   SpatialMap(const std::string & fitsFile, const std::string & extension="")
      : m_proj(0) {
      init();
      readFitsFile(fitsFile, extension);
   }

   SpatialMap::SpatialMap(const SpatialMap &);

   SpatialMap & operator=(const SpatialMap &);

   virtual ~SpatialMap();

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &, const std::string &) const {
      return 0;
   }

   void readFitsFile(const std::string & fitsFile,
                     const std::string & extension="");

   virtual optimizers::Function *clone() const {
      return new SpatialMap(*this);
   }

   const std::string & fitsFile() const {
      return m_fitsFile;
   }

private:

   std::string m_fitsFile;
   std::string m_extension;

   std::vector<float> m_image;

   int m_naxis1; 
   int m_naxis2; 
   int m_naxis3;

   astro::SkyProj * m_proj;

   void init();

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // Likelihood_SpatialMap_h
