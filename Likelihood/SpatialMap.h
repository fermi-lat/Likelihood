/** 
 * @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.13 2005/02/17 23:22:31 jchiang Exp $
 *
 */

#ifndef Likelihood_SpatialMap_h
#define Likelihood_SpatialMap_h

#include "optimizers/Function.h"

#include "Likelihood/FitsImage.h"

namespace Likelihood {
/** 
 * @class SpatialMap
 *
 * @brief A Function object that returns the value of a FITS image
 * at a SkyDir location.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.13 2005/02/17 23:22:31 jchiang Exp $
 *
 */
    
class SpatialMap : public optimizers::Function {

public:

   SpatialMap() {init();}

   SpatialMap(const std::string &fitsFile) {
      init();
      readFitsFile(fitsFile);
   }

   virtual ~SpatialMap() {}

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &, const std::string &) const
      {return 0;}

   void readFitsFile(const std::string &fitsFile);

   virtual optimizers::Function *clone() const {
      return new SpatialMap(*this);
   }

   const std::string & fitsFile() const {return m_fitsFile;}

private:

   std::string m_fitsFile;
   std::vector<double> m_ra;
   double m_raMin, m_raMax;
   std::vector<double> m_dec;
   double m_decMin, m_decMax;

   std::string m_coordSys;

   std::vector<double> m_image;

   void init();

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // Likelihood_SpatialMap_h
