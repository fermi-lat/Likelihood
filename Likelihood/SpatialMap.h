/** 
 * @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.11 2004/09/28 04:32:25 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.11 2004/09/28 04:32:25 jchiang Exp $
 *
 */
    
class SpatialMap : public optimizers::Function {

public:

   SpatialMap() {init();}

   SpatialMap(const std::string &fitsFile) {
      init();
      readFitsFile(fitsFile);
   }

   SpatialMap(const SpatialMap &rhs) : optimizers::Function(rhs) {
      m_ra = rhs.m_ra;
      m_dec = rhs.m_dec;
      m_image.resize(rhs.m_image.size());
      m_image = rhs.m_image;
      m_fitsFile = rhs.m_fitsFile;
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
   std::vector<double> m_dec;
   std::vector<double> m_image;

   void init();

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // Likelihood_SpatialMap_h
