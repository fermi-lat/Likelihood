/** 
 * @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.6 2003/07/19 04:38:02 jchiang Exp $
 *
 */

#ifndef Likelihood_SpatialMap_h
#define Likelihood_SpatialMap_h

#include "optimizers/Function.h"
#include "Likelihood/SkyDirArg.h"
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.6 2003/07/19 04:38:02 jchiang Exp $
 *
 */
    
class SpatialMap : public optimizers::Function, public FitsImage {
public:

   SpatialMap(std::string fitsfile);
   SpatialMap(const SpatialMap &rhs) : optimizers::Function(rhs), FitsImage(rhs) {
      m_ra = rhs.m_ra;
      m_dec = rhs.m_dec;
   }

   virtual ~SpatialMap() {}

   double value(optimizers::Arg&) const;

   double derivByParam(optimizers::Arg &, const std::string &) const
      {return 0;}

   virtual optimizers::Function *clone() const {
      return new SpatialMap(*this);
   }

private:

   std::vector<double> m_ra;
   std::vector<double> m_dec;

   // disable this
   double integral(optimizers::Arg &, optimizers::Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // Likelihood_SpatialMap_h
