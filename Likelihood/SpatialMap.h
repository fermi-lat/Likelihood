/** @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.3 2003/05/02 19:02:14 jchiang Exp $
 *
 */

#ifndef SpatialMap_h
#define SpatialMap_h

#include "Likelihood/Function.h"
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
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SpatialMap.h,v 1.3 2003/05/02 19:02:14 jchiang Exp $
 *
 */
    
class SpatialMap : public Function, public FitsImage {
public:

   SpatialMap(std::string fitsfile);
   SpatialMap(const SpatialMap &rhs) : Function(rhs), FitsImage(rhs) {
      m_ra = rhs.m_ra;
      m_dec = rhs.m_dec;
   }

   virtual ~SpatialMap() {}

   double value(Arg&) const;

   double derivByParam(Arg &, const std::string &) const
      {return 0;}

   virtual Function *clone() const {
      return new SpatialMap(*this);
   }

private:

   std::vector<double> m_ra;
   std::vector<double> m_dec;

   // disable this
   double integral(Arg &, Arg &) const {return 0;}

};

} // namespace Likelihood

#endif // SpatialMap_h
