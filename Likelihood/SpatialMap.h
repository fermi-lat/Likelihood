/** @file SpatialMap.h
 * @brief Declaration for the SpatialMap Function class
 * @author J. Chiang
 *
 * $Header$
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
 * $Header$
 *
 */
    
class SpatialMap : public Function, public FitsImage {
public:

   SpatialMap(std::string fitsfile);
   virtual ~SpatialMap() {}

   double value(Arg&) const;

   double derivByParam(Arg &x, const std::string &paramName) const
      {return 0;}

   virtual Function *clone() const {
      return new SpatialMap(*this);
   }

private:

   std::vector<double> m_ra;
   std::vector<double> m_dec;

   // disable this
   double integral(Arg &xmin, Arg &xmax) const {return 0;}

};

} // namespace Likelihood

#endif // SpatialMap_h
