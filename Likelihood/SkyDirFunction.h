/** @file SkyDirFunction.h
 * @brief Declaration of SkyDirFunction class
 * $Header:
 */

#ifndef SkyDirFunction_h
#define SkyDirFunction_h

#include "../Likelihood/Function.h"
#include "astro/SkyDir.h"

namespace Likelihood {
/** 
 * @class SkyDirFunction
 *
 * @brief A class that encapsulates sky location information in a
 * Function context.  This allows sky coordinates, e.g., (ra, dec) or
 * (l, b), to be treated by the Likelihood::SourceModel class as 
 * Function Parameters.
 *
 * @author J. Chiang
 *    
 * $Header: 
 */
    
class SkyDirFunction : public Function, public astro::SkyDir {
public:

   SkyDirFunction() {m_init(0., 0.);}
   SkyDirFunction(double lon, double lat) {m_init(lon, lat);}
   SkyDirFunction(const astro::SkyDir &dir);

   astro::SkyDir getDir() const {return m_dir;}

   virtual double value(double) const {return 0.;}
   virtual double operator()(double x) const {return value(x);}

   //! overloaded setParam methods to include updating of m_dir
   void setParam(const std::string &paramName, double paramValue, 
		 bool isFree) {
      setParameter(paramName, paramValue, isFree);
      update_m_dir(paramName, paramValue);
   }
   void setParam(const std::string &paramName, double paramValue) {
      setParameter(paramName, paramValue);
      update_m_dir(paramName, paramValue);
   }

//   virtual double derivByParam(double, const std::string &paramName) const;

private:

   void m_init(double lon, double lat);
   
   void update_m_dir(std::string paramName, double paramValue);

   astro::SkyDir::CoordSystem m_coord_type;
   double m_lon;
   double m_lat;

   astro::SkyDir m_dir;

};

} // namespace Likelihood

#endif // SkyDirFunction_h
