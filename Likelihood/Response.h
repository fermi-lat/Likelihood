/** @file Response.h
 * @brief Response base class declaration
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Response.h,v 1.13 2003/03/25 23:22:02 jchiang Exp $
 */

#ifndef Response_h
#define Response_h

#include <iostream>
#include <valarray>
#include "astro/SkyDir.h"
#include "Likelihood/ScData.h"

namespace Likelihood {

/** 
 * @class Response
 *
 * @brief Base class for the LAT instrument response functions.
 * Spacecraft info is shared among all derived classes via static data
 * members.
 *
 * Each derived class -- Aeff, Psf, Edisp -- shall be Singleton.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Response.h,v 1.13 2003/03/25 23:22:02 jchiang Exp $
 */

class Response {
    
public:
    
   virtual ~Response(){}

   //! type-fields for specifying response file HDUs
   enum HDU {Front = 2, Back, Combined};

   //! return the maximum allowed value of the source inclination wrt
   //! the instrument z-axis (in degrees)
   static double incMax() {return s_incMax;}

   //! my own zeroth order bilinear interpolater
   static double bilinear(const std::vector<double> &xx, double x,
                          const std::vector<double> &yy, double y, 
                          const std::valarray<double> &z);

protected:

   Response();

   //! maximum inclination for response files
   static double s_incMax;

   //! share the spacecraft data among all response functions
   ScData * scData;
   
};

} // namespace Likelihood

#endif // Response_h
