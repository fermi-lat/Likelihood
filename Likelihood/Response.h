/** @file Response.h
 * @brief Response base class declaration
 * @author J. Chiang
 * 
 * $Header$
 */

#ifndef Response_h
#define Response_h

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
 * $Header$
 */

class Response {
    
public:
    
   virtual ~Response(){}

   //! type-fields for specifying response file HDUs
   enum HDU {Front = 2, Back, Combined};

   //! return the maximum allowed value of the source inclination wrt
   //! the instrument z-axis (in degrees)
   static double incMax() {return s_incMax;}

protected:

   Response();

   //! maximum inclination for response files
   static double s_incMax;

   //! share the spacecraft data among all response functions
   ScData * scData;
   
   //! the NR hunt routine 
   static void m_hunt(double *xx, int nx, double x, int *i);

   //! and my own zeroth order bilinear interpolater
   static double m_bilinear(int nx, double *xx, int i, double x,
                            int ny, double *yy, int j, double y, double *z);

};

} // namespace Likelihood

#endif // Response_h
