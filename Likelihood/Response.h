#ifndef Response_h
#define Response_h

#include "astro/SkyDir.h"
#include "../Likelihood/ScData.h"

namespace Likelihood {

/** 
 * @class Response
 *
 * @brief Base class for the LAT instrument response functions.
 * Spacecraft info is shared among all derived classes via static data
 * members.
 *
 * Should each derived class -- Aeff, Psf, Edisp -- be Singleton?
 *
 * @author J. Chiang
 *    
 * $Header: */

class Response {
    
public:
    
   virtual ~Response(){};

   //! type-fields for specifying response file HDUs
   enum HDU {Front = 2, Back, Combined};

protected:

   Response();

   //! share the spacecraft data among all response functions
   ScData * scData;
   
   //! the NR hunt routine 
   static void m_hunt(double xx[], int nx, double x, int *i);

   //! and my own zeroth order bilinear interpolater
   static double m_bilinear(int nx, double *xx, int i, double x,
			    int ny, double *yy, int j, double y, double *z);

};

} // namespace Likelihood
#endif // Response_h
