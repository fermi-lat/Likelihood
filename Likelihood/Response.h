#ifndef Response_h
#define Response_h

#include "astro/SkyDir.h"
#include "../Likelihood/ScNtuple.h"

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
    
   Response();
   Response(const std::string &scFile, int scHdu);
   virtual ~Response(){};

protected:

   //! share the spacecraft data among all response functions
   static std::vector<ScNtuple> m_scData;
   static bool m_haveScData;
        
   //! the NR hunt routine 
   static void m_hunt(double xx[], int nx, double x, int *i);

   //! and my own zeroth order bilinear interpolater
   static double m_bilinear(int nx, double *xx, int i, double x,
			    int ny, double *yy, int j, double y, double *z);

private:

   std::string m_scFile;
   int m_scHdu;

   //! accessing the spacecraft data
   void m_readScData(const std::string &file, int hdu);

};

} // namespace Likelihood
#endif // Response_h
