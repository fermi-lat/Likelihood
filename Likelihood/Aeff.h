/** @file Aeff.h
 * @brief Interface definition for the LAT Effective Area class
 * @author J. Chiang
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Aeff.h,v 1.7 2003/03/17 00:53:42 jchiang Exp $
 *
 */

#ifndef Aeff_h
#define Aeff_h

#include "astro/SkyDir.h"
#include "Likelihood/Response.h"
#include "Likelihood/Table.h"

namespace Likelihood {

/** 
 * @class Aeff
 *
 * @brief LAT effective area class.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Aeff.h,v 1.7 2003/03/17 00:53:42 jchiang Exp $
 */

class Aeff : public Response {
    
public:

   virtual ~Aeff(){};

   //! return the effective area in instrument coordinates
   double value(double energy, double inclination);
   double operator()(double energy, double inclination)
      {return value(energy, inclination);};

   //! effective area in sky coordinates
   double value(double energy, astro::SkyDir srcDir, double time);
   double operator()(double energy, astro::SkyDir srcDir, double time)
      {return value(energy, srcDir, time);};

   //! returns the Singleton object pointer
   static Aeff * instance();

   //! method to read in the aeff data
   void readAeffData(const std::string &aeffFile, int hdu);

protected:

   Aeff(){};

private:

   //! pointer to the Singleton instantiation
   static Aeff * s_instance;

   //! effective area stored in straw-man CALDB format
   std::string m_aeffFile;
   int m_aeffHdu;
   Table m_aeffData;

   std::vector<double> m_energy;
   std::vector<double> m_theta;
   std::valarray<double> m_aeff;

};

} // namespace Likelihood

#endif // Aeff_h
