#ifndef Aeff_h
#define Aeff_h

#include "astro/SkyDir.h"
#include "Response.h"
#include "Table.h"

namespace Likelihood {

/** 
 * @class Aeff
 *
 * @brief LAT effective area class.
 *
 * @author J. Chiang
 *    
 * $Header:
 */

class Aeff : public Response {
    
public:

   Aeff();
   Aeff(const std::string &scFile, int hdu) : Response() 
      {m_readAeffData(scFile, hdu);};
   virtual ~Aeff(){};

   //! return the effective area in instrument coordinates
   double value(double energy, double inclination);
   double operator()(double energy, double inclination)
      {return value(energy, inclination);};

   //! effective area in sky coordinates
   double value(double energy, astro::SkyDir &srcDir, double time);
   double operator()(double energy, astro::SkyDir &srcDir, double time)
      {return value(energy, srcDir, time);};

private:

   //! effective area stored in straw-man CALDB format
   void m_readAeffData(const std::string &file, int hdu);
   std::string m_aeffFile;
   int m_aeffHdu;
   Table m_aeffData;

   std::vector<double> m_energy;
   std::vector<double> m_theta;

   //! need to find a better way to store the effective area table,
   //! perhaps with std::valarray....for now, use a pointer
   double *m_aeff;

};

} // namespace Likelihood
#endif // Aeff_h
