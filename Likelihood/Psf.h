#ifndef Psf_h
#define Psf_h

#include "astro/SkyDir.h"
#include "Response.h"
#include "Table.h"

namespace Likelihood {

/** 
 * @class Psf
 *
 * @brief LAT point spread function class
 *
 * @author J. Chiang
 *    
 * $Header:
 */

class Psf : public Response {
    
public:

   Psf();
   Psf(const std::string &scFile, int hdu) : Response() 
      {m_readPsfData(scFile, hdu);};
   virtual ~Psf(){};

   //! PSF in instrument coordinates
   double value(double separation, double energy, double inc);
   double operator()(double separation, double energy, double inc)
      {return value(separation, energy, inc);};

   //! PSF in sky coordinates
   double value(astro::SkyDir &appDir, double energy, 
		astro::SkyDir &srcDir, double time);
   double operator()(astro::SkyDir &appDir, double energy, 
		     astro::SkyDir &srcDir, double time)
      {return value(appDir, energy, srcDir, time);};

   //! retrieve PSF parameters (sig1, sig2, wt) in instrument coordinates
   void fillPsfParams(double energy, double inclination,
		      std::vector<double> &psf_params);

private:

   //! effective area stored in straw-man CALDB format
   void m_readPsfData(const std::string &file, int hdu);
   std::string m_psfFile;
   int m_psfHdu;
   Table m_psfData;

   std::vector<double> m_energy;
   std::vector<double> m_theta;
   std::vector<double> m_wt;

   //! need to find a better way to store the psf width table data
   //! perhaps with std::valarray....for now, use pointers
   double *m_sig1, *m_sig2;

};

} // namespace Likelihood
#endif // Psf_h
