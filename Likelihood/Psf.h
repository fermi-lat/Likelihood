/** @file Psf.h
 * @brief Psf class implementation
 * @author J. Chiang
 * 
 * $Header$
 */

#ifndef Psf_h
#define Psf_h

#include "astro/SkyDir.h"
#include "Likelihood/Response.h"
#include "Likelihood/Table.h"

namespace Likelihood {

/** 
 * @class Psf
 *
 * @brief LAT point spread function class
 *
 * @author J. Chiang
 *    
 * $Header$
 */

class Psf : public Response {
    
public:

   virtual ~Psf(){}

   //! PSF in instrument coordinates
   double value(double separation, double energy, double inc);
   double operator()(double separation, double energy, double inc)
      {return value(separation, energy, inc);}

   //! PSF in sky coordinates
   double value(astro::SkyDir appDir, double energy, 
                astro::SkyDir srcDir, double time);
   double operator()(astro::SkyDir appDir, double energy, 
                     astro::SkyDir srcDir, double time)
      {return value(appDir, energy, srcDir, time);}

   //! retrieve PSF parameters (sig1, sig2, wt) in instrument coordinates
   void fillPsfParams(double energy, double inclination,
                      std::vector<double> &psf_params);

   //! returns the Singleton object pointer
   static Psf * instance();

   //! method to read in the psf data
   void readPsfData(const std::string &psfFile, int hdu);

protected:

   Psf(){}

private:

   //! pointer to the Singleton instantiation
   static Psf * s_instance;

   //! PSF stored in straw-man CALDB format
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
