/** 
 * @file Psf.h
 * @brief Psf class implementation
 * @author J. Chiang
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Psf.h,v 1.12 2003/06/10 18:18:29 burnett Exp $
 */

#ifndef Psf_h
#define Psf_h

#include "Likelihood/Response.h"
#include "Likelihood/Table.h"

namespace Likelihood {

class astro::SkyDir;
class LikelihoodException;

/** 
 * @class Psf
 *
 * @brief LAT point spread function class
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Psf.h,v 1.12 2003/06/10 18:18:29 burnett Exp $
 */

class Psf : public Response {
    
public:

   virtual ~Psf(){}

   //! PSF in instrument coordinates
   double value(double separation, double energy, double inc);
   double operator()(double separation, double energy, double inc)
      {return value(separation, energy, inc);}

   //! PSF in sky coordinates
   double value(const astro::SkyDir &appDir, double energy, 
                const astro::SkyDir &srcDir, double time);
   double operator()(const astro::SkyDir &appDir, double energy, 
                     const astro::SkyDir &srcDir, double time)
      {return value(appDir, energy, srcDir, time);}

   //! retrieve PSF parameters (sig1, sig2, wt) in instrument coordinates
   void fillPsfParams(double energy, double inclination,
                      std::vector<double> &psf_params);

   //! returns the Singleton object pointer
   static Psf * instance();

   //! method to read in the psf data
   void readPsfData(const std::string &psfFile, int hdu)
      throw(LikelihoodException);

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
   std::valarray<double> m_sig1, m_sig2;

};

} // namespace Likelihood

#endif // Psf_h
