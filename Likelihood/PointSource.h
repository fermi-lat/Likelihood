#ifndef PointSource_h
#define PointSource_h

#include "Source.h"
#include "Function.h"
#include "SkyDirFunction.h"

namespace Likelihood {

/** 
 * @class PointSource
 *
 * @brief A gamma-ray point source, e.g., pulsars, AGNs, X-ray binaries.
 *
 * @author J. Chiang
 *    
 * $Header:
 */

class PointSource : public Source {

public:

   PointSource() : m_spectrum(0){setDir(0., 0.);}
   PointSource(double ra, double dec) : m_spectrum(0){setDir(ra, dec);}
   PointSource(const PointSource &rhs);
   virtual ~PointSource(){}

   //! returns photons/cm^2-s-sr-MeV having been convolved through
   //! the LAT instrument response
   double fluxDensity(double energy, double time,
                      const astro::SkyDir &dir) const;

   //! the predicted number of photons for a given ROI
   //! the vector<double> is a place holder for an RoiCuts class object
   double Npred(const std::vector<double>) const {return 0.;}

   //! set source location using J2000 coordinates
   void setDir(double ra, double dec) {
      m_dir = SkyDirFunction(astro::SkyDir(ra, dec));
      m_functions["Position"] = &m_dir;
   }

   //! set source location via SkyDir class
   void setDir(const astro::SkyDir &dir) {
      m_dir = SkyDirFunction(dir);
      m_functions["Position"] = &m_dir;
   }

   astro::SkyDir getDir() const {return m_dir.getDir();}

   //! angular separation between the source direction and dir in radians
   double getSeparation(const astro::SkyDir &dir) 
      {return dir.SkyDir::difference(m_dir.getDir());};

   //! set the spectral model (do we want to clone spectrum rather than 
   //! pass a pointer?)
   void setSpectrum(Function *spectrum) {
      m_spectrum = spectrum;
      m_functions["Spectrum"] = m_spectrum;
   }

protected:

   //! location on the Celestial sphere 
   SkyDirFunction m_dir;

   //! spectral model
   Function *m_spectrum;
        
};

} //namespace Likelihood

#endif // PointSource_h
