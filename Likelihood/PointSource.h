#ifndef PointSource_h
#define PointSource_h

#include "Source.h"
#include "Function.h"

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

   PointSource() : m_spectrum(0){setDir(0., 0.);};
   PointSource(const PointSource &rhs);
   virtual ~PointSource(){};

   //! returns photons/cm^2-s-sr-GeV having been convolved through
   //! the LAT instrument response
   virtual double fluxDensity(double energy, double time,
			      const astro::SkyDir &dir) const;

   //! set source location using J2000 coordinates
   void setDir(double ra, double dec);

   //! set source location via SkyDir class
   void setDir(const astro::SkyDir &dir) {m_dir = dir;};
   astro::SkyDir getDir() const {return m_dir;};
   double getSeparation(const astro::SkyDir &dir) const 
      {return m_dir.SkyDir::difference(dir);};

   //! set the spectral model
   void setSpectrum(Function *spectrum) {m_spectrum = spectrum;};
private:

   //! location on the Celestial sphere 
   astro::SkyDir m_dir;

   //! spectral model
   Function *m_spectrum;
        
};

} //namespace Likelihood

#endif // PointSource_h
