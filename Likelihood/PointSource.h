#ifndef PointSource_h
#define PointSource_h

#include "../Likelihood/Source.h"
#include "../Likelihood/Function.h"
#include "../Likelihood/SkyDirFunction.h"
#include "../Likelihood/Event.h"
#include "../Likelihood/Arg.h"

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

   PointSource() : m_spectrum(0) {
      setDir(0., 0.);
      m_makeEnergyVector();
   }
   PointSource(double ra, double dec) : m_spectrum(0) {
      setDir(ra, dec);
      m_makeEnergyVector();
   }
   PointSource(const PointSource &rhs);
   virtual ~PointSource(){}

   //! returns photons/cm^2-s-sr-MeV having been convolved through
   //! the LAT instrument response
   virtual double fluxDensity(double energy, double time,
                              const astro::SkyDir &dir) const;

   double fluxDensity(const Event &evt) const
      {return fluxDensity(evt.getEnergy(), evt.getArrTime(), evt.getDir());}

   //! predicted number of photons given RoiCuts and ScData
   double Npred();

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
   //! pass a pointer?...should also check that the Parameter names do
   //! not conflict with "longitude" and "latitude" of m_dir)
   void setSpectrum(Function *spectrum) {
      m_spectrum = spectrum;
      m_functions["Spectrum"] = m_spectrum;
   }

protected:

   //! location on the Celestial sphere 
   SkyDirFunction m_dir;

   //! spectral model
   Function *m_spectrum;

private:

   //! vector of energy values for Npred spectrum quadrature
   std::vector<double> m_energies;

   //! method to create a logrithmically spaced grid given RoiCuts
   void m_makeEnergyVector(int nee = 100);
        
};

} //namespace Likelihood

#endif // PointSource_h
