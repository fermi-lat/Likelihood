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
      if (!m_haveStaticMembers) {
         m_makeEnergyVector();
         m_makeSigmaVector();
         m_haveStaticMembers = true;
	 computeGaussFractions();
         computeExposure();
      }
   }
   PointSource(double ra, double dec) : m_spectrum(0) {
      setDir(ra, dec);
      if (!m_haveStaticMembers) {
         m_makeEnergyVector();
         m_makeSigmaVector();
         m_haveStaticMembers = true;
	 computeGaussFractions();
         computeExposure();
      }
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

   //! derivative of Npred wrt named Parameter
   double NpredDeriv(const std::string &paramName);

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

   //! computes the enclosed fraction of a Gaussian as a function of 
   //! Gaussian width sigma and stores these data for later
   //! interpolation. These data depend on the ROI radius and the
   //! source offset, psi, from the ROI center.
   void computeGaussFractions();

   //! computes the integrated exposure at the PointSource sky location
   void computeExposure();

protected:

   //! location on the Celestial sphere 
   SkyDirFunction m_dir;

   //! spectral model
   Function *m_spectrum;

private:

   //! flag to indicate that static member data has been computed
   static bool m_haveStaticMembers;

   //! vector of energy values for Npred spectrum quadrature
   static std::vector<double> m_energies;

   //! integrated exposure at PointSource sky location
   std::vector<double> m_exposure;

   //! method to create a logrithmically spaced grid given RoiCuts
   static void m_makeEnergyVector(int nee = 100);

   //! Gaussian widths in units of radians
   static std::vector<double> m_sigGauss;

   //! method to create the sigma grid for m_gaussFraction
   static void m_makeSigmaVector(int nsig = 100);

   //! storage of the ROI contained fraction of a 2D "Gaussian"
   std::vector<double> m_gaussFraction;

   //! nested class that returns the desired integrand
   class Gint : public Function {
   public:
      Gint(double sig, double cr, double cp, double sp) : 
	 m_sig(sig), m_cr(cr), m_cp(cp), m_sp(sp) {}
      virtual ~Gint(){};
      double value(Arg &mu) const;
      double derivByParam(Arg &, const std::string &) const {return 0;}
   private:
      double m_sig;
      double m_cr;
      double m_cp;
      double m_sp;
   };
};

} //namespace Likelihood

#endif // PointSource_h
