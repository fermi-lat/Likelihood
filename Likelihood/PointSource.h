/** 
 * @file PointSource.h
 * @brief PointSource class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PointSource.h,v 1.20 2003/07/19 04:38:01 jchiang Exp $
 */

#ifndef Likelihood_PointSource_h
#define Likelihood_PointSource_h

#include "Likelihood/Source.h"
#include "optimizers/Function.h"
#include "Likelihood/SkyDirFunction.h"
#include "Likelihood/Event.h"
#include "optimizers/dArg.h"

namespace Likelihood {

/** 
 * @class PointSource
 *
 * @brief A gamma-ray point source, e.g., pulsars, AGNs, X-ray binaries.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/PointSource.h,v 1.20 2003/07/19 04:38:01 jchiang Exp $
 */

class PointSource : public Source {

   friend class ExposureMap;

public:

   //! The default constructor does not compute exposure since 
   //! it assumes that the spacecraft data is not available or
   //! that the ROI cuts have not been specified.
   // A later setDir(ra, dec, true) will force the exposure to
   // be computed.
   PointSource() : m_spectrum(0) {setDir(0., 0., false); m_srcType = "Point";}

   //! This constructor does ask for exposure to be computed and 
   //! therefore *requires* the spacecraft data to be available and the 
   //! ROI cuts to be specified beforehand.
   PointSource(double ra, double dec) : m_spectrum(0) 
      {setDir(ra, dec);  m_srcType = "Point";}

   PointSource(const PointSource &rhs);

   virtual ~PointSource() {
      delete m_spectrum;
   }

   //! Returns photons/cm^2-s-sr-MeV having been convolved through
   //! the LAT instrument response
   double fluxDensity(const Event &evt) const
      {return fluxDensity(evt.getEnergy(), evt.getArrTime(), evt.getDir());}

   double fluxDensity(double energy, double time,
                      const astro::SkyDir &dir) const;

   //! Returns the derivative wrt to the named Parameter
   double fluxDensityDeriv(const Event &evt, std::string &paramName) const
      {return fluxDensityDeriv(evt.getEnergy(), evt.getArrTime(), 
                               evt.getDir(), paramName);}

   double fluxDensityDeriv(double energy, double time,
                           const astro::SkyDir &dir,
                           std::string &paramName) const;

   //! Predicted number of photons given RoiCuts and ScData
   virtual double Npred();

   //! Derivative of Npred wrt named Parameter
   virtual double NpredDeriv(const std::string &paramName);

   //! Set source location using J2000 coordinates
   void setDir(double ra, double dec, bool updateExposure = true) {
      m_dir = SkyDirFunction(astro::SkyDir(ra, dec));
      m_functions["Position"] = &m_dir;
      if (updateExposure) computeExposure();
   }

   //! Set source location via SkyDir class
   void setDir(const astro::SkyDir &dir, bool updateExposure = true) {
      m_dir = SkyDirFunction(dir);
      m_functions["Position"] = &m_dir;
      if (updateExposure) computeExposure();
   }

   astro::SkyDir getDir() const {return m_dir.getDir();}

   //! Angular separation between the source direction and dir in radians
   double getSeparation(const astro::SkyDir &dir) 
      {return dir.SkyDir::difference(m_dir.getDir());}

   //! Set the spectral model (should also check that the Parameter
   //! names do not conflict with "longitude" and "latitude" of m_dir)
   void setSpectrum(optimizers::Function *spectrum) {
      m_spectrum = spectrum->clone();
      m_functions["Spectrum"] = m_spectrum;
   }

   virtual Source *clone() const {
      return new PointSource(*this);
   }

protected:
   //! Computes the integrated exposure at the PointSource sky location.
   void computeExposure(int verbose = 1);

   //! Compute the integrated exposure using the provided 
   //! vector of energy values
   void computeExposure(std::vector<double> &energies,
                        std::vector<double> &exposure,
                        int verbose = 1);

//    //! provide access to the exposure values
//    void getExposure(std::vector<double> &exposure) const
//       {exposure = m_exposure;}

   //! location on the Celestial sphere 
   SkyDirFunction m_dir;

   //! spectral model
   optimizers::Function *m_spectrum;

private:

   //! flag to indicate that static member data has been computed
   static bool s_haveStaticMembers;

   //! vector of energy values for Npred spectrum quadrature
   static std::vector<double> s_energies;

   //! integrated exposure at PointSource sky location
   std::vector<double> m_exposure;

   //! method to create a logrithmically spaced grid given RoiCuts
   static void makeEnergyVector(int nee = 100);

   //! Gaussian widths in units of radians
   static std::vector<double> s_sigGauss;

   //! storage of the ROI contained fraction of a 2D "Gaussian"
   std::vector<double> m_gaussFraction;

   //! method to create the sigma grid for m_gaussFraction
   static void makeSigmaVector(int nsig = 100);

   //! Computes the enclosed fraction of a Gaussian as a function of 
   //! Gaussian width sigma and stores these data for later
   //! interpolation. These data depend on the ROI radius and the
   //! source offset, psi, from the ROI center.
   void computeGaussFractions();

   //! fraction of the psf that is contained within the ROI
   double psfFrac(double energy, double inc);

   //! nested class that returns the integrand for the m_gaussFraction
   //! integrals
   class Gint : public optimizers::Function {
   public:
      Gint() {};
      Gint(double sig, double cr, double cp, double sp) : 
           m_sig(sig), m_cr(cr), m_cp(cp), m_sp(sp) {}
      virtual ~Gint(){};
      double value(optimizers::Arg &mu) const;
      double derivByParam(optimizers::Arg &, const std::string &) const {return 0;}
   private:
      double m_sig;
      double m_cr;
      double m_cp;
      double m_sp;
   };

   //! a static object needed to compute the m_gaussFraction integrals
   //! using the DGAUS8 integrator
   static Gint s_gfunc;

   //! a static member function to provide the interface to s_gfunc
   //! that is required by DGAUS8
   static double gfuncIntegrand(double *mu) {
      optimizers::dArg muarg(*mu);
      return s_gfunc(muarg);
   }
};

} //namespace Likelihood

#endif // Likelihood_PointSource_h
