/** 
 * @file DiffuseSource.h
 * @brief DiffuseSource class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DiffuseSource.h,v 1.26 2005/02/27 06:42:24 jchiang Exp $
 */

#ifndef Likelihood_DiffuseSource_h
#define Likelihood_DiffuseSource_h

#include <stdexcept>

#include "Likelihood/Source.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/Exception.h"

namespace optimizers {
   class Function;
}

namespace Likelihood {

class Event;
class Observation;

/** 
 * @class DiffuseSource
 *
 * @brief Non-point-like gamma-ray source.  The extragalactic and
 * Galactic diffuse emission components are the obvious examples, but
 * this class also includes discrete diffuse sources such as the LMC
 * or supernova remnants.
 *
 * This representation assumes that a single spectral model describes
 * the emission over the entire angular extent of the source.
 * Therefore, in order to have spectral variations across a source, it
 * must comprise sub-components that can be represented using this
 * class.
 *
 * Note that the member function double spatialDist(astro::SkyDir &)
 * returns the spatial distribution of the emission as a function of
 * direction.
 *
 * @author J. Chiang
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DiffuseSource.h,v 1.26 2005/02/27 06:42:24 jchiang Exp $ 
 *  
 */

class DiffuseSource : public Source {

public:

   /// @param spatialDist A Function object describing the spatial
   ///        distribution of the emission.
   /// @param requireExposure If true, then an ExposureMap must have
   ///        been read in so that the spatially integrated for the
   ///        spectral calculation; if false, then the map is not
   ///        integrated.
   DiffuseSource(optimizers::Function *spatialDist,
                 const Observation & observation,
                 bool requireExposure = true);

   DiffuseSource(const DiffuseSource &rhs);

   virtual ~DiffuseSource() {
      delete m_spatialDist;
      delete m_spectrum;
   }

   /// Returns photons/cm^2-s-sr-MeV having been convolved through
   /// the LAT instrument response
   virtual double fluxDensity(const Event &evt) const;

   /// Returns the derivative wrt to the named Parameter
   virtual double fluxDensityDeriv(const Event &evt, 
                                   const std::string &paramName) const;

   virtual double fluxDensity(double inclination, double phi, double energy,
                              const astro::SkyDir &appDir, int evtType) const {
      (void)(inclination);
      (void)(phi);
      (void)(energy);
      (void)(appDir);
      (void)(evtType);
      return 0;
   }


   virtual double fluxDensityDeriv(double inclination, double phi, 
                                   double energy,
                                   const astro::SkyDir &appDir,
                                   int evtType,
                                   const std::string & paramName) const {
      (void)(inclination);
      (void)(phi);
      (void)(energy);
      (void)(appDir);
      (void)(evtType);
      (void)(paramName);
      return 0;
   }

   /// Predicted number of photons given RoiCuts and ScData
   virtual double Npred();
   
   /// Derivative of Npred wrt named Parameter
   virtual double NpredDeriv(const std::string &paramName);

   /// Predicted number of counts within a specified energy range
   virtual double Npred(double emin, double emax);

   /// Return the spatial distribution of the gamma-ray emission
   double spatialDist(const astro::SkyDir & dir) const {
      SkyDirArg SDarg(dir);
      return (*m_spatialDist)(SDarg);
   }

#ifndef SWIG
   double spatialDist(SkyDirArg dir) const {
      return (*m_spatialDist)(dir);
   }
#endif

   /// Set the spectral model (should also check that the Parameter
   /// names do not conflict with "longitude" and "latitude" of m_dir)
   void setSpectrum(optimizers::Function *spectrum) {
      m_spectrum = spectrum->clone();
      m_functions["Spectrum"] = m_spectrum;
   }

   virtual Source *clone() const {
      return new DiffuseSource(*this);
   }

   virtual double pixelCounts(double emin, double emax,
                              double wtMin, double wtMax) const;

   virtual double pixelCountsDeriv(double emin, double emax,
                                   double wtMin, double wtMax,
                                   const std::string & paramName) const;

protected:

   /// spatial model
   optimizers::Function *m_spatialDist;

   /// spectral model
   optimizers::Function *m_spectrum;

private:

   /// flag to indicate that static member data have been computed
   static bool s_haveStaticMembers;

   /// vector of energy values for Npred spectrum quadrature
   static std::vector<double> s_energies;

   /// method to create a logrithmically spaced grid of energies
   static void makeEnergyVector(double emin, double emax, int nee = 100);

   /// vector of angle integrated diffuse exposure as as a function of
   /// energy
   std::vector<double> m_exposure;
};

} //namespace Likelihood

#endif // Likelihood_DiffuseSource_h
