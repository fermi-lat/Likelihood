/** 
 * @file DiffuseSource.h
 * @brief DiffuseSource class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/DiffuseSource.h,v 1.34 2006/12/20 02:00:28 jchiang Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/DiffuseSource.h,v 1.34 2006/12/20 02:00:28 jchiang Exp $ 
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
   DiffuseSource(optimizers::Function * spatialDist,
                 const Observation & observation,
                 bool requireExposure = true);

   DiffuseSource(const DiffuseSource &rhs);

   virtual ~DiffuseSource() {
      delete m_spatialDist;
   }

   /// Returns photons/cm^2-s-sr-MeV having been convolved through
   /// the LAT instrument response
   virtual double fluxDensity(const Event &evt,
			      CachedResponse* cResp = 0) const;

   /// Returns the derivative wrt to the named Parameter
   virtual double fluxDensityDeriv(const Event &evt, 
                                   const std::string &paramName,
				   CachedResponse* cResp = 0) const;

   virtual double fluxDensity(double inclination, double phi, double energy,
                              const astro::SkyDir &appDir, int evtType,
			      CachedResponse* cResp = 0) const {
      (void)(inclination);
      (void)(phi);
      (void)(energy);
      (void)(appDir);
      (void)(evtType);
      (void)(cResp);
      return 0;
   }

   virtual double fluxDensityDeriv(double inclination, double phi, 
                                   double energy,
                                   const astro::SkyDir &appDir,
                                   int evtType,
                                   const std::string & paramName,
				   CachedResponse* cResp = 0) const {
      (void)(inclination);
      (void)(phi);
      (void)(energy);
      (void)(appDir);
      (void)(evtType);
      (void)(paramName);
      (void)(cResp);
      return 0;
   }

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

   virtual Source *clone() const {
      return new DiffuseSource(*this);
   }

   virtual double pixelCounts(double emin, double emax,
                              double wtMin, double wtMax) const;

   virtual double pixelCountsDeriv(double emin, double emax,
                                   double wtMin, double wtMax,
                                   const std::string & paramName) const;

   virtual const std::vector<double> & exposure() const {
      return m_exposure;
   }

   virtual const optimizers::Function & spectrum() const {
      return * m_spectrum;
   }

private:

   /// spatial model
   optimizers::Function * m_spatialDist;

};

} //namespace Likelihood

#endif // Likelihood_DiffuseSource_h
