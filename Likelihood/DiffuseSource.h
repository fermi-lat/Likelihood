/** 
 * @file DiffuseSource.h
 * @brief DiffuseSource class declaration
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/DiffuseSource.h,v 1.55 2016/09/21 22:44:47 echarles Exp $
 */

#ifndef Likelihood_DiffuseSource_h
#define Likelihood_DiffuseSource_h

#include <cmath>

#include <stdexcept>

#include "optimizers/dArg.h"

#include "Likelihood/Source.h"
#include "Likelihood/SkyDirArg.h"
#include "Likelihood/TrapQuad.h"

namespace optimizers {
   class Function;
}

namespace Likelihood {

   class Event;
   class MapBase;
   class Observation;
   class WcsMap2;
   class HealpixProjMap;

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
		   bool requireExposure=true,
		   bool mapBasedIntegral=false);
     
     /// Copy c'tor
     DiffuseSource(const DiffuseSource &rhs);
     
     /// D'tor
     virtual ~DiffuseSource() {
       delete m_spatialDist;
     }

     /// Clone operator
     virtual Source *clone() const {
       return new DiffuseSource(*this);
     }

     /* ------------- Simple access functions ---------------------- */

     /// The spatial distribution
     inline const optimizers::Function * spatialDist() const {
       return m_spatialDist;
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

     /// Return the Map object for this diffuse source.
     /// This will throw an exception if the spatial distribution
     /// can not be cast to a MapBase
     const MapBase * mapBaseObject() const;

     /// Return the Map object for this diffuse source.
     /// This will throw an exception if the spatial distribution
     /// can not be cast to a MapBase
     MapBase * mapBaseObject();

     /// Return the flag that flux must be intergrated over a map
     inline bool mapBasedIntegral() const { m_mapBasedIntegral; }

     /// Return the size of the source, computing it if needed
     double mapRadius() const {
       if ( m_mapRadius < 0. ) {
	 m_mapRadius = computeMapRadius();
       }
       return m_mapRadius;
     }


     /* ------------------ Functions for unbinned likelihood --------------- */

     /// @return photons/cm^2-s-sr-MeV having been convolved through
     /// the LAT instrument response for a particular event
     /// @param evt container for event data
     /// @param cResp Cached instrument response
     virtual double fluxDensity(const Event & evt,
				CachedResponse * cResp = 0) const;
     
     /// @return derivative of fluxDensity wrt a model Parameter
     /// @param evt container for event data
     /// @param paramName name of the parameter in question
     /// @param cResp Cached instrument response   
     virtual double fluxDensityDeriv(const Event & evt, 
				     const std::string & paramName,
				     CachedResponse * cResp = 0) const;

     /// @return derivative of fluxDensity wrt a model Parameter
     /// This version is disabled
     virtual double fluxDensity(double inclination, double phi, double energy,
				const astro::SkyDir & appDir, int evtType,
				double time, CachedResponse* cResp = 0) const;
     
     virtual double fluxDensityDeriv(double inclination, double phi, double energy,
				     const astro::SkyDir & appDir, int evtType,
				     double time,
				     const std::string & paramName,
				     CachedResponse* cResp = 0) const;

     /// FIXME, what exactly does this do?
     virtual void computeExposure(const std::vector<double> & energies,
                                bool verbose=false);

     /* ------------ Functions for both binned and unbinned likelihood ------------ */


     /// @return Photon flux integrated over the ROI energy bounds. 
     /// Units are #/cm^2/s
     virtual double flux() const;
     
     /// @return Derivative of integrated photon flux wrt the named parameter
     virtual double fluxDeriv(const std::string & parName) const;
     
     /// @return Photon flux integrated over the given energy range.
     /// Units are #/cm^2/s
     virtual double flux(double emin, double emax, size_t npts=100) const;
     
     /// @return Derivative of integrated photon flux wrt the named parameter
     /// over the given energy range.
     virtual double fluxDeriv(const std::string & parName, 
			      double emin, double emax, size_t npts=100) const;
     
     /// @return Energy flux integrated over the ROI energy bounds. 
     /// Units are MeV/cm^2/s
     virtual double energyFlux() const;

     /// @return Derivative of integrated energy flux wrt the named parameter
     virtual double energyFluxDeriv(const std::string & parName) const;
     
     /// @return Energy flux integrated over the given energy range.
     /// Units are MeV/cm^2/s
     virtual double energyFlux(double emin, double emax, size_t npts=100) const;
     
     /// @return Derivative of integrated energy flux wrt the named parameter
     /// over the given energy range.
     virtual double energyFluxDeriv(const std::string & parName, 
				    double emin, double emax, size_t npts=100) const;

     /// @return the normalization factor for the spatial distribution
     virtual double angularIntegral(double energy) const;
     
     /// @return the diffuse response for a particular event
     double diffuseResponse(const Event & evt) const;

   protected:
     
     
     /* Version of diffuse response calculation for WCS maps */
     double diffuseResponse_wcs(const Event & evt, const WcsMap2& wcsmap) const;
     
     /* Version of diffuse response calculation for HEALPix maps */
     double diffuseResponse_healpix(const Event & evt, const HealpixProjMap& healmap) const;
     
     /* Compute the size of this source */
     double computeMapRadius() const;
     
     
   private:
     
     /// spatial model
     optimizers::Function * m_spatialDist;
     
     /// Flag to show that flux must be intergrated over a map
     bool m_mapBasedIntegral;
     
     /// Cache the value of the map radius, as it is kinda a pain to compute
     mutable double m_mapRadius;     
     
     // This function fills the Source::m_exposure array which contains
     // the integral over the source region of the spatial distribution
     // of the source times the unbinned exposure map.
     void integrateSpatialDist();
     
     /// @return Approximate energy-dependent outer angle (in radians) 
     ///         for diffuse response calculation.
     double psfRange(double energy) const;
   };
  
} //namespace Likelihood

#endif // Likelihood_DiffuseSource_h
