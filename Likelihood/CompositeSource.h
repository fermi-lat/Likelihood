/** 
 * @file CompositeSource.h
 * @brief A source made of a collection of sources
 *
 * @author E. Charles
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Source.h,v 1.54 2016/09/09 21:11:52 echarles Exp $
 */

#ifndef Likelihood_CompositeSource_h
#define Likelihood_CompositeSource_h

#include "Likelihood/Source.h"

#include <string>
#include "Likelihood/SourceModel.h"
#include "Likelihood/BinnedConfig.h"


namespace astro {
  class SkyDir;
}

namespace Likelihood {

  class Observation;
  class SourceMapCache;
  class SourceMap;
  class BinnedCountsCache;
  class Drm;

  /** 
   * @class CompositeSource
   *
   * @brief Source made of other sources
   *
   */
  
  class CompositeSource : public Source {

  public:

    /* Default c'tor */
    CompositeSource(const Observation& observation=0);

    /* Standard c'tor

       observation  :   The container for the observation data
       name         :   Name of this composite
       srcs         :   Vector of the sources in this composite
       specFuncName :   Spectral function for this composite 
    */
    CompositeSource(const Observation& observation,
		    const std::string& name,
		    const std::string& specFuncName);

    
    /* Copy c'tor */
    CompositeSource(const  CompositeSource&rhs);
    
    /* D'tor */
    virtual ~CompositeSource();
  
    /* Clone operation */
    virtual Source * clone() const;
  
    /* -------------------------  Trival access --------------------------- */ 

    /* The configuration */
    inline const BinnedLikeConfig& config() const { return m_config; }

    /* The xml file to write this source to */
    inline const std::string& xmlFile() const { return m_xmlFile; }


    /* ---------------- Access to the SourceModel ------------------------- */

    inline const SourceModel& sourceModel() const { return m_sourceModel; }
    //inline SourceModel& sourceModel() { return m_sourceModel; }    
    //inline const std::map<std::string, Source *>& sources() const { return m_sourceModel.sources(); }

    /* Add a source to this composite.
       
       src       : The Source in question
       srcMap    : If not NULL, add this SourceMap to the SourceMapCache
       fromClone : If true, add a cloned version of the source
     */
    void addSource(Source *src, SourceMap* srcMap = 0, 
		   bool fromClone=false);
        
    /* Remove a source from this composite
       
       srcName   : Name of the source to remove
       srcMap    : Filled with the SourceMap for this source
       
    */
    Source * deleteSource(const std::string &srcName,
			  SourceMap*& srcMap);

    /* Remove everything and clears out the Source Map cache */
    void deleteAllSources();

    /* Get a particular source by name 
       
       This will return 0 if the source does not exist.
     */
    Source * getSource(const std::string &srcName) {
      return m_sourceModel.getSource(srcName);
    }
    
    /* Get a particular source by name 

       This will throw an exception if the source does not exist */
    const Source & source(const std::string & srcName) const {
      return m_sourceModel.source(srcName);
    }

    /* Fill a vector with pointer to sources, given their names 

       This calls getSource, so it will put NULL pointers on the vector
       for sources that are not in the composite
     */
    inline void getSources(const std::vector<std::string>& srcNames, 
			   std::vector<const Source*>& srcs) const {
      return  m_sourceModel.getSources(srcNames,srcs); 
    }

    /* return the number of sources in this composite */
    inline unsigned int getNumSrcs() const { return m_sourceModel.getNumSrcs(); }
    
    /* Fill a vector with the names of the source in this composite */
    inline void getSrcNames(std::vector<std::string> & srcNames) const {
      return m_sourceModel.getSrcNames(srcNames); 
    }

    /* Check if a particular source is part of this composite */
    bool hasSrcNamed(const std::string & srcName) const {
      return m_sourceModel.hasSrcNamed(srcName); 
    }

    
    /// Steal a source a SourceModel
    Source* steal_source(SourceModel& other,
			 const std::string& srcName);
    
    /// Give a source to a SourceModel
    Source* give_source(SourceModel& other,
			const std::string& srcName);
  
 

    /* ---------------- Interaction with the SourceMaps ------------------------- */

    /*  Build the source map cache */
    SourceMapCache* buildSourceMapCache(const BinnedCountsCache& dataCache,
					const std::string & srcMapsFile,
					const Drm* drm = 0);

 
    /* Access to the source map cache */
    inline const SourceMapCache* sourceMapCache() const { return m_srcMapCache; }


    /* ---------------- Methods inherited from Source class -------------*/
    
    /* ---------------- Methods for Unbinned Analysis ------------------ */

    /// @return photons/cm^2-s-sr-MeV having been convolved through
    /// the LAT instrument response for a particular event
    /// @param evt container for event data
    /// @param cResp Cached instrument response
    virtual double fluxDensity(const Event & evt, 
			       CachedResponse* cResp = 0) const;
    
    /// @return fluxDensity in instrument coordinates (photons/cm^2-s-sr-MeV)
    /// @param inclination angle of source direction wrt the instrument
    ///        z-axis (degrees)
    /// @param phi azimuthal angle of source direction wrt instrument
    ///        z- and x-axes (degrees)
    /// @param energy True energy of photon (MeV)
    /// @param appDir Apparent photon direction
    /// @param evtType Event type, i.e., front- vs back-converting event, 
    ///        0 vs 1
    /// @param time Event arrival time
    virtual double fluxDensity(double inclination, double phi, double energy, 
			       const astro::SkyDir & appDir, 
			       int evtType, double time, 
			       CachedResponse* cResp = 0) const;
    
    
    /// @return derivative of fluxDensity wrt a model Parameter
    /// @param evt container for event data
    /// @param paramName name of the parameter in question
    /// @param cResp Cached instrument response   
    virtual double fluxDensityDeriv(const Event & evt, 
				    const std::string & paramName,
				    CachedResponse* cResp = 0) const;
    
    /// @return derivative of fluxDensity wrt a model Parameter
    /// @param inclination angle of source direction wrt the instrument
    ///        z-axis (degrees)
    /// @param phi azimuthal angle of source direction wrt instrument
    ///        z- and x-axes (degrees)
    /// @param energy True energy of photon (MeV)
    /// @param appDir Apparent photon direction
    /// @param evtType Event type, i.e., front- vs back-converting event, 
    ///        0 vs 1
    /// @param time Event arrival time
    /// @param paramName name of the parameter in question
    /// @param cResp Cached instrument response
    virtual double fluxDensityDeriv(double inclination, double phi, 
				    double energy, const astro::SkyDir & appDir,
				    int evtType, double time, 
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
    virtual double fluxDeriv(const std::string & parName, double emin,
			     double emax, size_t npts=100) const;
    
    /// @return Energy flux integrated over the ROI energy bounds. 
    /// Units are MeV/cm^2/s
    virtual double energyFlux() const;
    
    /// @return Derivative of integrated energy flux wrt the named parameter
    virtual double energyFluxDeriv(const std::string & parName) const;
    
    /// @return Energy flux integrated over the given energy range.
    /// Units are MeV/cm^2/s
    virtual double energyFlux(double emin, double emax,
			      size_t npts=100) const;
    
    /// @return Derivative of integrated energy flux wrt the named parameter
    /// over the given energy range.
    virtual double energyFluxDeriv(const std::string & parName, double emin,
				   double emax, size_t npts=100) const;



    /* ----------------------- Specific methods for this class ---------------- */
    
    /* Compute the full source map, summed over all the sources, and including the spectra.
       
       model :  The vector begin filled.
     */
    void fillSummedSourceMap(std::vector<float>& model) const;

  protected:

    
  private:
    
    /// All the sources that comprise this composite
    SourceModel m_sourceModel;

    /// The parameters used to make source maps for this composite
    BinnedLikeConfig m_config;

    /// The source maps for this composite
    SourceMapCache* m_srcMapCache;

    /// Where to write the xml for this composite.  Empty string -> same file as source
    std::string m_xmlFile;

  };

} // namespace Likelihood

#endif // Likelihood_Source_h
