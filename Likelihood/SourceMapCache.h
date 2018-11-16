/**
 * @file SourceMapCache.h
 * @brief Functionality to deal with source maps extracted from BinnedLikelihood
 * @author E. Charles, (from BinnedLikelihood by J. Chiang)
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMapCache.h,v 1.6 2017/09/29 01:38:06 echarles Exp $
 */

#ifndef Likelihood_SourceMapCache_h
#define Likelihood_SourceMapCache_h

#include <map>
#include <vector>
#include <string>

#include "Likelihood/CountsMapBase.h"
#include "Likelihood/BinnedConfig.h"
#include "Likelihood/FileUtils.h"
#include "Likelihood/Accumulator.h"

namespace tip {
   class Extension;
}

namespace Likelihood {

   class BinnedCountsCache;
   class BinnedLikeConfig;
   class CountsMapBase;
   class Drm;
   class Drm_Cache;
   class SourceMap;
   class ProjMap;
   class WeightMap;

   /*
    * @class BinnedLikelihood
    * @brief Binned version of the log-Likelihood function.
    *
    */

   class SourceMapCache {
     
   public:

     /* ------------- Static Utility Functions -------------------- */

     /* Get the model counts contribution to a particular pixel 

	srcMap     : The SourceMap for the source in question
	spec       : The specturm (or spectral derivative) for the source in question
	jmin       : The unrolled index for the source map at the lower edge of the energy bin
	jmax       : The unrolled index for the source map at the upper edge of the energy bin
	kref       : The energy bin
	emin       : The energy at the low edge of the bin
	emax       : The energy at the high edge of the bin
	log_ratio  : The log of the energy ratios in the bin

	returns the contribution
     */
     static double model_counts_contribution(SourceMap& srcMap,
					     const std::vector<double> & spec,
					     size_t jmin,
					     size_t jmax,
					     size_t kref,
					     const double& emin,
					     const double& emax,
					     const double& log_ratio);

    /* Get the model counts contribution to a energy layer

	npred_vals     : The total npred for that energy layer
	npred_weights  : The weight factors for that energy layer
	spec       : The specturm (or spectral derivative) for the source in question
	kref       : The energy bin
	emin       : The energy at the low edge of the bin
	emax       : The energy at the high edge of the bin
	log_ratio  : The log of the energy ratios in the bin
	counts     : Filled with the counts contribution
	counts_wt  : Filled with the weighted counts contribution
     */
     static void npred_contribution(const std::vector<double>& npred_vals,
				    const std::vector<std::pair<double,double> >& npred_weights,
				    const std::vector<double>& spec,
				    size_t kref,
				    const double& emin,
				    const double& emax,
				    const double& log_ratio,
				    double& counts,
				    double& counts_wt);


     /* Add (or subtract) the counts for a source onto a vector 	
	This is used by several functions.

	modelCounts: The vector being added to.
	srcMap     : The SourceMap for the source in question
	npix       : Number of pixel in the map, used for indexing
	filledPixels : Vector with the indices of the filled pixels,
	dataCache  : Object with info about the binning
	use_edisp_val : Apply the energy dispersion
	subtract   : If true, subtract from the vector.  	
     */     

     static void addSourceCounts_static(std::vector<double> & modelCounts,
					SourceMap& srcMap,
					const BinnedCountsCache& dataCache,
					bool use_edisp_val,
					bool subtract);
 
     /* Add (or subtract) the counts for a fxied source onto the vectors that 
	collect that info for fixed sources.
	This is used by several functions.

	fixed_counts_spec : The vector being added to.
	fixed_counts_spec_wt : The vector being added to.
	fixed_counts_spec_edisp : The vector being added to.
	fixed_counts_spec_edisp_wt : The vector being added to.
	srcMap     : The SourceMap for the source in question
	dataCache  : Object with info about the binning
	use_edisp_val : Apply the energy dispersion
	subtract   : If true, subtract from the vector.  	
     */     
     static void addFixedNpreds_static(std::vector<double>& fixed_counts_spec,
				       std::vector<double>& fixed_counts_spec_wt,
				       std::vector<double>& fixed_counts_spec_edisp,
				       std::vector<double>& fixed_counts_spec_edisp_wt,
				       SourceMap& srcMap,
				       const BinnedCountsCache& dataCache,
				       bool use_edisp_val,
				       bool subtract);

#ifndef SWIG
     /* Add (or subtract) the contributions to the derivatives from a single source.

	posDerivs : The vector being added to.
	negDerivs : The vector being added to.
	freeIndex : The overall index of the first parameter for this source
	srcMap     : The SourceMap for the source in question
	dataCache  : Object with info about the binning
	use_edisp_val : Apply the energy dispersion
     */     
     static void addFreeDerivs_static(std::vector<Kahan_Accumulator>& posDerivs,
				      std::vector<Kahan_Accumulator>& negDerivs,
				      long freeIndex,
				      SourceMap& srcMap,
				      const std::vector<float> & data,
				      const std::vector<double>& model, 
				      const BinnedCountsCache& dataCache,
				      bool use_edisp_val,
				      size_t kmin, size_t kmax);
#endif //SWIG

          
    /* Add (or subtract) the counts for a source onto a vector.
       This version loops over all the pixels, not just the filled ones.

	modelMap   : The vector being added to.
	srcMap     : The SourceMap for the source in question
	dataCache  : Object with info about the binning
	use_mask   : Use the Weights mask
	use_edisp_val : Apply the energy dispersion
     */     
     static void updateModelMap_static(std::vector<float> & modelMap,
				       SourceMap& srcMap,
				       const BinnedCountsCache& dataCache,				       
				       bool use_mask,
				       bool use_edisp_val);

     

   public:
     
     /* Regular c'tor 
	
	dataCache     : Observed data.  Defines the binning used for the analysis.
	observation   : Wrapper containing information about the observation
	srcMapsFile   : Name of file containing Source Maps (i.e., the output of gtsrcmaps)
	config        : Object with configuration parameters
	drm           : The detector response matrix, NULL if energy dispersion is off for this source.
     */	
     SourceMapCache(const BinnedCountsCache& dataCache,
		    const Observation & observation,
		    const std::string & srcMapsFile,
		    const BinnedLikeConfig& config,
		    const Drm* drm = 0);


     /// Copy c'tor
     SourceMapCache(const SourceMapCache& other);

     /// Clone function
     virtual SourceMapCache* clone() const {
       return new SourceMapCache(*this);
     }


     ///
     virtual ~SourceMapCache();

     /* --------------- Simple Access functions ----------------------*/

     /// Return the binned observed data
     const BinnedCountsCache& dataCache() const { return m_dataCache; }
  
     /// Return the configuration
     const BinnedLikeConfig& config() const { return m_config; }
    
     /// Return the name of the file with the source maps
     inline const std::string& srcMapsFile() const { return m_srcMapsFile; }

     /// Return the number of cached sources
     inline size_t n_srcs() const { return m_srcMaps.size(); }


     /* ----------------- Simple setter functions ------------------------ */
      
     /// Turn on verbose mode
     void setVerbose(bool verbose) {
       m_config.psf_integ_config().set_verbose(verbose);
     }
     
     /// Turn on energy dispersion
     void set_edisp_flag(bool use_edisp) { 
       m_config.set_use_edisp(use_edisp);
     }


     /* ---------------- Methods Used by SourceModel ---------- */
     
     /* Create a counts map based on the current model.
	FIXME, this should make sure that map being filled 
	has the same shape as the template binning */
     virtual CountsMapBase * createCountsMap(const std::vector<const Source*>& srcs, 
					     CountsMapBase & dataMap) const;

     /// Create a counts map based on the current model.
     virtual CountsMapBase * createCountsMap(const std::vector<const Source*>& srcs) const;


     /* ---------------- Methods Used by LogLike ---------- */

     /* Return the total predicted number of counts in the ROI for a particular source

	srcName  : Name of the source
	weighted : If true returns the weights counts
     */
     double NpredValue(const Source& src, size_t kmin, size_t kmax, bool weighted=false) const;


     /* Return the total predicted number of counts in the ROI for a particular source

	This version forces the recalculation of Npred, whether the source is
	fixed or not. It is also called from buildFixedModelWts.

	It is not inherited from LogLike.
     */
     double NpredValue(const Source& src, SourceMap & srcMap, size_t kmin, size_t kmax, bool weighted=false) const;

     
     /* --------------- Functions for dealing with source maps -------------- */
     
     /// Returns true if a SourceMap has been build for a source
     bool hasSourceMap(const std::string & name) const;
      
     /* Returns a pointer to the SourceMap corresponding to a particular source
	
	If the source does not exist this will throw an exception
	If the sources exits, but the SourceMap does not, 
	this will create and return the SourceMap for that source */
     SourceMap * getSourceMap(const Source& src,			      
			      bool verbose=true,
			      const BinnedLikeConfig* config = 0) const;

     /* Create a new SourceMap corresponding to a particular source
	
	If the source does not exist this will throw an exception */
     SourceMap * createSourceMap(const Source& src, const BinnedLikeConfig* config = 0) const;

     
     /* Remove the SourceMap corresponding to a particular source */
     void eraseSourceMap(const std::string & srcName);


     /* Insert a SourceMap into this cache */
     void insertSourceMap(const std::string & srcName,
			  SourceMap& srcMap);

     /* Remove SourceMap into from cache */
     SourceMap* removeSourceMap(const std::string & srcName);


     /* Instantiate or retrieve a SourceMap for a list of sources and
	populate the internal map of SourceMap objects.  

	recreate  : If true new source maps will be generated for all components.  
	saveMaps  : If true the current maps will be written to the source map file. 
     */
     void loadSourceMaps(const std::vector<const Source*>& srcs,
			 bool recreate=false, bool saveMaps=false);
     
     /* Instantiate or retrieve a SourceMap for a single sources and
	add it to the internal map of SourceMap objects.  

	recreate          : If true new source map will be generated.  
	buildFixedWeights : If true the fixed component weights will be recomputed
     */
     void loadSourceMap(const Source& src,
			bool recreate=false,
			const BinnedLikeConfig* config=0);

     /* Directly set the image for a particular SourceMap
	
	name  : name of the source in question
	image : The image data.  Must be the same size as the SourceMap
     */
     void setSourceMapImage(const Source& src,
			    const std::vector<float>& image);
  
     /* Write source maps to a file 
	
	filename : The name of the file.  If empty use the current source maps file
	srcs     : The sources to write
	replace  : If true replace the SourceMaps for already in that file
     */
     void saveSourceMaps(const std::string & filename,
			 const std::vector<const Source*>& srcs,
			 bool replace=false);

     
      /* Write partial source maps to a file 
	
	filename : The name of the file.  If empty use the current source maps file
	source   : The source to write
	replace  : If true replace the SourceMaps for already in that file
	
     */
     void saveSourceMap_partial(const std::string & filename,
				const Source& source,
				int kmin = 0, int kmax = -1,
				bool replace=false);

     
     /* --------------- Functions for dealing with model maps -------------- */

      

     /* Compute a model map for an individual source

	This will temporarily create (and then delete ) SourceMaps for 
	sources that don't have them
     */
     void computeModelMap(const Source& src, 
			  std::vector<float> & modelMap,
			  bool use_mask = true) const;

     /* Compute a model map summing over a list of sources
	
	This will temporarily create (and then delete ) SourceMaps for 
	sources that don't have them
     */
     void computeModelMap(const std::vector<const Source*>& srcs, 
			  std::vector<float> & modelMap,
			  bool use_mask = true) const;

     
     /* This function adds the Model counts for a source map
	to a vector.  it is used by the various computeModelMap 
	functions */
     void updateModelMap(std::vector<float> & modelMap, 
			 const Source& src, 
			 SourceMap* srcMap,
			 bool use_mask = true) const;

     /* This function adds the Model counts for a source map
	to a vector.  it is used by the various computeModelMap 
	functions */
     void updateModelMap_fromSrcMap(std::vector<float> & modelMap, 
				    const Source& src, 
				    SourceMap* srcMap,
				    bool use_mask = true) const;
     


     /* --------------- Functions for dealing the NPreds -------------- */

     
     /* Get or compute the npreds() vector from a particular source map.  
	If the SourceMap doesn't exist it is temporarily created and deleted

	srcName   : Name of the source in question
	npreds    : Filled with the npreds for that source

	Note, the npreds() vector is actually the differential 
	number of predicted counts at the energy bin edges.  
	It must be integrated over the enerby bin to get the number
	of predicted counts.
     */
     void getNpreds(const Source& src,
		    std::vector<double> & npreds) const;
     
     /* Return the model counts spectrum from a particular source

	The counts spectra for the various source in the model are
	cached.   This will update the cache if needed. 
     */
     const std::vector<double>& modelCountsSpectrum(const Source& src,
						    bool weighted=false) const;
 

     /* Set the weights map */
     void setWeightsMap(const ProjMap* wmap);

 
     /* Compute the full model for all the fixed source.

	sources :  All the source to add to the summed map
	model   :  The vector begin filled
	kmin    :  Minimum energy layer
	kmax    :  Maximum energy layer	
     */
     void fillSummedSourceMap(const std::vector<const Source*>& sources, 
			      std::vector<float>& model, 
			      int kmin=0, int kmax=-1);
   

     /* Compute the source map for a single source.

	sourceName :  The source in question
	model   :  The vector begin filled
	kmin    :  Minimum energy layer
	kmax    :  Maximum energy layer	
     */
     void fillSingleSourceMap(const Source& src,
			      std::vector<float>& model, 
			      FileUtils::SrcMapType& mapType,
			      int kmin=0, int kmax=-1) const;


     /* ---------------------- other functions ------------------------ */

   
     /* Return true if we use energy dispersion for a particular source */
     bool use_edisp(const Source* src = 0) const;


     /* --------------------- Debugging -------------------- */
     size_t memory_size() const;

   protected:
     
     /// Disable assignement operator
     SourceMapCache & operator=(const SourceMapCache& rhs) {
       throw std::runtime_error("Copy-assignment operator of SourceMapCache not implemented");
     }

     
   private:


     /* ------------- Dealing with SourceMaps -------------------- */

    
     tip::Extension* replaceSourceMap(const Source & src, 
				      const std::string & fitsFile) const;
     
  
     tip::Extension* appendSourceMap(const Source & src, 
				     const std::string & fitsFile) const;
     
     tip::Extension* replaceSourceMap_partial(const Source & src, 
					      const std::string & fitsFile,
					      int kmin=0, int kmax=-1) const;
       
     tip::Extension* appendSourceMap_partial(const Source & src, 
					     const std::string & fitsFile,
					     int kmin=0, int kmax=-1) const;
    

     /* --------------- Computing Counts Spectra ------------------- */
 

      void computeFixedCountsSpectrum();
   
     
     /* Add (or subtract) the weights for a source onto a vector 	
	This is used by several functions.

	modelWts   : The vector being added to.
	srcName    : The name of the source in question
	srcMap     : The SourceMap for the source in question
	subtract   : If true, subtract from the vector.  	
	latchParams : If true, the parameters are latched in the SourceMap
     void addSourceWts(std::vector<std::pair<double, double> > & modelWts,
		       const Source & src,
		       SourceMap * srcMap=0, 
		       bool subtract=false,
		       bool latchParams=false) const;
     */
     
     
     /* --------------- Dealing with Energy Dispersion ------------------- */
     void updateCorrectionFactors(const Source & src, SourceMap & sourceMap) const;
     
     
    

     /* ---------------- Data Members --------------------- */

     /* ---------------- Data and binning --------------------- */

     /// The observed data.  Used to provide the binning.
     const BinnedCountsCache& m_dataCache;

     /// The obvseration
     const Observation& m_observation;

     /* ---------------- The current model ------------------------ */

     /// The set of source maps, keyed by source name
     mutable std::map<std::string, SourceMap *> m_srcMaps;


     /* ---------For keeping track of energy dispersion ----------- */

     /// Detector response matrix for energy dispersion.  Null pointer -> no energy dispersion
     const Drm * m_drm;

     /* ------------- configuration parameters -------------------- */
     std::string m_srcMapsFile;   //! Where the SourceMaps are stored
     BinnedLikeConfig m_config;   //! All of the options

};

}

#endif // Likelihood_BinnedLikelihood_h
