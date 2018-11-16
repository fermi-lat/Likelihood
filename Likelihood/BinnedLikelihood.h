/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.100 2017/10/06 01:30:59 echarles Exp $
 */

#ifndef Likelihood_BinnedLikelihood_h
#define Likelihood_BinnedLikelihood_h

#include <map>
#include <stdexcept>

#include "tip/Image.h"

#include "optimizers/dArg.h"

#include "Likelihood/Accumulator.h"
#include "Likelihood/CountsMapBase.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/Pixel.h"
#include "Likelihood/BinnedConfig.h"
#include "Likelihood/BinnedCountsCache.h"
#include "Likelihood/SourceMapCache.h"
#include "Likelihood/FileUtils.h"


namespace tip {
   class Extension;
}


namespace Likelihood {

   class Drm;
   class Drm_Cache;
   class SourceMap;
   class WeightMap;

   /*
    * @class BinnedLikelihood
    * @brief Binned version of the log-Likelihood function.
    *
    */

   class BinnedLikelihood : public LogLike {
     
   public:
     
     /* Regular c'tor 
	
	dataMap      : Observed data.  Defines the binning used for the analysis.
	observation  : Wrapper containing information about the observation
	srcMapsFile  : Name of file containing Source Maps (i.e., the output of gtsrcmaps)
	computePointSources : Pre-computes the Source Maps for point source.   
	                      Faster, but makes gtsrcmaps file much larger
	applyPsfCorrections : Apply psf integral corrections
	performConvolution  : Perform convolution with psf
	resample      : Resample input counts map for convolution
	resamp_factor : Factor used from resample maps for diffuse sources.
	minbinsz      : Minimum pixel size for rebinning fine 
     */	
     BinnedLikelihood(CountsMapBase & dataMap, 
		      const Observation & observation,
		      const std::string & srcMapsFile="",
		      bool computePointSources=true,
		      bool applyPsfCorrections=true,
		      bool performConvolution=true,
		      bool resample=true,
		      double resamp_factor=2,
		      double minbinsz=0.1);

     /* C'tor for weighted likelihood
	
	dataMap      : Observed data.  Defines the binning used for the analysis.
	weightMap    : Map with weights to use for analysis.
	observation  : Wrapper containing information about the observation
	srcMapsFile  : Name of file containing Source Maps (i.e., the output of gtsrcmaps)
	computePointSources : Pre-computes the Source Maps for point source.   
	                      Faster, but makes gtsrcmaps file much larger
	applyPsfCorrections : Apply psf integral corrections
	performConvolution  : Perform convolution with psf
	resample      : Resample input counts map for convolution
	resamp_factor : Factor used from resample maps for diffuse sources.
	minbinsz      : Minimum pixel size for rebinning fine maps 
	overwriteWeights : Overwrite the existing weights file
     */	
     BinnedLikelihood(CountsMapBase & dataMap, 
		      const ProjMap& weightMap,
		      const Observation & observation,
		      const std::string & srcMapsFile="",
		      bool computePointSources=true,
		      bool applyPsfCorrections=true,
		      bool performConvolution=true,
		      bool resample=true,
		      double resamp_factor=2,
		      double minbinsz=0.1,
		      bool overwriteWeights=false);
     
     /*  C'tor for weighted likelihood
	
	dataMap      : Observed data.  Defines the binning used for the analysis.
	observation  : Wrapper containing information about the observation
	config       : Configuration object
	srcMapsFile  : Name of file containing Source Maps (i.e., the output of gtsrcmaps)
	weightMap    : Map with weights to use for analysis.
	overwriteWeights : Overwrite the existing weights file
     */	
     BinnedLikelihood(CountsMapBase & dataMap, 
		      const Observation & observation,	
		      const BinnedLikeConfig& config,
		      const std::string & srcMapsFile = "",
		      const ProjMap* weightMap = 0,
		      bool overwriteWeights=false);
     
     /* D'tor */
     virtual ~BinnedLikelihood() throw();

     /* --------------- Simple Access functions ----------------------*/

     /// Return the binned data cache
     inline const BinnedCountsCache& dataCache() const { return m_dataCache; }

     /// Return the binned observed data
     inline const CountsMapBase & countsMap() const { return m_dataCache.countsMap(); }

     /// Return the energy bin edges
     inline const std::vector<double> & energies() const { return m_dataCache.energies(); }
  
     /// Return the configuration
     inline const BinnedLikeConfig& config() const { return m_config; }

     /// Return the number of energy bins
     inline size_t num_ebins() const { return m_dataCache.num_ebins(); }

     /// Return the number of pixels
     inline size_t num_pixels() const { return m_dataCache.num_pixels(); }

     /// Return the size of the data map
     inline size_t data_map_size() const { return m_dataCache.data_map_size(); }

     /// Return the size of the scours maps
     inline size_t source_map_size() const { return m_dataCache.source_map_size(); }
    

     /// Return the min and max energy bins to use
     std::pair<int, int> klims() const {
       return std::make_pair(static_cast<int>(m_kmin), static_cast<int>(m_kmax));
     }

     /// Return the observed counts spectrum
     inline const std::vector<double> & countsSpectrum(bool weighted=false) const { 
       return m_dataCache.countsSpectrum(weighted);
     }
     

     /// Return the weights in their original projection
     inline const ProjMap* weightMap_orig() const { return m_dataCache.weightMap_orig(); }

     /// Return the weights reprojected into counts map binning
     inline const WeightMap* weightMap() const { return m_dataCache.weightMap(); }

     /// Return the weighted counts map
     inline const CountsMapBase* weightedCounts() const { return m_dataCache.weightedCounts(); }
 
     /// Return true if the weighted counts map exists
     inline bool has_weights() const { return m_dataCache.has_weights(); }

     /// Return the name of the file with the source maps
     inline const std::string& srcMapsFile() const { return m_srcMapsFile; }

     /// Return the list of fixed sources
     inline const std::vector<std::string> & fixedSources() const { return m_fixedSources; }

     /// Return the predicted counts for all the fixed sources, summed together
     inline const std::vector<double> & fixedModelSpectrum() const { return m_fixed_counts_spec; }

     /// Return the weighted predicted counts for all the fixed sources, summed together
     inline const std::vector<double> & fixedModelSpectrum_wt() const { return m_fixed_counts_spec_wt; }
      
     /// Return the predicted counts for all the fixed sources, summed together, energy dispersion applied
     inline const std::vector<double> & fixedModelSpectrum_edisp() const { return m_fixed_counts_spec_edisp; }
      
     /// Return the weighted predicted counts for all the fixed sources, summed together, energy dispersion applied
     inline const std::vector<double> & fixedModelSpectrum_edisp_wt() const { return m_fixed_counts_spec_edisp_wt; }
        
     /// Return the fixed model counts
     inline const std::vector<double> & fixedModelCounts() const { return m_fixedModelCounts; }

     /// Check if updating the fixed model is allowed
     inline bool updateFixedWeights() const { return m_updateFixedWeights; }

     /// Set flag to enable or disable updating the fixed model 
     void setUpdateFixedWeights(bool update) {
       m_updateFixedWeights = update;
     }

     /* ----------------- Simple setter functions ------------------------ */

     /// Set the min and max energy bins to use
     void set_klims(size_t kmin, size_t kmax) {
       m_modelIsCurrent = false;
       m_kmin = kmin;
       m_kmax = kmax;
       // buildFixedModelWts();
     }
     
     /// Turn on verbose mode
     void setVerbose(bool verbose) {
       m_config.psf_integ_config().set_verbose(verbose);
       m_srcMapCache.setVerbose(verbose);
     }

     /// Turn on energy dispersion
     void set_edisp_flag(bool use_edisp) {
       if ( use_edisp == m_config.use_edisp() ) return;
       m_config.set_use_edisp(use_edisp);
       m_srcMapCache.set_edisp_flag(use_edisp);
       m_fixedModelCounts.clear();
       m_fixed_counts_spec.clear();   
       m_fixed_counts_spec_wt.clear();
       m_fixed_counts_spec_edisp.clear();
       m_fixed_counts_spec_edisp_wt.clear();
     }

     /// Set flag to use a single map for all the fixed sources
     void set_use_single_fixed_map(bool use_sfm) {
       m_config.set_use_single_fixed_map(use_sfm);
     }

     /// Set flag to same all source maps
     void set_save_all_srcmaps(bool val) {
       m_config.set_save_all_srcmaps(val);
     }

     /// Directly set the data in the counts map
     void setCountsMap(const std::vector<float> & counts);

     /// Set the weights map
     void setWeightsMap(const ProjMap* wmap);


     /* ---------------- Methods inherited from optimizers package ---------- 
	specifically optimizers::Function and optimizers::Statistic ---------------*/
     
     /// Return the current value of the log-likelihood
     /// Note that the argument is ignored.
     virtual double value(optimizers::Arg &) const;

     /// Return the current value of the log-likelihood
     virtual double value() const {
       optimizers::dArg dummy(0);
       return value(dummy);
     }

     /// Calculate the derivitives of the log-likelihood w.r.t. the free parameters
     virtual void getFreeDerivs(std::vector<double> & derivs) const;

     /// Set the parameter values
     virtual std::vector<double>::const_iterator setParamValues_(std::vector<double>::const_iterator);
     
     /// Set the free parameter values
     virtual std::vector<double>::const_iterator setFreeParamValues_(std::vector<double>::const_iterator);

     /* ---------------- Methods inherited from SourceModel ---------- */
     
     /* Create a counts map based on the current model.
	FIXME, this should make sure that map being filled 
	has the same shape as the template binning */
     virtual CountsMapBase * createCountsMap(CountsMapBase & dataMap) const {
       std::vector<float> map;
       computeModelMap(map);
       dataMap.setImage(map);
       return &dataMap;
     }

     /// Create a counts map based on the current model.
     virtual CountsMapBase * createCountsMap() const;

     /// Create the source model by reading an XML file.
     virtual void readXml(std::string xmlFile, 
			  optimizers::FunctionFactory & funcFactory,
			  bool requireExposure=true, 
			  bool addPointSources=true,
			  bool loadMaps=true);
  
     /// Add a source to the source model 
     /// This is the version inherited from SourceModel
     virtual void addSource(Source * src, bool fromClone=true, SourceMap* srcMap=0, bool loadMap=true);

     /// Add a source to the soruce model
     /// This version is particular to this class
     virtual void addSource(Source * src, BinnedLikeConfig* config, bool fromClone=true);
 
     /// Remove a source from the source model and return it
     virtual Source * deleteSource(const std::string & srcName);

     
     /// These are required since this inherits from LogLike rather than
     /// for SourceModel.  The inheritance hierarchy for this class and
     /// LogLike should be refactored.
     virtual void set_ebounds(double emin, double emax) {
       throw std::runtime_error("BinnedLikelihood::set_ebounds "
				"not implemented.");
     }
     
     virtual void unset_ebounds() {
       throw std::runtime_error("BinnedLikelihood::unset_ebounds "
				"not implemented.");
     }

     /* Synchronize parameter vector owned by this class with 
	parameters owned by the sources */
     virtual void syncParams();


     /* ---------------- Methods inherited from LogLike ---------- */

     /* Synchronize parameter vector owned by this class with 
	parameters owned by one of the sources.
	
	This actually just calls syncParams() and sets m_modelIsCurrent to false
     */
     virtual void syncSrcParams(const std::string & srcName);


     /* Return the total predicted number of counts in the ROI for a particular source

	srcName  : Name of the source
	weighted : If true returns the weights counts
     */
     virtual double NpredValue(const std::string & srcName, bool weighted=false) const;


     /* Return the total predicted number of counts in the ROI for a particular source

	This version forces the recalculation of Npred, whether the source is
	fixed or not. It is also called from buildFixedModelWts.

	It is not inherited from LogLike.
     */
     double NpredValue(const std::string & name, SourceMap & srcMap, bool weighted=false) const;

     
     /* --------------- Functions for dealing with source maps -------------- */
     
     /// Returns true if a SourceMap has been build for a source
     bool hasSourceMap(const std::string & name) const;
   
     /* Returns a reference to the SourceMap corresponding to a particular source
	
	If the source does not exist this will throw an exception
	If the sources exits, but the SourceMap does not, 
	this will create and return the SourceMap for that source */
     SourceMap & sourceMap(const std::string & name) const;
   
     /* Returns a pointer to the SourceMap corresponding to a particular source
	
	If the source does not exist this will throw an exception
	If the sources exits, but the SourceMap does not, 
	this will create and return the SourceMap for that source */
     SourceMap * getSourceMap(const std::string & srcName,
			      bool verbose=true) const;

     /* Create a new SourceMap corresponding to a particular source
	
	If the source does not exist this will throw an exception */
     SourceMap * createSourceMap(const std::string & srcName);

     
     /* Remove the SourceMap corresponding to a particular source */
     void eraseSourceMap(const std::string & srcName);
     
     /* Insert a SourceMap into this cache */
     void insertSourceMap(const std::string & srcName,
			  SourceMap& srcMap) {
       m_srcMapCache.insertSourceMap(srcName,srcMap);
     }

     /* Remove SourceMap into from cache */
     SourceMap* removeSourceMap(const std::string & srcName) {
       return m_srcMapCache.removeSourceMap(srcName);
     }

     /* Instantiate or retrieve a SourceMap for all sources and
	populate the internal map of SourceMap objects.  

	recreate  : If true new source maps will be generated for all components.  
	saveMaps  : If true the current maps will be written to the source map file. 
     */
     void loadSourceMaps(bool recreate=false, bool saveMaps=false);
     
     /* Instantiate or retrieve a SourceMap for a list of sources and
	populate the internal map of SourceMap objects.  

	recreate  : If true new source maps will be generated for all components.  
	saveMaps  : If true the current maps will be written to the source map file. 
     */
     void loadSourceMaps(const std::vector<std::string>& srcNames,
			 bool recreate=false, bool saveMaps=false);
     
     /* Instantiate or retrieve a SourceMap for a single sources and
	add it to the internal map of SourceMap objects.  

	recreate          : If true new source map will be generated.  
	buildFixedWeights : If true the fixed component weights will be recomputed
     */
     void loadSourceMap(const std::string & srcName, 
			bool recreate=false,
			bool buildFixedWeights=true);

     /* Directly set the image for a particular SourceMap
	
	name  : name of the source in question
	image : The image data.  Must be the same size as the SourceMap
     */
     void setSourceMapImage(const std::string & name,
			    const std::vector<float>& image);
  
     /* Write all of the source maps to a file 
	
	filename : The name of the file.  If empty use the current source maps file
	replace  : If true replace the SourceMaps for already in that file
     */
     void saveSourceMaps(const std::string & filename="",
			 bool replace=false,
			 bool includecountsmap=true);

     /* Write some of the source maps to a fie 
	
	filename : The name of the file.  If empty use the current source maps file
	replace  : If true replace the SourceMaps for already in that file
     */
     void saveSourceMap_partial(const std::string & filename,
				const Source& source,
				int kmin=0, int kmax=-1,
				bool replace=false);


     /* Write a single combinded source map for all the fixed sources
	
	filename : The name of the file.  If empty use the current source maps file
     */
     tip::Extension* saveTotalFixedSourceMap(bool replace = false);

     
     /* --------------- Functions for dealing with model maps -------------- */

     /* Compute a model map summing over all the sources it the model.
	Return the total number of model counts.
	
	This version uses the cached information about the fixed sources,
	and will update the m_model data member and set m_modelIsCurrent.

	weighted  : If true return the weighted counts
     */     
     double computeModelMap_internal(bool weighted=false) const;

     
     /* Compute a model map summing over all the sources it the model 
	
	This will temporarily create (and then delete ) SourceMaps for 
	sources that don't have them
      */
     void computeModelMap(std::vector<float> & modelMap,
			  bool use_mask = true) const;

     /* Compute a model map for an individual source

	This will temporarily create (and then delete ) SourceMaps for 
	sources that don't have them
     */
     void computeModelMap(const std::string& srcName, 
			  std::vector<float> & modelMap,
			  bool use_mask = true) const;

     /* Compute a model map summing over a list of sources
	
	This will temporarily create (and then delete ) SourceMaps for 
	sources that don't have them
     */
     void computeModelMap(const std::vector<std::string>& srcNames, 
			  std::vector<float> & modelMap,
			  bool use_mask = true) const;

     
     /* This function adds the Model counts for a source map
	to a vector.  it is used by the various computeModelMap 
	functions */
     void updateModelMap(std::vector<float> & modelMap, 
			 SourceMap * srcMap,
			 bool use_mask = true) const;

     /* This function adds the Model counts for a source map
	to a vector.  it is used by the various computeModelMap 
	functions */
     void updateModelMap_fromSrcMap(std::vector<float> & modelMap, 
				    const Source& src,
				    SourceMap* srcMap,
				    bool use_mask = true) const;

     /* --------------- Functions for dealing the NPreds -------------- */

     /* Return the total number of predicted counts.  
	This just calls computeModelMap.

	weighted  : If true return the weighted npred
     */
     double npred(bool weighted=false) { 
       return computeModelMap_internal(weighted);
     }
     
     /* Get or compute the npreds() vector from a particular source map.  
	If the SourceMap doesn't exist it is temporarily created and deleted

	srcName   : Name of the source in question
	npreds    : Filled with the npreds for that source

	Note, the npreds() vector is actually the differential 
	number of predicted counts at the energy bin edges.  
	It must be integrated over the enerby bin to get the number
	of predicted counts.
     */
     void getNpreds(const std::string & srcName,
		    std::vector<double> & npreds) const;
     
     /* Return the model counts spectrum from a particular source

	The counts spectra for the various source in the model are
	cached.   This will update the cache if needed. 
     */
     const std::vector<double>& modelCountsSpectrum(const std::string &srcname,
						    bool weighted=false) const;

  
     /* Return true if any fixed source has changed. 
	This will also return true if the list of fixed sources has chagned. */
     bool fixedModelUpdated() const;

     /* Compute the full source map, summed over a list of sources, and including the spectra.
	
     */
     void fillSummedSourceMap(const std::vector<std::string>& sources, std::vector<float>& model);

     /* Compute a single source map, optionally including the spectrum
	
	sourceName :  The source in question
	model   :  The vector begin filled
	mapType :  Indicates type of map begin filled
	kmin    :  Minimum energy layer
	kmax    :  Maximum energy layer	     */
     void fillSingleSourceMap(const std::string& sourceName, 
			      std::vector<float>& model,
			      FileUtils::SrcMapType& mapType,
			      int kmin=0, int kmax=-1);
   
     /* Compute the model for all the fixed source.

	process_all  : If true, build SourceMaps for all the sources.
     */
     void buildFixedModelWts(bool process_all=false);

     /* Add a source to the set of fixed sources

	This will call addSourceCounts to add the source contribution to m_fixedModelCounts
     */
     void addFixedSource(const std::string & srcName);
   
     /* Remove a source to the set of fixed sources

	This will call addSourceCounts (with subtract=true) to subtract 
	the source contribution from m_fixedModelCounts
     */     
     void deleteFixedSource(const std::string & srcName);


     /* ---------------------- other functions ------------------------ */

     /* The source component model-weighted counts spectrum, i.e.,
	based on the source model weights per pixel, this is the counts
	spectrum weighted by the probability of attributing counts in
	each pixel to this source. */
     std::vector<double> countsSpectrum(const std::string & srcName, 
					bool use_klims=true) const;

     
     /* Return true if we use energy dispersion for a particular source */
     bool use_edisp(const std::string & srcname="") const;

     /* Return the DRM (detector response matrix), building it if needed */
     Drm & drm();


     /* ---------------------- Debugging functions ---------------------- */

     /* Return the sum of the model (possibly weighted) , by actually summing the model.
      This is used to test the optimized implimentation using the npreds */
     static float npred_explicit(const BinnedLikelihood& like, 
				 bool weighted=false);
     
     /* Return the sum of the model (possibly weighted) for a single source, 
	by actually summing the model.
	This is used to test the optimized implimentation using the npreds */
     static float npred_explicit_src(const BinnedLikelihood& like, 
				     const std::string& srcName, 
				     bool weighted=false);
     
     /* Return the negative log-likelihood value (possibly weighted) by doing the full summation
	This is used to test the optimized implimentation using the npreds and 
	the filled pixels */
     static double nll_explict(const BinnedLikelihood& like, 
			       bool weighted=false);

   protected:
     
     /// Disable assignement operator
     BinnedLikelihood & operator=(const BinnedLikelihood & rhs) {
       throw std::runtime_error("Copy-assignment operator of BinnedLikelihood not implemented");
     }

     /// Disable clone function
     virtual BinnedLikelihood * clone() const {
       return new BinnedLikelihood(*this);
     }

     /// Hook to transfer information to a composite source
     virtual void initialize_composite(CompositeSource& comp) const;

     /* Save the weights SourceMap to the SourceMap file
	
	replace : if true, replace the current version 
     */
     tip::Extension* saveWeightsMap(bool replace=false) const;
     

   private:

     /* ------------- Static Utility Functions -------------------- */

     /// Calls the specturm function of a source at a given energy
     static double spectrum(const Source * src, double energy);
     
     /// Set the dimensions on a tip image
     static void setImageDimensions(tip::Image * image, long * dims);
 

   

     /* --------------- Computing Counts Spectra ------------------- */
 

     /// Integrates weights over a pixel to get the counts
     double pixelCounts(double emin, double emax, double y1, double y2, double log_ratio) const;
    
     void computeFixedCountsSpectrum();
   
     
     /* Add (or subtract) the weights for a source onto a vector 	
	This is used by several functions.

	modelWts   : The vector being added to.
	srcName    : The name of the source in question
	srcMap     : The SourceMap for the source in question
	subtract   : If true, subtract from the vector.  	
	latchParams : If true, the parameters are latched in the SourceMap
     void addSourceWts(std::vector<std::pair<double, double> > & modelWts,
		       const std::string & srcName,
		       SourceMap * srcMap=0, 
		       bool subtract=false,
		       bool latchParams=false) const;
     */

     
     /* Add (or subtract) the counts for a source onto a vector 	
	This is used by several functions.

	modelCounts: The vector being added to.
	srcName    : The name of the source in question
	srcMap     : The SourceMap for the source in question
	subtract   : If true, subtract from the vector.  	
	latchParams : If true, the parameters are latched in the SourceMap
     */     
     void addSourceCounts(std::vector<double> & modelCounts,
			  const std::string & srcName,
			  SourceMap * srcMap=0, 
			  bool subtract=false,
			  bool latchParams=false) const;

     /* Add (or subtract) the weights for a source onto the pool for fixed sources

	srcName    : The name of the source in question
	srcMap     : The SourceMap for the source in question
	subtract   : If true, subtract from the vector.  	
     */     
     void addFixedNpreds(const std::string & srcName,
			 SourceMap * srcMap=0, 
			 bool subtract=false);
     

     /* ---------------- Data Members --------------------- */

     /* ---------------- Data and binning --------------------- */

     /// This keeps track of the counts map and associated stuff
     BinnedCountsCache m_dataCache;
 
     /// Minimum and maximum energy plane indexes to use in likelihood 
     /// calculations.
     size_t m_kmin, m_kmax;     

     
     /* ------------- configuration parameters -------------------- */
     std::string m_srcMapsFile;   //! Where the SourceMaps are stored
     BinnedLikeConfig m_config;   //! All of the options

      /* ---------For keeping track of energy dispersion ----------- */

     /// Detector response matrix for energy dispersion.  Null pointer -> no energy dispersion
     Drm * m_drm;

     /* ---------------- The current model ------------------------ */

     /// The set of source maps, keyed by source name
     mutable SourceMapCache m_srcMapCache;

     /// The total model of the ROI, summed over sources, but only 
     /// for the filled pixels.
     mutable std::vector<double> m_model;
   
     /// Flag that the model is up to data
     mutable bool m_modelIsCurrent;


     /* ---------For keeping track of fixed source -------------- */

     /// List of fixed sources
     std::vector<std::string> m_fixedSources;   

     /// Summed counts for all fixed sources.
     /// This vector has the size of m_filledPixels.size()
     std::vector<double> m_fixedModelCounts;  

     /// Summed counts spectra for fixed sources
     /// This is the summed counts for all pixels for each energy bin.
     /// This vector has the size of m_energies.size() - 1
     std::vector<double> m_fixed_counts_spec;
  
     /// Summed weighted counts spectra for fixed sources
     /// This is the summed counts for all pixels for each energy bin.
     /// This vector has the size of m_energies.size() - 1
     std::vector<double> m_fixed_counts_spec_wt;
  
     /// Summed counts spectra for fixed sources with energy dispersion
     /// This is the summed counts for all pixels for each energy bin.
     /// This vector has the size of m_energies.size() - 1
     std::vector<double> m_fixed_counts_spec_edisp;
  
     /// Summed weighted counts spectra for fixed sources with energy dispersion
     /// This is the summed counts for all pixels for each energy bin.
     /// This vector has the size of m_energies.size() - 1
     std::vector<double> m_fixed_counts_spec_edisp_wt;
  
     /// Flag to allow updating of Fixed model weights
     bool m_updateFixedWeights;

};

}

#endif // Likelihood_BinnedLikelihood_h
