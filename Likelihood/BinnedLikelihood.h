/**
 * @file BinnedLikelihood.h
 * @brief Binned version of the log-likelihood function.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedLikelihood.h,v 1.84 2016/07/12 02:20:43 mdwood Exp $
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

namespace Likelihood {

   class Drm;
   class SourceMap;

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
		      double minbinsz=0.1);
     
     ///
     BinnedLikelihood(const BinnedLikelihood & other);
     
     ///
     virtual ~BinnedLikelihood() throw();

     /* --------------- Simple Access functions ----------------------*/

     /// Return the binned observed data
     const CountsMapBase & countsMap() const { return m_dataMap; }

     /// Return the weights in their original projection
     const ProjMap* weightMap() const { return m_weightMap; }

     /// Return the weights reprojected into counts map binning
     const SourceMap* weightSrcMap() const { return m_weightSrcMap; }

     /// Return the energy bin edges
     const std::vector<double> & energies() const { return m_energies; }

     /// Return the observed counts spectrum
     const std::vector<double> & countsSpectrum() const { return m_countsSpectrum; }

     /// Return the value of the flag to combined the maps for all fixed sources
     bool use_single_fixed_map() const { return m_use_single_fixed_map; }

     /// Return the min and max energy bins to use
     std::pair<int, int> klims() const {
       return std::make_pair(static_cast<int>(m_kmin), static_cast<int>(m_kmax));
     }

     /// Return the predicted counts for all the fixed sources, summed together
     const std::vector<double> & fixedModelSpectrum() const { return m_fixed_counts_spec; }
      
     /// Set flag to force updating of the fixed model 
     void setUpdateFixedWeights(bool update) {
       m_updateFixedWeights = update;
     }


     /* ----------------- Simple setter functions ------------------------ */

     /// Set the min and max energy bins to use
     void set_klims(size_t kmin, size_t kmax) {
       m_modelIsCurrent = false;
       m_kmin = kmin;
       m_kmax = kmax;
       buildFixedModelWts();
     }
     
     /// Turn on verbose mode
     void setVerbose(bool verbose) {
       m_verbose = verbose;
     }

     /// Turn on energy dispersion
     void set_edisp_flag(bool use_edisp) { 
       m_use_edisp = use_edisp;
     }

     /// Set flag to use a single map for all the fixed sources
     void set_use_single_fixed_map(bool use_sfm) {
       m_use_single_fixed_map = use_sfm;
     }

     /// Directly set the data in the counts map
     void setCountsMap(const std::vector<float> & counts);


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
     virtual void addSource(Source * src, bool fromClone=true);

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
     double NpredValue(const std::string & name, const SourceMap & srcMap, bool weighted=false) const;

     
     /* --------------- Functions for dealing with source maps -------------- */
     
     /// Returns true if a SourceMap has been build for a source
     bool hasSourceMap(const std::string & name) const;
   
     /* Returns a reference to the SourceMap corresponding to a particular source
	
	If the source does not exist this will throw an exception
	If the sources exits, but the SourceMap does not, 
	this will create and return the SourceMap for that source */
     const SourceMap & sourceMap(const std::string & name) const;
   
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
			 bool replace=false);
     
     /* --------------- Functions for dealing with model maps -------------- */

     /* FIXME: sort out the version of computeModelMap */

     /* Compute a model map summing over all the sources it the model.
	Return the total number of model counts.
	
	This version uses the cached information about the fixed sources

	weighted  : If true return the weighted counts
     */     
     double computeModelMap(bool weighted=false) const;

     
     /* Compute a model map summing over all the sources it the model 
	
	This will temporarily create (and then delete ) SourceMaps for 
	sources that don't have them
      */
     void computeModelMap(std::vector<float> & modelMap) const;

     /* Compute a model map for an individual source

	This will temporarily create (and then delete ) SourceMaps for 
	sources that don't have them
     */
     void computeModelMap(const std::string& srcName, 
			  std::vector<float> & modelMap) const;

     /* Compute a model map summing over a list of sources
	
	This will temporarily create (and then delete ) SourceMaps for 
	sources that don't have them
     */
     void computeModelMap(const std::vector<std::string>& srcNames, 
			  std::vector<float> & modelMap) const;

     
     /* This function adds the Model counts for a source map
	to a vector.  it is used by the various computeModelMap 
	functions */
     void updateModelMap(std::vector<float> & modeMap, 
			 const SourceMap * srcMap) const;


     /* --------------- Functions for dealing the NPreds -------------- */

     /* Return the total number of predicted counts.  
	This just calls computeModelMap.

	weighted  : If true return the weighted npred
     */
     double npred(bool weighted=false) { 
       return computeModelMap(weighted);
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
     const std::vector<double>& modelCountsSpectrum(const std::string &srcname) const;

  
     /* Return true if any fixed source has changed. 
	This will also return true if the list of fixed sources has chagned. */
     bool fixedModelUpdated() const;

   
     /* Compute the model for all the fixed source.

	process_all  : If true, build SourceMaps for all the sources.
     */
     void buildFixedModelWts(bool process_all=false);

     /* Add a source to the set of fixed sources

	This will add an entry in the map of m_fixedModelNpreds,
	subtract the source contribution to m_fixedNpreds,
	and call addSourceWts to add the source contribution to m_fixedModelWts
     */
     void addFixedSource(const std::string & srcName);
   
     /* Remove a source to the set of fixed sources

	This will remove an entry in the map of m_fixedModelNpreds,
	subtract the source contribution to m_fixedNpreds,
	and call addSourceWts (with subtract=true) to subtract 
	the source contribution from m_fixedModelWts
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

   protected:
     
     /// Disable assignement operator
     BinnedLikelihood & operator=(const BinnedLikelihood & rhs) {
       throw std::runtime_error("Copy-assignment operator not implemented");
     }

     /// Disable clone function
     virtual BinnedLikelihood * clone() const {
       return new BinnedLikelihood(*this);
     }
     

     /* Save the weights SourceMap to the SourceMap file
	
	replace : if true, replace the current version 
     */
     void saveWeightsMap(bool replace=false) const;
     
     /// Fill the map of weighted counts
     void fillWeightedCounts();


   private:

     /* ------------- Static Utility Functions -------------------- */

     /// Calls the specturm function of a source at a given energy
     static double spectrum(const Source * src, double energy);
     
     /// Integrates weights over a pixel to get the counts
     static double pixelCounts(double emin, double emax, double y1, double y2);

     /// Set the dimensions on a tip image
     static void setImageDimensions(tip::Image * image, long * dims);
 

     /* ------------- Dealing with SourceMaps -------------------- */

     bool fileHasSourceMap(const std::string & srcName, 
			   const std::string & fitsFile) const;
     
     void replaceSourceMap(const std::string & srcName, 
			   const std::string & fitsFile) const;
     
     void replaceSourceMap_wcs(const SourceMap& srcMap, 
			       const std::string & fitsFile) const;
     
     void replaceSourceMap_healpix(const SourceMap& srcMap, 
				   const std::string & fitsFile) const;
     
     void appendSourceMap(const std::string & srcName, 
			  const std::string & fitsFile,
			  bool isWeights = false) const;
     
     void appendSourceMap_wcs(const SourceMap& srcMap,
			      const std::string & fitsFile,
			      bool isWeights = false) const;
     
     void appendSourceMap_healpix(const SourceMap& srcMap, 
				  const std::string & fitsFile,
				  bool isWeights = false) const;
     

     /* --------------- Computing Counts Spectra ------------------- */

     void computeCountsSpectrum();
     
     void computeCountsSpectrum_wcs();
     
     void computeCountsSpectrum_healpix();
   
     void computeFixedCountsSpectrum();
   
     
     /* Add (or subtract) the weights for a source onto a vector 	
	This is used by several functions.

	modelWts   : The vector being added to.
	srcName    : The name of the source in question
	srcMap     : The SourceMap for the source in question
	subtract   : If true, subtract from the vector.  	
     */     
     void addSourceWts(std::vector<std::pair<double, double> > & modelWts,
		       const std::string & srcName,
		       const SourceMap * srcMap=0, 
		       bool subtract=false) const;
     
     /* Fills the m_filledPixels data member with only the pixels with data counts */
     void identifyFilledPixels();
      
     
     /* --------------- Dealing with Energy Dispersion ------------------- */

     void updateCorrectionFactors(const std::string & srcName, const SourceMap & sourceMap) const;
     
     void edisp_correction_factors(const std::string & srcName,
				   const std::vector<double> & true_counts_spec,
				   std::vector<double> &);
     
    

     /* ---------------- Data Members --------------------- */

     /* ---------------- Data and binning --------------------- */

     /// The observed data.  Used to provide the binning.
     CountsMapBase& m_dataMap;

     /// Set of the pixels in the ROI
     const std::vector<Pixel> & m_pixels;

     /// Energy bin edges
     std::vector<double> m_energies;

     /// Minimum and maximum energy plane indexes to use in likelihood 
     /// calculations.
     size_t m_kmin, m_kmax;     

     /// The observed counts spectrum.  
     /// I.e., the data summed of the ROI for each energy bin
     std::vector<double> m_countsSpectrum;

     /// These are the indices of the pixels with counts     
     /// This is used to speed up the evaluation of the log-likelihood
     std::vector<unsigned int> m_filledPixels;

     /* ---------------- The current model ------------------------ */

     /// The set of source maps, keyed by source name
     mutable std::map<std::string, SourceMap *> m_srcMaps;

     /// The total model of the ROI, summed over sources. 
     mutable std::vector<double> m_model;
   
     /// Flag that the model is up to data
     mutable bool m_modelIsCurrent;


     /* ---------------- Stuff for weighted likelihood --------------- */
     
     /// Weights map.  Null ptr -> don't use weights 
     const ProjMap* m_weightMap;
     
     /// Weights map reprojected into counts map binning
     SourceMap* m_weightSrcMap;
     
     /// Vector of the weighted counts (FIXME, replace with CountsMapBase*)
     std::vector<float> m_weightedCounts;




     /* ---------For keeping track of fixed source -------------- */

     /// List of fixed sources
     std::vector<std::string> m_fixedSources;   

     /// Map of model parameters, to be used to determine if fixed
     /// sources have changed parameter values.
     std::map<std::string, std::vector<double> > m_modelPars;
     
     /// FIXME, what are the differences between the next three
     /// Summed npred values at each energy boundary value for fixed sources.
     std::vector<double> m_fixedNpreds;
     /// Summed weights for all fixed sources
     std::vector<std::pair<double, double> > m_fixedModelWts;  
     /// Summed counts spectra for fixed sources
     mutable std::vector<double> m_fixed_counts_spec;
  
     /// Maps of Npreds for fixed sources, keyed by source name
     std::map<std::string, double> m_fixedModelNpreds;         //! Npreds 
     std::map<std::string, double> m_fixedModelWeightedNpreds; //! Weights NPreds 

     /// Flag to force updating of Fixed model weights
     bool m_updateFixedWeights;

     /* ---------For keeping track of energy dispersion ----------- */

     /// Detector response matrix for energy dispersion.  Null pointer -> no energy dispersion
     Drm * m_drm;

     /// Maps of true and measured energy spectra for sources 
     mutable std::map<std::string, std::vector<double> > m_true_counts;
     mutable std::map<std::string, std::vector<double> > m_meas_counts;

     /// Maps of references pixels to use when computing energy dispersion
     std::map<std::string, std::map<size_t, size_t> > m_krefs;


     /* ------------- configuration parameters -------------------- */
     std::string m_srcMapsFile;   //! Where the SourceMaps are stored
     bool m_computePointSources;  //! Pre-compute and store SourceMaps for point sources
     bool m_applyPsfCorrections;  //! Apply PSF integral corrections
     bool m_performConvolution;   //! Do PSF convolution when making SourceMaps
     bool m_resample;             //! turn on Resample when projecting Diffuse source maps
     double m_resamp_factor;      //! Resample factor for projecting Diffuse source maps
     double m_minbinsz;           //! Minimum pixel size for rebinning fine maps 
     bool m_verbose;              //! Turn on verbose output
     bool m_use_edisp;            //! Use energy dispersion
     bool m_use_single_fixed_map; //! Use a single model for all fixed components


     /// Accumulators for derivatives, FIXME not sure why these are data members.
     mutable std::map<long, Kahan_Accumulator> m_posDerivs;
     mutable std::map<long, Kahan_Accumulator> m_negDerivs;
     





};

}

#endif // Likelihood_BinnedLikelihood_h
