/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.71 2016/09/29 00:25:55 echarles Exp $
 */

#ifndef Likelihood_SourceMap_h
#define Likelihood_SourceMap_h

#include <vector>
#include <map>

#include "st_facilities/libStApiExports.h"

#include "Likelihood/BinnedConfig.h"
#include "Likelihood/FileUtils.h"
#include "Likelihood/SparseVector.h"

namespace astro {
   class HealpixProj;
}

namespace st_stream {
   class StreamFormatter;
}

namespace Likelihood {
  
   class BinnedCountsCache;
   class PsfIntegConfig;
   class Drm;
   class Drm_Cache;
   class DiffuseSource;
   class MeanPsf;
   class PointSource;
   class Source;
   class WcsMap2;
   class WeightMap;

/*
 * @class SourceMap
 *
 */

#ifdef SWIG
class SourceMap {
#else
class  SCIENCETOOLS_API SourceMap {
#endif

public:

  static void fill_sparse_model(const std::vector<float>& vect,
				SparseVector<float>& sparse);
  
  static void fill_full_model(const SparseVector<float>& sparse,
			      std::vector<float>& vect);
  
public:

   /* Standar c'tor 
      
      src           : The source making this source map for
      dataCache     : The counts map used as a template for the binning
      observation   : Object with data about the observation
      psf_config    : Object with PSF integration parameters
      drm           : The detector response matrix, NULL if energy dispersion is off for this source.
      weights       : A weights map, for weighted likelihood analysis.  Null -> no weighting.
      save_model    : Flag to indicate that we should save the model (e.g., when the source is fixed)
                      This avoid reload the model map if the source is subsequently freed.
    */
   SourceMap(const Source & src, 
	     const BinnedCountsCache * dataCache,
             const Observation & observation,
	     const PsfIntegConfig & psf_config, 
	     const Drm* drm = 0,
	     const WeightMap* weights = 0,
	     bool save_model = false);


   /* C'tor to re-read this from a file
      
      sourceMapsFile : The path to the file we are reading from
      src           : The source in question
      dataMap       : The counts map used as a template for the binning
      observation   : Object with data about the observation
      weights       : A weights map, for weighted likelihood analysis.  Null -> no weighting.
      drm           : The detector response matrix, NULL if energy dispersion is off for this source.
      save_model    : Flag to indicate that we should save the model (e.g., when the source is fixed)
                      This avoid reload the model map if the source is subsequently freed.
   */
   SourceMap(const std::string & sourceMapsFile,
             const Source & src,
	     const BinnedCountsCache * dataCache,
	     const Observation & observation,
	     const WeightMap* weights = 0,
	     const Drm* drm = 0,
	     bool save_model = false);

   /* Copy c'tor */
   SourceMap(const SourceMap& other);

   /* d'tor, does clean up */
   ~SourceMap();

   /* --------------- Operators ----------------------------------- */

   /* Get access to the model value for a particular pixel */
   float operator[](size_t idx) const;

   /* --------------- Simple Access functions ----------------------*/

   /* The source in question */
   inline const Source* src() const { return m_src; }

   /* The name of the source */
   inline const std::string & name() const { return m_name; }

   /* The source type, either 'Point' or 'Diffuse' */
   inline const std::string & srcType() const { return m_srcType; }

   /* The parameters for the PSF Integration */
   inline const PsfIntegConfig& psf_config() const { return m_psf_config; }

   /* The detector response matrix */
   inline const Drm* drm() const { return m_drm; }

   /* The weights for the weighted log-likelihood.  Null-> no weights */
   inline const WeightMap* weights() const { return m_weights; } 

   /* How the source map is stored */
   inline FileUtils::SrcMapType mapType() const { return m_mapType; }

   /* Flag to indicat that we should save the model */
   inline bool save_model() const { return m_save_model; }

   /* --------------- Class Methods ----------------------*/

   /* Clear out the model, but save the cached data.
      This is useful for dealing with fixed sources, as it frees
      up a lot of memory.       
   */
   void clear_model(bool force=false) {
     if ( force || !m_save_model ) {
       m_model.clear();
       m_sparseModel.clear();
     }
   }      

   /* Set the source associated with this source map, this is useful for
      functions that add & remove source from the source model */
   void setSource(const Source& src);

   /* Extract a vector of spectral normalization values from a Source object
      and latch it in this class.
      
      energies:  The energies at which to evalute the spectrum
      latch_params : If true, the model parameters are latched
   */   
   void setSpectralValues(const std::vector<double>& energies,
			  bool latch_params = false);

   /* Extract the spectral derivatives from the Source object and 
      latch them in this class.
      
      energies:  The energies at which to evalute the derivatives
      paraNames: The names of the params w.r.t. which to evaluate the derivaties
    */
   void setSpectralDerivs(const std::vector<double>& energies,
			  const std::vector<std::string>& paramNames);

   
   /* Test to see if the spectrum has changes w.r.t. the cached values */
   bool spectrum_changed() const;


   /* Compute the total model counts summed between two energy bins *
      This uses the Drm_Cache, so that must be up to date.

      kmin         : index of lowest energy bin
      kman         : index of highest energy bin
      use_edisp    : return model counts with energy dispersion applied
      use_weighted : return weighted model counts */
   double summed_counts(size_t kmin, size_t kmax,
			bool use_edisp = true,
			bool use_weighted = false);

   /* --------------- Cached Data ----------------------*/
 
   /* These functions return cached data.  
      This might not be calculated yet, or might be out of date.
      Use with caution.
      The version of the functions without cached_ will always return 
      up-to-date information. */

   /* The source map model.  This must be multiplied by the spectrum for each pixel 
      and integrated over the energy bin to obtain the predicted counts */
   inline const std::vector<float> & cached_model() const { return m_model; }
 

   /* The sparse version of source map model.  This must be multiplied by the spectrum for each pixel 
      and integrated over the energy bin to obtain the predicted counts */
   inline const SparseVector<float> & cached_sparse_model() const { return m_sparseModel; }
   

   /* These are the derivatives of the 'spectrum' values.  I.e., the derivatives evaluated 
      at the energy pointsThese are the spectral derivatives.  */
   inline const std::vector<std::vector<double> >& cached_specDerivs() const { return m_derivs; }

   /* This contains both the convolved and un-convolved counts spectra */
   inline const Drm_Cache* cached_drm_Cache() const { return m_drm_cache; }


   /* --------------- Cached Data ----------------------*/

   /* These functions will either return the cached value or compute it if needed or if force = true */

   /* The source map model.  This must be multiplied by the spectrum for each pixel 
      and integrated over the energy bin to obtain the predicted counts */
   const std::vector<float> & model(bool force=false);

   /* These are the 'spectrum' values, I.e., the spectrum evaluated at the energy points */
   const std::vector<double> & specVals(bool force=false);

   /* These are the derivatives of the 'spectrum' values.  I.e., the derivatives evaluated 
      at the energy pointsThese are the spectral derivatives.  */
   const std::vector<std::vector<double> >& specDerivs(const std::vector<std::string>& paramNames, bool force = false);
   
   /* These are the npreds, i.e., the model summed over each energy plane. */
   const std::vector<double> & npreds(bool force=false);

   /* These are the npred weights, i.e., the weights to apply to the npreds
      to correctly reproduce the weighted counts */
   const std::vector<std::pair<double,double> > & npred_weights(bool force=false);
   
   /* This contains both the convolved and un-convolved counts spectra */
   const Drm_Cache* drm_cache(bool force=false);

   /* Add the model to an extermal vector */
   void addToVector(std::vector<float>& vect, bool includeSpec = false);
   
   /* Subtract the model from an extermal vector */
   void subtractFromVector(std::vector<float>& vect, bool includeSpec = false);

   /// Explicitly set the image data
   void setImage(const std::vector<float>& model);

   /// Explicitly set the weights data
   void setWeights(const WeightMap* weights);

   /// Explicity set the flag to save the model
   inline void setSaveModel(bool val) { m_save_model = val; }

   /// Set the filename (e.g., b/c we are writing the source map)
   inline void setFilename(const std::string& filename) { m_filename = filename; }


   /* --------------------- Debugging -------------------- */
   size_t memory_size() const;
   
   void test_sparse(const std::string& prefix) const;


protected:

   /* Read the model from a file */
   int readModel(const std::string& sourceMapFile);

   /* Read an image from a FITS file */
   int readImage(const std::string& sourceMapFile);

   /* Read HEALPix data from a table in a FITS file */
   int readTable_healpix(const std::string& sourceMapFile);

   /* Make the model */
   int make_model();

   /* Sparsify the full model */
   void sparsify_model(bool clearFull = true);

   /* Fill the full model */
   void expand_model(bool clearSparse = true);


private:

   /* Implementation of addToVector for full-model storage */
   void addToVector_full(std::vector<float>& vect, bool includeSpec=false) const;
 
   /* Implementation of addToVector for sparse-model storage */
   void addToVector_sparse(std::vector<float>& vect, bool includeSpec=false) const;

   /* Implementation of subtractFromVector for full-model storage */
   void subtractFromVector_full(std::vector<float>& vect, bool includeSpec=false) const;
 
   /* Implementation of subtractFromVector for sparse-model storage */
   void subtractFromVector_sparse(std::vector<float>& vect, bool includeSpec=false) const;
 
   /* Compute the NPreds */
   void computeNpredArray();

   /* Adjust the source map for a phased exposure correction       
      This is just multiplying the source map by the phased exposure map projected into the counts map frame
   */
   void applyPhasedExposureMap();


   /* Get the value from the sparse model */
   float find_value(size_t idx) const {
     return m_sparseModel[idx];
   }


   /* ---------------- Data Members --------------------- */

   /// The source in question.
   const Source* m_src;   

   /// Source name.
   const std::string m_name;

   /// The name of the file where the source map is stored
   std::string m_filename;

   /// This is either "Diffuse" or "Point"
   const std::string m_srcType;

   /// Pointer to the data map that is the template for the binning
   const BinnedCountsCache * m_dataCache;
   
   /// Object holding data about the observation
   const Observation & m_observation;

   /// PSF for this particular source ( NULL for diffuse sources ) 
   MeanPsf* m_meanPsf; 

   /// For Progress messages
   st_stream::StreamFormatter * m_formatter;

   /// Options for treatment of PSF
   PsfIntegConfig m_psf_config;

   /// The detector response matrix.  Null -> ignore energy disperison
   const Drm* m_drm;

   /// A weights map, for weighted likelihood analysis.  Null -> no weighting.   
   const WeightMap* m_weights;

   /// Flag to indicate that we should save the model
   bool m_save_model;

   /// This is the 'source map' data. 
   /// It is not the counts map, but rather the coefficients
   /// That must be multiplied by the spectrum for each pixel 
   /// and integrated over the energy bin to obtain the predicted counts
   std::vector<float> m_model;

   /// This is the "sparse" version of the source map data.
   SparseVector<float> m_sparseModel;

   /// What type of source map data do we have
   FileUtils::SrcMapType m_mapType;

   /// These are the 'spectrum' values
   /// I.e., the spectrum evaluated at the energy points
   std::vector<double> m_specVals;

   /// These are the spectral parameters for this source.
   /// They are used to determine if the spectrum has changed
   std::vector<double> m_modelPars;

   /// These are the 'latched' spectral parameters for this source.
   /// I.e., the version that the binned likelihood 
   /// Thinks the source has
   std::vector<double> m_latchedModelPars;

   /// These are the derivatives of the 'spectrum' values
   /// I.e., the derivatives evaluated at the energy points
   std::vector<std::vector<double> > m_derivs;

   /// These are the npreds, i.e., the model summed over each energy plane.
   std::vector<double> m_npreds;

   /// These are the npred weights, i.e., the weights to apply to the npreds
   /// to correctly reproduce the weighted counts
   std::vector<std::pair<double,double> > m_npred_weights;

   /// Caches of the true and measured energy spectra for sources
   Drm_Cache* m_drm_cache;


};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
