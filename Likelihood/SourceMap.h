/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.79 2017/09/29 01:38:05 echarles Exp $
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
      config        : Object with PSF integration parameters
      drm           : The detector response matrix, NULL if energy dispersion is off for this source.
      weights       : A weights map, for weighted likelihood analysis.  Null -> no weighting.
      save_model    : Flag to indicate that we should save the model (e.g., when the source is fixed)
                      This avoid reload the model map if the source is subsequently freed.
    */
   SourceMap(const Source & src, 
	     const BinnedCountsCache * dataCache,
             const Observation & observation,
	     const BinnedLikeConfig & config, 
	     const Drm& drm,
	     const WeightMap* weights = 0,
	     bool save_model = true);


   /* C'tor to re-read this from a file
      
      sourceMapsFile : The path to the file we are reading from
      src           : The source in question
      dataMap       : The counts map used as a template for the binning
      observation   : Object with data about the observation
      weights       : A weights map, for weighted likelihood analysis.  Null -> no weighting.
      drm           : The detector response matrix, NULL if energy dispersion is off for this source.
      save_model    : Flag to indicate that we should save the model (e.g., when the source is fixed)
                      This avoids reload the model map if the source is subsequently freed.
   */
   SourceMap(const std::string & sourceMapsFile,
             const Source & src,
	     const BinnedCountsCache * dataCache,
	     const Observation & observation,
	     const BinnedLikeConfig& config,
	     const Drm& drm,
	     const WeightMap* weights = 0,
	     bool save_model = true);

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
   inline const BinnedLikeConfig& config() const { return m_config; }

   /* The detector response matrix */
   inline const Drm* drm() const { return m_drm; }

   /* How to apply the energy dispersion */
   inline int edisp_val() const { return m_edisp_val; }

   /* The number of extra energy bins in the SourceMap, as compared to the CountsMap */
   inline size_t edisp_bins() const { return m_edisp_bins; }

   /* The offset between the extra energy bins in the DRM and in this source map */
   inline int edisp_offset() const { return m_edisp_offset; }

   /* The weights for the weighted log-likelihood.  Null-> no weights */
   inline const WeightMap* weights() const { return m_weights; } 

   /* How the source map is stored */
   inline FileUtils::SrcMapType mapType() const { return m_mapType; }

   /* Flag to indicate that we should save the model */
   inline bool save_model() const { return m_save_model; }

   /* Flag to indicated that model is out of sync with file */
   inline bool model_is_local() const { return m_model_is_local; }


   /* --------------- Class Methods ----------------------*/

   /* Clear out the model, but save the cached data.
      This is useful for dealing with fixed sources, as it frees
      up a lot of memory.       
   */
   void clear_model(bool force=false) {
     if ( force || ( !m_save_model && !m_model_is_local ) ) {
       // std::cout << "clear_model() called on " << m_name << " " << std::boolalpha;
       // std::cout << force << " " << m_save_model << " " << m_model_is_local << std::endl;
       // This deallocates the memory used by the model
       // in C++-11 there is a function shrink_to_fit that we could use.
       std::vector<float> nullVect;
       m_model.swap(nullVect);
       m_sparseModel.clear();
       m_model_is_local = false;
       m_dataCleared = true;
     }
     return;
   }      

   /* Set the source associated with this source map, this is useful for
      functions that add & remove source from the source model */
   void setSource(const Source& src);

   /* Update the DRM cache in this SourceMap */
   const Drm_Cache* update_drm_cache(const Drm* drm, bool force = false);

   /* The energies for this particular source map */
   inline const std::vector<double>& energies() const { return m_energies; }

   /* The log of the energy ratios for the particular source map */
   inline const std::vector<double>& log_energy_ratios() const { return m_logEnergyRatios; }

   /* The number of energies in the particular source map */
   inline size_t n_energies() const { return m_energies.size(); }

   /* The number of energy bins in this particular source map */
   inline size_t n_energy_bins() const { return m_logEnergyRatios.size(); }

   
   /* Extract a vector of spectral normalization values from a Source object
      and latch it in this class.
      
      latch_params : If true, the model parameters are latched
   */   
   void setSpectralValues(bool latch_params = false);

   /* Extract the spectral derivatives from the Source object and 
      latch them in this class.
      
      energies:  The energies at which to evalute the derivatives
      paraNames: The names of the params w.r.t. which to evaluate the derivaties
    */
   void setSpectralDerivs(const std::vector<std::string>& paramNames);
       

   /* Test to see if the spectrum has changes w.r.t. the cached values */
   bool spectrum_changed() const;


   /* Get a particular version of the counts specturm */
   const std::vector<double>& counts_spectra(int edisp_val,
					     bool use_weighted) const;


   /* Compute the total model counts summed between two energy bins *
      This uses the Drm_Cache, so that must be up to date.

      kmin         : index of lowest energy bin
      kman         : index of highest energy bin
      dataCache    :  Object with info about the binning
      edisp_val    : how to apply the energy dispersion
      use_weighted : return weighted model counts */
   double summed_counts(size_t kmin, size_t kmax,
			int edisp_val = 0,
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
   
   /* These are the 'spectrum' values.  I.e., the spectrum evaluated at the energy points */
   inline const std::vector<double> & cached_specValues() const { return m_specVals; }

   /* These are the 'spectrum' weights.  They are the quantity that appear in the log-log quadrature formula */
   inline const std::vector<std::pair<double,double> >& cached_specWts() const { return m_specWts; }

   /* These are the derivatives of the 'spectrum' values.  I.e., the derivatives evaluated 
      at the energy pointsThese are the spectral derivatives.  */
   inline const std::vector<std::vector<double> >& cached_specDerivs() const { return m_derivs; }

   /* These are the npreds, i.e., the model summed over each energy plane. */
   inline const std::vector<double>& cached_npreds() const { return  m_npreds; }

   /// These are the weighted npreds
   inline const std::vector<std::vector<std::pair<double,double> > >& cached_weighted_npreds() const { return m_weighted_npreds; }

   /* This contains both the convolved and un-convolved counts spectra */
   inline const Drm_Cache* cached_drm_Cache() const { return m_drm_cache; }

   /* The DRM used for the current cache, null -> no energy dispersion */
   inline const Drm* cached_drm() const { return m_drm; }

   /* This is the cached MeanPSF object */
   inline const MeanPsf* cached_meanPsf() const { return m_meanPsf; }

   /* This is the binning object */
   inline const BinnedCountsCache* dataCache() const { return m_dataCache; }
      


   /* --------------- Cached Data ----------------------*/

   /* These functions will either return the cached value or compute it if needed or if force = true */

   /* The source map model.  This must be multiplied by the spectrum for each pixel 
      and integrated over the energy bin to obtain the predicted counts */
   std::vector<float> & model(bool force=false);

   /* These are the 'spectrum' values, I.e., the spectrum evaluated at the energy points */
   const std::vector<double> & specVals(bool force=false, bool latch_params=false);

   /* These are the 'spectrum' weights.  They are the quantity that appear in the log-log quadrature formula */
   const std::vector<std::pair<double,double> >&  specWts(bool force=false);

   /* These are the derivatives of the 'spectrum' values.  I.e., the derivatives evaluated 
      at the energy pointsThese are the spectral derivatives.  */
   const std::vector<std::vector<double> >& specDerivs(const std::vector<std::string>& paramNames, bool force = false);
   
   /* These are the npreds, i.e., the model summed over each energy plane. */
   const std::vector<double> & npreds(bool force=false);

   /* These are weighted npreds */
   const std::vector<std::vector<std::pair<double,double> > >& weighted_npreds(bool force=false);
   
   /* This contains both the convolved and un-convolved counts spectra */
   const Drm_Cache* drm_cache(bool force=false);


   /* Add the model to an extermal vector */
   void addToVector(std::vector<float>& vect, bool includeSpec = false, int kmin=0, int kmax=-1);
   
   /* Subtract the model from an extermal vector */
   void subtractFromVector(std::vector<float>& vect, bool includeSpec = false, int kmin=0, int kmax=-1);

   /// Explicitly set the image data
   void setImage(const std::vector<float>& model);

   /// Explicitly set the weights data
   void setWeights(const WeightMap* weights);

   /// Explicity set the flag to save the model
   inline void setSaveModel(bool val) { 
      m_save_model = val; 
      if (val) reloadIfCleared();
//      std::cout << "setSaveModel(), m_save_model = " << std:: boolalpha << m_save_model << std::endl;
   }

   /// Set the filename (e.g., b/c we are writing the source map)
   inline void setFilename(const std::string& filename) { m_filename = filename; }

   /// Set value of model_is_local (e.g., b/c we are writing the source map)
   inline void setModelIsLocal(bool val) { m_model_is_local = val; }

   /// reload the modelfile if it has been cleared
   void reloadIfCleared();


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

   /* Latch the energy vector */
   void set_energies();

private:

   /* Implementation of addToVector for full-model storage */
   void addToVector_full(std::vector<float>& vect, bool includeSpec=false, int kmin=0, int kmax=-1) const;
 
   /* Implementation of addToVector for sparse-model storage */
   void addToVector_sparse(std::vector<float>& vect, bool includeSpec=false, int kmin=0, int kmax=-1) const;

   /* Implementation of subtractFromVector for full-model storage */
   void subtractFromVector_full(std::vector<float>& vect, bool includeSpec=false, int kmin=0, int kmax=-1) const;
 
   /* Implementation of subtractFromVector for sparse-model storage */
   void subtractFromVector_sparse(std::vector<float>& vect, bool includeSpec=false, int kmin=0, int kmax=-1) const;
 
   /* Compute the NPreds */
   void computeNpredArray();

   /* Compute the NPreds */
   void computeNpredArray_sparse();


   /* Adjust the source map for a phased exposure correction       
      This is just multiplying the source map by the phased exposure map projected into the counts map frame
   */
   void applyPhasedExposureMap();

   /* Adjust the source map for a phased exposure correction       
      This is just multiplying the source map by the phased exposure map projected into the counts map frame
      
      This version is used if we have a sparse map
   */   
   void applyPhasedExposureMap_sparse();


   /* Get the value from the sparse model */
   float find_value(size_t idx) const {
     return m_sparseModel[idx];
   }

   /* Reset data values after loading a source */
   void resetSourceData();

   /* Load or build source map data */
   void getSourceData();


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

   /// Options for treatment of PSF and energy disperson
   const BinnedLikeConfig& m_config;

   /// The detector response matrix. 
   const Drm* m_drm;

   /// How to apply the energy dispersion:
   //   m_edisp_val < 0 -> use rescaling
   //   m_edisp_val = 0 -> not applied
   //   m_edisp_val > 0 -> use complete method
   int m_edisp_val;

   /// The number of extra energy bins in the SourceMap, as compared to the CountsMap
   size_t m_edisp_bins;

   /// The offset between this source map energy binning and the Drm
   int m_edisp_offset;

   /// A weights map, for weighted likelihood analysis.  Null -> no weighting.   
   const WeightMap* m_weights;

   /// Flag to indicate that we should save the model
   bool m_save_model = true;

   /// Flag to indicated that model is out of sync with file
   bool m_model_is_local;

   /// The energies for this particular source map
   std::vector<double> m_energies;

   /// The log of the energy ratios
   std::vector<double> m_logEnergyRatios;   

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

   /// These are the 'spectrum' weights
   /// These are the quantities that appear in the log-log quadrature forumal
   std::vector<std::pair<double,double> > m_specWts;

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

   /// These are the weighted npreds
   std::vector<std::vector<std::pair<double,double> > > m_weighted_npreds;

   /// Caches of the true and measured energy spectra for sources
   Drm_Cache* m_drm_cache;

   /// flag indicating that data model has been cleared
   bool m_dataCleared = false;

   /// flag indicating that the data has been loaded at least once
   bool m_loaded = false;

};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
