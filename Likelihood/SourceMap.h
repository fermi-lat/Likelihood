/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.64 2016/09/13 19:26:21 echarles Exp $
 */

#ifndef Likelihood_SourceMap_h
#define Likelihood_SourceMap_h

#include "st_facilities/libStApiExports.h"

#include "Likelihood/BinnedConfig.h"
#include "Likelihood/BinnedExposure.h"
#include "Likelihood/MeanPsf.h"
#include "Likelihood/Pixel.h"

#include "healpix_map.h"

namespace astro {
   class HealpixProj;
}

namespace st_stream {
   class StreamFormatter;
}

namespace Likelihood {
  
   // EAC, switch to using CountsMapBase and projection specific sub-classes
   class CountsMapBase;
   class CountsMap;
   class CountsMapHealpix;
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

   /* Standar c'tor 
      
      src           : The source making this source map for
      dataMap       : The counts map used as a template for the binning
      observation   : Object with data about the observation
      psf_config    : Object with PSF integration parameters
      drm           : The detector response matrix, NULL if energy dispersion is off for this source.
      weights       : A weights map, for weighted likelihood analysis.  Null -> no weighting.
    */
   SourceMap(const Source& src, 
	     const CountsMapBase * dataMap,
             const Observation & observation,
	     const PsfIntegConfig & psf_config, 
	     const Drm* drm = 0,
	     const WeightMap* weights = 0);


   /* C'tor to re-read this from a file
      
      sourceMapsFile : The path to the file we are reading from
      src           : The source in question
      dataMap       : The counts map used as a template for the binning
      observation   : Object with data about the observation
      weights       : A weights map, for weighted likelihood analysis.  Null -> no weighting.
      drm           : The detector response matrix, NULL if energy dispersion is off for this source.
    */
   SourceMap(const std::string & sourceMapsFile,
             const Source & src,
	     const CountsMapBase * dataMap,
	     const Observation & observation,
	     const WeightMap* weights = 0,
	     const Drm* drm = 0);

   /* d'tor, does clean up */
   ~SourceMap();

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


   /* --------------- Class Methods ----------------------*/

   /* Clear out the model, but save the cached data.
      This is useful for dealing with fixed sources, as it frees
      up a lot of memory.       
   */
   void clear_model() {
     m_model.clear();
   }      

   /* Extract a vector of spectral normalization values from a Source object
      and latch it in this class.
      
      energies:  The energies at which to evalute the spectrum
   */   
   void setSpectralValues(const std::vector<double>& energies);

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

   /// Explicitly set the image data
   void setImage(const std::vector<float>& model);

   /// Explicitly set the weights data
   void setWeights(const WeightMap* weights);

   /// Set the filename (e.g., b/c we are writing the source map)
   inline void setFilename(const std::string& filename) { m_filename = filename; }

protected:

   /* Read the model from a file */
   void readModel(const std::string& sourceMapFile);

   /* Read an image from a FITS file */
   void readImage(const std::string& sourceMapFile);

   /* Read HEALPix data from a table in a FITS file */
   void readTable_healpix(const std::string& sourceMapFile);

   /* Make the model */
   int make_model();

private:

   /* Compute the NPreds */
   void computeNpredArray();

   /* Fill this by projecting a weight map into the counts map frame 
      
      If any of the counts map pixels or energies are outside the weight_map
      they will be set to one and extrapolated will be set to true.
    */
   void makeProjectedMap(const ProjMap& weight_map, bool& extrapolated);

   /* Adjust the source map for a phased exposure correction       
      This is just multiplying the source map by the phased exposure map projected into the counts map frame
   */
   void applyPhasedExposureMap();



   /* ---------------- Data Members --------------------- */

   /// The source in question. This will by NULL if the SourceMap is the weights map
   const Source* m_src;   

   /// Source name.
   const std::string m_name;

   /// The name of the file where the source map is stored
   std::string m_filename;

   /// This is either "Diffuse" or "Point"
   const std::string m_srcType;

   /// Pointer to the data map that is the template for the binning
   const CountsMapBase * m_dataMap;
   
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

   /// This is the 'source map' data. 
   /// It is not the counts map, but rather the coefficients
   /// That must be multiplied by the spectrum for each pixel 
   /// and integrated over the energy bin to obtain the predicted counts
   std::vector<float> m_model;

   /// These are the 'spectrum' values
   /// I.e., the spectrum evaluated at the energy points
   std::vector<double> m_specVals;

   /// These are the spectral parameters for this source.
   /// They are used to determine if the spectrum has changed
   std::vector<double> m_modelPars;

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

   /// @brief Scaling factor between the true and projected angular separation
   //  for the fast PSF integration
   std::vector< std::vector< double > > m_pixelOffset;

};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
