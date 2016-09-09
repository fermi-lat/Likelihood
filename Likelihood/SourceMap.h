/**
 * @file SourceMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/SourceMap.h,v 1.62 2016/04/16 00:50:25 echarles Exp $
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
   class DiffuseSource;
   class MeanPsf;
   class PointSource;
   class Source;
   class WcsMap2;

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
      weights       : A weights map, for weighted likelihood analysis.  Null -> no weighting.
    */
   SourceMap(Source * src, const CountsMapBase * dataMap,
             const Observation & observation,
	     const PsfIntegConfig & psf_config, 
	     const SourceMap* weights = 0);

   /* C'tor used to build a weights map

      This will reproject the weights map into the binning used for the analysis

      weight_map    : The weights map (in an arbitraty binning and projection )
      dataMap       : The counts map used as a template for the binning
      observation   : Object with data about the observation
      verbose       : Print out progress messages
    */
   SourceMap(const ProjMap& weight_map,
	     const CountsMapBase * dataMap,
             const Observation & observation,
	     bool verbose=true);


   /* C'tor to re-read this from a file
      
      sourceMapsFile : The path to the file we are reading from
      srcName       : Name of the source, also name of the HDU in the file
      observation   : Object with data about the observation
      weights       : A weights map, for weighted likelihood analysis.  Null -> no weighting.
      isWeights     : If true this is a weights source map
    */
   SourceMap(const std::string & sourceMapsFile,
             const std::string & srcName,
             const Observation & observation,
	     const SourceMap* weights = 0,
	     bool isWeights = false);

   /* d'tor, does clean up */
   ~SourceMap();

   /* --------------- Simple Access functions ----------------------*/

   /* The name of the source */
   inline const std::string & name() const { return m_name; }

   /* The source type, either 'Point' or 'Diffuse' */
   inline const std::string & srcType() const { return m_srcType; }

   /* The parameters for the PSF Integration */
   inline const PsfIntegConfig& psf_config() const { return m_psf_config; }

   /* The weights for the weighted log-likelihood.  Null-> no weights */
   inline const SourceMap* weights() const { return m_weights; } 


   /* --------------- Class Methods ----------------------*/

   /// These functions return cached data.  
   /// This might not be calculated yet, or might be out of date.
   /// Use with caution.
   /// The version of the functions without cached_ will always return 
   /// up-to-date information.
   inline const std::vector<float> & cached_model() const { return m_model; }

   inline const std::vector<double> & cached_npreds() const { return m_npreds; }

   inline const std::vector<std::pair<double,double> > & cached_npred_weights() const { return m_npred_weights; }
   
   /// These functions force the calculation of the quanitities they return;  
   /// FIXME, implement this
   const std::vector<float> & model() const;

   const std::vector<double> & npreds() const;

   const std::vector<std::pair<double,double> > & npred_weights() const;
   
   /// Explicitly set the image data
   void setImage(const std::vector<float>& model);

   /// Explicitly set the weights data
   void setWeights(const SourceMap* weights);


protected:

   /* Read an image from a FITS file */
   void readImage(const std::string& sourceMapFile);

   /* Read HEALPix data from a table in a FITS file */
   void readTable_healpix(const std::string& sourceMapFile);


private:

   /* Compute the NPreds */
   void computeNpredArray(bool isWeights=false);

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

   /// Source name.
   const std::string m_name;

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

   /// A weights map, for weighted likelihood analysis.  Null -> no weighting.   
   const SourceMap* m_weights;

   /// Flag to indicate that the data map is owned and should be deleted in d'tor
   bool m_deleteDataMap;

   /// This is the 'source map' data. 
   /// It is not the counts map, but rather the coefficients
   /// That must be multiplied by the spectrum for each pixel 
   /// and integrated over the energy bin to obtain the predicted counts
   mutable std::vector<float> m_model;

   /// These are the npreds, i.e., the model summed over each energy plane.
   mutable std::vector<double> m_npreds;

   /// These are the npred weights, FIXME (what are those?)
   mutable std::vector<std::pair<double,double> > m_npred_weights;

   /// These are for the fast PSF integration

   /// @brief Scaling factor between the true and projected angular separation
   std::vector< std::vector< double > > m_pixelOffset;

};

} // namespace Likelihood

#endif // Likelihood_SourceMap_h
