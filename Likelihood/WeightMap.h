/**
 * @file WeightMap.h
 * @brief Spatial distribution of a source folded through the
 *        instrument response.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/WeightMap.h,v 1.64 2016/09/13 19:26:21 echarles Exp $
 */

#ifndef Likelihood_WeightMap_h
#define Likelihood_WeightMap_h

#include <vector>
#include <string>

#include "st_facilities/libStApiExports.h"

namespace st_stream {
   class StreamFormatter;
}

namespace Likelihood {
  
   // EAC, switch to using CountsMapBase and projection specific sub-classes
   class CountsMapBase;
   class Observation;
   class ProjMap;
   class Source;
   class WcsMap2;

   /*
    * @class WeightMap
    *
    */
   
#ifdef SWIG
   class WeightMap {
#else
   class  SCIENCETOOLS_API WeightMap {
#endif

   public:
     
     /* C'tor used to build a weights map
	
	This will reproject the weights map into the binning used for the analysis
	
	weight_map    : The weights map (in an arbitraty binning and projection )
	dataMap       : The counts map used as a template for the binning
	observation   : Object with data about the observation
	verbose       : Print out progress messages
     */
     WeightMap(const ProjMap& weight_map,
	       const CountsMapBase * dataMap,
	       const Observation & observation,
	       bool verbose=true);
     
     
     /* C'tor to re-read this from a file
	
	sourceMapsFile : The path to the file we are reading from
	dataMap       : The counts map used as a template for the binning
	observation   : Object with data about the observation
     */
     WeightMap(const std::string & sourceMapsFile,
	       const CountsMapBase * dataMap,
	       const Observation & observation);
     
     /* d'tor, does clean up */
     ~WeightMap();
     
     
     
     /* --------------- Cached Data ----------------------*/
     
     /* The source map model.  This must be multiplied by the spectrum for each pixel 
	and integrated over the energy bin to obtain the predicted counts */
     inline const std::vector<float> & model() const { return m_model; }
     
     
     /// Explicitly set the image data
     void setImage(const std::vector<float>& model);
         
   protected:
     
     /* Read an image from a FITS file */
     void readImage(const std::string& sourceMapFile);
     
     /* Read HEALPix data from a table in a FITS file */
     void readTable_healpix(const std::string& sourceMapFile);
     
   private:
     
     /* Fill this by projecting a weight map into the counts map frame 
	
	If any of the counts map pixels or energies are outside the weight_map
	they will be set to one and extrapolated will be set to true.
     */
     void makeProjectedMap(const ProjMap& weight_map, bool& extrapolated);
     
     
     /* ---------------- Data Members --------------------- */
     
     /// The name of the file where the weights map is store
     std::string m_srcMapFile; 

     /// Pointer to the data map that is the template for the binning
     const CountsMapBase * m_dataMap;
     
     /// Object holding data about the observation
     const Observation & m_observation;
     
     /// For Progress messages
     st_stream::StreamFormatter * m_formatter;
     
     /// This is the 'source map' data. 
     /// It is not the counts map, but rather the coefficients
     /// That must be multiplied by the spectrum for each pixel 
     /// and integrated over the energy bin to obtain the predicted counts
     std::vector<float> m_model;

   };

} // namespace Likelihood

#endif // Likelihood_WeightMap_h
