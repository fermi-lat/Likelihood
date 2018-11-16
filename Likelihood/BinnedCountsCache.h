/**
 * @file BinnedCountsCache.h
 * @brief Small class to encapsulate stuff that needed for BinnedLikelihood.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/BinnedCountsCache.h,v 1.3 2017/10/06 01:30:59 echarles Exp $
 */

#ifndef Likelihood_BinnedCountsCache_h
#define Likelihood_BinnedCountsCache_h

#include <string>
#include <vector>

#include "Likelihood/CountsMapBase.h"



namespace tip {
   class Extension;
}


namespace Likelihood {

  class WeightMap;
  class ProjMap;
  class Observation;

  /*
   * @class BinnedCountsCache
   * @brief Cache of CountsMap, weights maps and other things we pass around for BinnedLikelihood
   *
   */
  
  class BinnedCountsCache {
    
  public:
    
    /* Regular c'tor 
       
       dataMap      : Observed data.  Defines the binning used for the analysis.
       observation  : Wrapper containing information about the observation
       weightMap    : Map with weights to use for analysis, NULL for no-weights
       srcMapsFile  : Name of file containing Source Maps (i.e., the output of gtsrcmaps)
       overwriteWeights : Overwrite the existing weights file
    */	
    BinnedCountsCache(CountsMapBase & dataMap, 
		      const Observation & observation,
		      const ProjMap* weightMap,
		      const std::string & srcMapsFile,
		      bool overwriteWeights = false);
    
    ///
    virtual ~BinnedCountsCache();
    
    /* --------------- Simple Access functions ----------------------*/
    
    /// Return the binned observed data
    inline const CountsMapBase & countsMap() const { return m_dataMap; }
    
    /// Return the energy bin edges
    inline const std::vector<double> & energies() const { return m_dataMap.energies(); }

    /// Return the log of the ratios of the energy bin edges
    inline const std::vector<double> & log_energy_ratios() const { return m_log_energy_ratios; }
    
    /// Return the number of energies
    inline size_t num_energies() const { return m_dataMap.num_energies(); }
    
    /// Return the number of energy bins
    inline size_t num_ebins() const { return m_dataMap.num_ebins(); }

    /// Return the number of pixels
    inline size_t num_pixels() const { return m_numPixels; }
    
    /// Return the size of the data map
    inline size_t data_map_size() const { return m_dataMap.data().size(); }
    
    /// Return the size of the scours maps
    inline size_t source_map_size() const { return m_dataMap.num_energies() * m_numPixels; }
    
    /// Return the number of filled pixels
    inline size_t nFilled() const { return m_filledPixels.size(); }

    /// Return the weights in their original projection
    inline const ProjMap* weightMap_orig() const { return m_weightMap_orig; }
    
    /// Return the weights reprojected into counts map binning
    inline const WeightMap* weightMap() const { return m_weightMap; }
    
    /// Return the weighted counts map
    inline const CountsMapBase* weightedCounts() const { return m_weightedCounts; }
    
    /// Return true if the weighted counts map exists
    inline bool has_weights() const { return m_weightedCounts != 0; }

    /// Get the data
    inline const std::vector<float>& data(bool weighted = false) const { 
      return weighted ? m_weightedCounts->data() : m_dataMap.data();
    }

    /// Return the observed counts spectrum
    inline const std::vector<double> & countsSpectrum(bool weighted=false) const { 
      return weighted ? m_countsSpectrum_wt : m_countsSpectrum; 
    }
    
    /// Get the filled pixel indices
    inline const std::vector<unsigned int>& filledPixels() const { return m_filledPixels; }

    /// Get the first pixels in each energy layer
    inline const std::vector<size_t>& firstPixels() const { return m_firstPixels; }

    /// Set the counts map by hand
    void setCountsMap(const std::vector<float> & counts);

    /// Set the weights map 
    void setWeightsMap(const ProjMap* wmap, 
		       const Observation & observation);


    /* Save the weights SourceMap to the SourceMap file
       
       replace : if true, replace the current version 
    */
    tip::Extension* saveWeightsMap(const std::string& srcMapsFile, 
				   bool replace=false) const;
    
    /// Fill the map of weighted counts
    void fillWeightedCounts();
 
    /* Fills the m_filledPixels data member with only the pixels with data counts */
    void identifyFilledPixels();
 

    void computeCountsSpectrum();
    
    void computeCountsSpectrum_wcs();
    
    void computeCountsSpectrum_healpix();

 
  protected:
    
    /// Disable assignement operator
    BinnedCountsCache & operator=(const BinnedCountsCache & rhs) {
      throw std::runtime_error("Copy-assignment operator of BinnedCountsCache not implemented");
    }
    
    /// Disable clone function
    virtual BinnedCountsCache * clone() const {
      return new BinnedCountsCache(*this);
    }
    
    
  private:
    
    /// Compute the log of rations between energy bin edges
    static void log_energy_ratios(const std::vector<double>& energies,
				  std::vector<double>& log_ratios);
    
    /* ---------------- Data Members --------------------- */
    
    /* ---------------- Data and binning --------------------- */
    
    /// The observed data.  Used to provide the binning.
    CountsMapBase& m_dataMap;

    /* ---------------- Counts Spectra ----------------------- */
    
    /// The observed counts spectrum.  
    /// I.e., the data summed of the ROI for each energy bin
    std::vector<double> m_countsSpectrum;
    
    /// The weighted observed counts spectrum.  
    /// I.e., the data summed of the ROI for each energy bin
    std::vector<double> m_countsSpectrum_wt;
    
    /// The number of pixels in the counts map
    size_t m_numPixels; 
    
    /// Log of ratios between energy bin edges
    std::vector<double> m_log_energy_ratios;
    
    /* ---------------- Stuff for weighted likelihood --------------- */
    
    /// Weights map.  Null ptr -> don't use weights.
    const ProjMap* m_weightMap_orig;
    
    /// Weights map reprojected into counts map binning.  Null ptr -> don't use weights 
    WeightMap* m_weightMap;
    
    /// Map of the weighted counts
    CountsMapBase* m_weightedCounts;

    /// These are the indices of the pixels with counts     
    /// This is used to speed up the evaluation of the log-likelihood
    std::vector<unsigned int> m_filledPixels;

    /// These are the indices of the first pixel in each energy layer m_filledPixels vector
    /// This is used to speed up the evaluation of the log-likelihood
    std::vector<size_t> m_firstPixels;
   
  };

}

#endif // Likelihood_BinnedCountsCache_h
