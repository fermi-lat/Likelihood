/**
 * @file FitScanner.h
 * @author Eric Charles
 *
 * A class to provide functionality that involves
 * scanning a test source over location or energy bins 
 * or normalization values.
 *
 * This is useful for making TSMaps, SEDs or likelihood v. flux maps 
 *
 * This class uses a number of approximations and shortcuts for speed.
 * 
 *   1) This only runs on binned data.
 *   2) Only the source normalizations are floated, not the other spectral 
 *  parameters.
 *   3) The fitting is done using Newton's method.  The normalizations 
 *  are not allowed to be negative, but this is simply enforced by not
 *  allowing the step between interations to be more negative than the 
 *  current value of the parameter.
 *   4) The gradiant (first derivatives) and Hessian (second derivatives) of the
 *  log-likelihood are calculated analytically.
 *   5) The convergence criterium is calculated using the estimated vertical 
 *  distance to the mininum.  This is computed as the inner product of the 
 *  current step sizes with the gradient vector, which is likely to be an 
 *  overestimate of the actual vertical distance to the minimum (by up to a 
 *  factor of 2)
 *   6) The expected counts model for the test source is computed at new locations
 *  by computing the counts model once at the center of the region of interest, 
 *  then translating that model by moving it a fixed number of pixels to each 
 *  new test location.   This will result in some projection inaccuracies.
 *
 *
 *  This file also contains two helper classes used by FitScanner:  
 *    TestSourceModelCache -> Used to cache the model image of the test source
 *    FitScanCache -> Used to cache the actual counts data and models for efficient fitting
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FitScanner.h,v 1.18 2016/07/02 01:07:07 echarles Exp $
 */

#ifndef Likelihood_FitScanner_h
#define Likelihood_FitScanner_h

// stl includes
#include <string>
#include <vector>
#include <map>

// CLHEP include

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

// Fermi includes
#include "astro/SkyDir.h"


// forward declarations
namespace astro {
  class SkyDir;
  //  class ProjBase;
  class SkyProj;
  class HealpixProj;
}

namespace evtbin {
  class Binner;
}

namespace tip {
  class Header;
}

namespace optimizers{
  class Optimizer;
  class FunctionFactory;
  class Function;
}


// This class lives in the Likelihood namespace
namespace Likelihood {

  /**
   * @class FitScanner
   * 
   */
   
  class BinnedLikelihood;
  class SummedLikelihood;
  class LogLike;
  class HistND;
  class Source;
  class PointSource;
  class AppHelpers;
  class CountsMapBase;  
  struct Snapshot_Status;
  class Snapshot;

  /* A utility class to cache the image the predicted counts map for the
     test source.  
     
     Simply moving the image around is much faster that recomputing it for each 
     point in the TS-map grid. 

     For small grids this is not to bad of an approximation.

     FIXME, add details on size of errors expected.
   */
  class TestSourceModelCache {

  public:

    /* Build from a BinnedLikelihood and a Source */
    TestSourceModelCache(const BinnedLikelihood& logLike,
			 const PointSource& source);

    /* D'tor, does nothing */
    ~TestSourceModelCache(){};
    
    /* Translate the cached map to a new location 

       newRef      : The new direction of the center of the model image
       out_model   : Filled with the values of the new model image

       returns 0 for success, error code for failure
     */
    int translateMap(const astro::SkyDir& newRef,
		     std::vector<float>& out_model) const;

    /* Write the current cached map to a FITS image
       
       fits_file    : Name of the fits file in question
       ext_name     : Nmae of the iamge extension
    */
    void writeTestSourceToFitsImage(const std::string& fits_file,
				    const std::string& ext_name) const; 
    

  protected:

    /* Translate the map using WCS projection by dx and dy pixels 

       dx          : Number of pixels offset in X
       dy          : Number of pixels offset in Y
       out_model   : Filled with the values of the new model image     

       returns 0 for success, error code for failure
    */ 
    int translateMap_Wcs(double dx, double dy, std::vector<float>& out_model) const;

    /* Translate the map using HEALPix projection by d_theta and d_phi degrees */
    int translateMap_Healpix(double d_theta, double d_phi, std::vector<float>& out_model) const {
      // FIXME, implement this
      return -1;
    }

  private:
          
    /* The reference image */
    std::vector<float> m_refModel;
    /* The projection used for both the reference and scanned images*/
    const astro::ProjBase& m_proj;
    /* The reference direction */
    const astro::SkyDir& m_refDir;
    

    /* The location of the reference pixel in image coordinates */
    std::pair<double,double> m_refPixel;

    /* Number of bins in x,y and energy */
    size_t m_nx;
    size_t m_ny;
    size_t m_ne;

    /* The current image */
    mutable std::vector<float> m_currentModel;

  };


  /* A utility class to store and apply a multivariate prior on fit parameters.

     The prior is expressed as a set of best fit values, and the covariance matrix between them.
     This corresponds to a set of Gaussian-distributed errors about the best-fit estimates.     
     
     In partical terms, the prior is applied by adding an addition term to the 
     loglikelihood, as therefore also to the gradiant and Hessian during the fitting.
   */
  class FitScanMVPrior {

  public:

    /* Construct from a vector of central values and a covariance matrix 
       
       This version reduces the matrix by cutting out all of the 
       unconstrainted paramters.
     */
    FitScanMVPrior(const CLHEP::HepVector& centralVals,
		   const CLHEP::HepSymMatrix& covariance,
		   const std::vector<bool>& constrainPars,
		   bool includeTestSource){
      update(centralVals,covariance,constrainPars,includeTestSource);
    }
    
    /* Construct from a vector of central values and a vector of errors 
       
       This version does not reduce the matrix, but rather sets 
       the non-constrained terms row & columns to zero in the hessian

    */
    FitScanMVPrior(const CLHEP::HepVector& centralVals,
		   const CLHEP::HepVector& uncertainties,
		   const std::vector<bool>& constrainPars,
		   bool includeTestSource){
      update(centralVals,uncertainties,constrainPars,includeTestSource);
    }
   
    /* D'tor, does nothing */
    ~FitScanMVPrior() {;}

    /* Update the cached values in this prior

       centralVals       : Central values from fit we are using to build this prior
       covariance        : Covariance matrix from fit we are using to build this prior
       constrainPars     : Flags showing which parameters to include in this prior     
       includeTestSource : If true, expand the dimension by one to allow for the test source
       
       This version reduces the matrix by cutting out all of the 
       unconstrainted paramters.
     */
    void update(const CLHEP::HepVector& centralVals,
		const CLHEP::HepSymMatrix& covariance,
		const std::vector<bool>& constrainPars,
		bool includeTestSource);

    
    /* Update the cached values in this prior

       centralVals       : Central values we are using to build this prior
       uncertainties     : Uncertainties we are we are using to build this prior
       constrainPars     : Flags showing which parameters to include in this prior     
       includeTestSource : If true, expand the dimension by one to allow for the test source

       This version does not reduce the matrix, but rather sets 
       the non-constrained terms row & columns to zero in the hessian
    */
    void update(const CLHEP::HepVector& centralVals,
		const CLHEP::HepVector& uncertainties,
		const std::vector<bool>& constrainPars,
		bool includeTestSource);


    // Calculate the contribution to the log-likelihood
    void logLikelihood(const CLHEP::HepVector& params, 
		       double& logLike) const;

    // Calculate the contribution to the gradient of the log-likelihood
    void gradient(const CLHEP::HepVector& params, 
		  CLHEP::HepVector& grad) const;

    // Calculate the contribution to the gradient 2nd derivativs of the log-likelihood.
    // By construction this is actually just the constant Hessian matrix
    inline const CLHEP::HepSymMatrix& hessian() const { 
      return m_hessian; 
    }

    // These are the central values that we are constraining towards
    inline const CLHEP::HepVector& centralVals() const { return m_centralVals; }    
    // This covariance matrix that gives the strength of the constraints
    inline const CLHEP::HepSymMatrix& covariance() const { return m_covariance; }
    /* This vector is the size of _all_ of the cached parameters,
       it tells use which should be constrained (all the other are FIXED).
       To have a free parameter just put in a very weak constraint */
    inline const std::vector<bool> constrainPars() const { return m_constrainPars; }
    /* Is the test source include in this constraint */
    inline bool includeTestSource() const { return m_includeTestSource; }

  protected:
    
    int latchReducedMatrix();
    
  private:
    
    CLHEP::HepVector m_centralVals;
    CLHEP::HepSymMatrix m_covariance;
    std::vector<bool> m_constrainPars;
    bool m_includeTestSource;

    CLHEP::HepSymMatrix m_hessian;
  };

  

   /* A utility class to wrap either a BinnedLikelihood or a SummedLikelihood
      and provide a common interface to the two classes
   */
  
  class FitScanModelWrapper {

  public:
    
    /* C'tor, trivial */
    FitScanModelWrapper()
      :m_npix(0),m_nebins(0),m_size(0){;}

    /* D'tor, trivial */
    ~FitScanModelWrapper() {;}    

    /* return a reference to the counts data   
       for the BinnedLikelihood this point to the CountMap in the BinnedLikelihood object
       for the SummedLikelihood this points to a local vector where we have merged the countsdata
    */
    virtual const std::vector<float>& data() const = 0;

    // Is this a SummedLikelihood or a BinnedLikelihood
    virtual bool isSummed() const = 0;

    /* Extract the predicted counts model given the name of a source
       for the BinnedLikelihood this just calls FitUtils::extractModelFromSource
       for the SummedLikelihood this merged together the models for the various components
    */       
    virtual void extractModelFromSource(const std::string& srcName,
					std::vector<float>& model,
					bool rescaleToNormOne = false) const;
    

    /* Extract the predicted counts model from a Source object
       for the BinnedLikelihood this just calls FitUtils::extractModelFromSource
       for the SummedLikelihood this merged together the models for the various components
    */       
    virtual void extractModelFromSource(Source& aSrc,
					std::vector<float>& model,
					bool rescaleToNormOne = false) const = 0;

    /* Extract the fitting templates 
       for the BinnedLikelihood this just calls FitUtils::extractModels
       for the SummedLikelihood this merged together the models for the various components
     */
    virtual void extractModels(const std::string& test_name,
			       std::vector<std::string>& freeSrcNames,
			       std::vector<std::vector<float> >& templates,		       
			       std::vector<float>& fixed,
			       std::vector<float>& test_source_model,
			       std::vector<float>& refPars,
			       std::vector<float>* weights = 0) const = 0;

    /* Extract the priors on the free source normalizations from the 
       underlying BinnedLikelihood object

       freeSrcNames      : Names of the sources 
       centralVals       : Filled with the central values from the priors
       uncertainties     : Filled with the untertainty values from the priors
       constrainPars     : Filled with the flags showing which parameters have priors

       return true if any priors were extract, false otherwise
    */ 
    virtual bool extractPriors(const std::vector<std::string>& freeSrcNames,
			       CLHEP::HepVector& centralVals,
			       CLHEP::HepVector& uncertainties,
			       std::vector<bool>& constrainPars) const;   

    
    
    /* Check to see if the fixed / free status of sources has changed
       w.r.t. the wrapped BinnedLikelihood or SummedLikelihood 
    */    
    bool fixed_changed(const std::vector<std::string>& srcNames,
		       const std::vector<bool>& currentFreeSources) const;

    /* Get the likelihood value */
    virtual double value() const = 0;

    /* Add a source to the model */
    virtual void addSource(Source* aSrc) = 0;

     /* Call syncParams on the model */
    virtual void syncParams() = 0;
  
    /* Remove a source from the model */
    virtual void removeSource(const std::string& sourceName) = 0;    

    /* Write the energy bins to a fits file */
    virtual int writeFits_EnergyBins(const std::string& fitsFile) const;

    /* write the Good time intervals */
    virtual int writeFits_GTIs(const std::string& fitsFile) const = 0;   

    /* get the "master" component
       for the BinnedLikelihood this is just the object itself
       for the SummedLikelihood this is just the first component 
    */
    virtual BinnedLikelihood& getMasterComponent() = 0;     

    /* get the "master" component
       for the BinnedLikelihood this is just the object itself
       for the SummedLikelihood this is just the first component 
    */
    virtual const BinnedLikelihood& getMasterComponent() const = 0;     

    /* get a component by index */
    virtual const size_t numComponents() const = 0;

    /* get a component by index */
    virtual const BinnedLikelihood* getComponent(size_t idx) const = 0;   

    /* get the energy bins */
    virtual const std::vector<double>& energies() const = 0;

    /* shift the test source */
    virtual int shiftTestSource(const std::vector<TestSourceModelCache*>& modelCaches,
				const astro::SkyDir& newDir,
				std::vector<float>& targetModel) const = 0;

    /* set the energy bins to use in the analysis */
    virtual void set_klims(size_t kmin, size_t kmax) = 0;

    /* cache the flux values at the energy bin edges */
    void cacheFluxValues(Source& aSrc);

    /* return the number of pixels,        
       for the SummedLikelihood this is summed over the components */
    inline size_t nPix() const { return m_npix; }

    /* return the number of energy bins, 
       for the SummedLikelihoods this is the union of the components  */    
    inline size_t nEBins() const { return m_nebins; }
    
    /* return the number of pixels*the number of energy bins */       
    inline size_t size() const { return m_size; }   

    /* return the reference energies (geometric mean of energy bins  */
    inline const std::vector<double>& ref_energies() const { return m_ref_energies; }

    /* return the differential fluxes for the test source model spectrum */
    inline const std::vector<double>& ref_dfdes() const { return m_ref_dfdes; }

    /* return the integral photon fluxes for the test source model spectrum */
    inline const std::vector<double>& ref_fluxes() const { return  m_ref_fluxes; }

    /* return the integral energy fluxes for the test source model spectrum */
    inline const std::vector<double> ref_energy_fluxes() const { return m_ref_energy_fluxes; }
  
    /* return the flux values at the energy bin edges */
    inline const std::vector<double>& nPreds() const { return m_nPreds; }

  protected:
    
    inline void setDims(size_t nPix, size_t nEBins) {
      m_npix = nPix; m_nebins = nEBins;
      m_size = m_npix*nEBins;
    }

    int writeFits_FluxTable(const std::string& fitsFile) const;

    /* set the TUNIT keyword */
    void setUnitKeyword(tip::Header& header,
			int icol,
			const std::string& unitString) const;

  private:
    
    size_t m_npix;
    size_t m_nebins;
    size_t m_size; 


    // Reference energies (geometric mean of energy bins)
    std::vector<double> m_ref_energies;    

    // differential fluxes for the test source model spectrum
    std::vector<double> m_ref_dfdes;

    // integral photon fluxes for the test source model spectrum
    std::vector<double> m_ref_fluxes;

    // integral energy fluxes for the test source model spectrum
    std::vector<double> m_ref_energy_fluxes;

    // nPreds for the test source model spectrum
    std::vector<double> m_nPreds;


  };

  
  class FitScanModelWrapper_Binned : public FitScanModelWrapper {

  public:

    /* C'tor from a BinnedLikelihood object */
    FitScanModelWrapper_Binned(BinnedLikelihood& binnedLike);

    /* D'tor, trivial */
    ~FitScanModelWrapper_Binned() {;}

    inline BinnedLikelihood& binnedLike() { return m_binnedLike; }

    /* return a reference to the counts data   
       for the BinnedLikelihood this point to the CountMap in the BinnedLikelihood object
    */
    virtual const std::vector<float>& data() const; 

    // Is this a SummedLikelihood or a BinnedLikelihood
    virtual bool isSummed() const { return false; }

    /* Extract the predicted counts model from a Source object
       for the BinnedLikelihood this just calls FitUtils::extractModelFromSource
    */       
    virtual void extractModelFromSource(Source& aSrc,
					std::vector<float>& model,
					bool rescaleToNormOne = false) const;

    /* Extract the fitting templates 
       for the BinnedLikelihood this just calls FitUtils::extractModels
     */
    virtual void extractModels(const std::string& test_name,
			       std::vector<std::string>& freeSrcNames,
			       std::vector<std::vector<float> >& templates,		       
			       std::vector<float>& fixed,
			       std::vector<float>& test_source_model,
			       std::vector<float>& refPars,
			       std::vector<float>* weights = 0) const;


    /* Get the likelihood value */
    virtual double value() const;

    /* Add a source to the model */
    virtual void addSource(Source* aSrc);

    /* Call syncParams on the model */
    virtual void syncParams();
    
    /* Remove a source from the model */
    virtual void removeSource(const std::string& sourceName);    

    /* write the Good time intervals */
    virtual int writeFits_GTIs(const std::string& fitsFile) const;   

    /* get the "master" component
       for the BinnedLikelihood this is just the object itself
    */
    virtual BinnedLikelihood& getMasterComponent() { return m_binnedLike; }
    virtual const BinnedLikelihood& getMasterComponent() const { return m_binnedLike; }

    /* get a component by index */
    virtual const size_t numComponents() const { return 1; }

    /* get a component by index */
    virtual const BinnedLikelihood* getComponent(size_t idx) const { return idx == 0 ? &m_binnedLike : 0; }
   
    /* get the energy bins */
    virtual const std::vector<double>& energies() const;

    /* shift the test source */
    virtual int shiftTestSource(const std::vector<TestSourceModelCache*>& modelCaches,
				const astro::SkyDir& newDir,
				std::vector<float>& targetModel) const;

    /* set the energy bins to use in the analysis */
    virtual void set_klims(size_t kmin, size_t kmax);

  private:

    BinnedLikelihood& m_binnedLike;

  };

  
  class FitScanModelWrapper_Summed : public FitScanModelWrapper {
     
  public:
    
    static double findMinAndMatches(const std::vector< const std::vector<double>* >& energyBins,
				    std::vector<size_t>& localIdx,
				    std::vector<int>& matches,
				    float tol = 1e-5);
    
    static void mergeVectors(const std::vector< const std::vector<float>* >& toMerge,
			     const std::vector<size_t>& npixelsByComp,
			     const std::vector< std::vector<int> >& energyBinLocal,
			     std::vector<float>& mergedData);
    
  public:
    
    /* C'tor from a SummedLikelihood object */
    FitScanModelWrapper_Summed(SummedLikelihood& summedLike);
    
    /* D'tor, trivial */
    ~FitScanModelWrapper_Summed() {;}
    
    inline SummedLikelihood& summedLike() { return m_summedLike; }    
     
    /* return a reference to the data,        
       for the BinnedLikelihood this point to the CountMap in the BinnedLikelihood object
       for the SummedLikelihood this points to a local vector where we have merged the countsdata
    */
    inline const std::vector<float>& data() const { return m_localData; }
    
    // Is this a SummedLikelihood or a BinnedLikelihood
    inline bool isSummed() const { return true; }
    
    /* Extract the predicted counts model from a Source object
       for the SummedLikelihood this merged together the models for the various components
    */       
    virtual void extractModelFromSource(Source& aSrc,
					std::vector<float>& model,
					bool rescaleToNormOne = false) const;
    
    /* Extract the fitting templates 
       for the SummedLikelihood this merged together the models for the various components
    */
    virtual void extractModels(const std::string& test_name,
			       std::vector<std::string>& freeSrcNames,
			       std::vector<std::vector<float> >& templates,		       
			       std::vector<float>& fixed,
			       std::vector<float>& test_source_model,
			       std::vector<float>& refPars,
			       std::vector<float>* weights = 0) const;
    
    /* Get the likelihood value */
    virtual double value() const;
    
    /* Add a source to the model */
    virtual void addSource(Source* aSrc);
    
    /* Call syncParams on the model */
    virtual void syncParams();
    
    /* Remove a source from the model */
    virtual void removeSource(const std::string& sourceName);    
    
    /* write the Good time intervals */
    virtual int writeFits_GTIs(const std::string& fitsFile) const;   
    
    /* get the "master" component
       for the SummedLikelihood this is just the first component 
    */
    virtual BinnedLikelihood& getMasterComponent() { return *m_master; }    
    virtual const BinnedLikelihood& getMasterComponent() const { return *m_master; }    
    
    /* get a component by index */
    virtual const size_t numComponents() const;
      
    /* get a component by index */
    virtual const BinnedLikelihood* getComponent(size_t idx) const;   

    /* get the energy bins */
    virtual const std::vector<double>& energies() const { return m_energiesMerged; }

    /* shift the test source */
    virtual int shiftTestSource(const std::vector<TestSourceModelCache*>& modelCaches,
				const astro::SkyDir& newDir,
				std::vector<float>& targetModel) const;
  
    /* set the energy bins to use in the analysis */
    virtual void set_klims(size_t kmin, size_t kmax);


  protected:
    
  private:
    
    SummedLikelihood& m_summedLike;
    
    BinnedLikelihood* m_master;
    
    // Merged version of the counts maps
    std::vector<float> m_localData;
    
    // Index of first pixels for each component
    std::vector<size_t> m_nPixelsByComp;

    // Number of enerby bins for each component
    std::vector<size_t> m_nEBinsByComp;

    // Size of each component
    std::vector<size_t> m_sizeByComp;    
    
    // Merged set of energy bins
    std::vector<double> m_energiesMerged;
    
    // True if a component has data for a particular energy bin
    std::vector< std::vector<int> > m_energyBinLocal;     
  };
  
  
  /* A utility class to extract the data need for fitting from a BinnedLikelihood object     
     and perform scans using that data, model and a test source.

     The extracted data are stored as series of std::vector<float>, where each vector
     corresponds to either a data map, or a predicted counts model for a specific model component.

     The fits are done with respect to the original source models extracted from the 
     BinnedLikelihood object.  The returned fit parameters are scale factors that need to 
     be applied to the normalization parameters in the XML model used by gtlike, with the
     exception of the test source, where the fit parameter is the normalization of the 
     test source.

     The FitScanCache::fitCurrent() function will fit the currently cached parameters and models.
     
     setTestSource() and shiftTestSource() can be used to move the location of the test source

     addTestSourceToCurrent() and removeTestSourceFromCurrent() can be used to control if the 
     test source is included in the fit

     refactorModel() can be used to fix and freeze other sources in the ROI.  Note that refactorModel()
     operates w.r.t. the baseline model of the region (i.e., the one that was copied in from the 
     BinnedLikelihood object)

     setEnergyBin() sets flags to only do the fit for a single energy bin

     scanNormalization() will do a series of fits, scanning the normalization parameter and refitting
     any other free parameters at each scan point (i.e., a "profile" likelihood)
     
  */  
  class FitScanCache {

  public:

    typedef enum { No_Prior, Local_Prior, Global_Prior, Init_Prior } Prior_Version; 
    typedef enum { No_Action=0,      // Model unchanged
		   Remake_Prior=0x1, // Need to remake the prior
		   Refactor=0x2,     // Model changed, but templates already latched
		   Update_Free=0x4,  // Free source changed, must relatch template
		   Update_Fixed=0x8, // Fixed source changed, must relatch fixed template
		   Rebuild=0x10       // Lots of stuff changed, rebuild the cache
    } Update_Action;

    static unsigned action_needed(Snapshot_Status stat);

  public:

    /* Build from a ModelWrapper object
       
       This copies over all the needed data from the ModelWrapper object,
       which wraps either a BinnedLikelihood or a SummedLikelihood
     */
    FitScanCache(FitScanModelWrapper& modelWrapper,
		 const std::string& testSourceName,
		 double tol, int maxIter, double initLambda,
		 bool useReduced, bool useWeights=false);

    /* D'tor */
    ~FitScanCache();

    /* Update stuff w.r.t. the wrapped BinnedLikelihood or SummedLikelihood 

       This first figures out the action needed, then calls update_with_action()
     */
    void update() {
      std::vector<std::string> changed_sources;
      unsigned action = find_action_needed(changed_sources);
      update_with_action(action,changed_sources);
    }

    /* Update stuff w.r.t. the wrapped BinnedLikelihood or SummedLikelihood 
       by taking a specific action.

       Don't use this unless you know what you are doing.
    */
    void update_with_action(unsigned action, const std::vector<std::string>& changed_sources);

    /* Decide the correct action needed to update this cache 
       w.r.t. the wrapped BinnedLikelihood or SummedLikelihood 
    */
    unsigned find_action_needed(std::vector<std::string>& changed_sources) const;    
    
    /* Reset everything to the initial master state */
    void setCache();

    /* Refactor the current model, fixing or freeing sources, and changing normalizations 
       
       freeSource     : Marks which source are free
       parScales      : Current values of the normalization parameters
       include_test   : If true, includes test source in fit model
     */
    void refactorModel(const std::vector<bool>& freeSources, 
		       const std::vector<float>& pars_scales,
		       bool include_test);

    /* Extract the scaling values for the normalization parameters for the current fit */
    void getParScales(std::vector<float>& pars_scales);
    
    /* Set the cache to only do the fit over a single energy bin */ 
    void setEnergyBin(int energyBin);

    /* Set the cache to only do the fit over a range of energy bins */ 
    void setEnergyBins(size_t firstEnergyBin, size_t lastEnergyBin);
    
    /* Update the model of the test source.
       This version recomputes the SourceMap and somewhat more expensive (and accurate)
       than shiftTestSource
    */
    void setTestSource(Source& aSrc);

    /* Update the model of the test source
       This version shift the SourceMap with respect to a precomputed version
       and in much less expensive, but not quite as accurate 
     */
    int shiftTestSource(const std::vector<TestSourceModelCache*>& modelCache,
			const astro::SkyDir& newDir);

    /* Set the cache to add in the test source with a specify normalization value */
    void addTestSourceToCurrent(double initNorm);

    /* Set the cache to remove the test source from the fit */
    void removeTestSourceFromCurrent();

    /* Set the prior */
    void buildPriorsFromExternal(const CLHEP::HepVector& centralVals,
				 const CLHEP::HepSymMatrix& covariance,
				 const std::vector<bool>& constrainPars,
				 bool globalPrior=false);

    /* Set the prior */
    void buildPriorsFromExternal(const std::vector<float>& centralVals,
				 const std::vector<float>& covariance,
				 const std::vector<bool>& constrainPars,
				 bool globalPrior=false);

    /* Set the prior from the current fit*/
    void buildPriorsFromCurrent(const std::vector<bool>& constrainPars,
				double covScaleFactor,
				bool globalPrior=false);
 
    /* Fit the currently cached values using Newton's method */
    int fitCurrent(Prior_Version whichPrior=No_Prior, int verbose=0);

    /* Calculate the log-likelihood for the currently cached values */
    int calculateLoglikeCurrent(double& logLike, Prior_Version whichPrior=No_Prior);

    /* Scan the log likelihood versus the normalizaton of the test source 
       
       nnorm      : Number of scan points
       normSigma  : Number of sigma +/- to scan over
       pos_errs   : Positive side errors on the scan
       neg_errs   : Negative side errors on the scan
       norms      : Filled with normalization values used in the scan
       logLikes   : Filled with the corresponding log likelihood values     
     */
    int scanNormalization(int nnorm, double normSigma,
			  double posErr, double negErr,
			  std::vector<double>& norms,
			  std::vector<double>& logLikes);

    /* Estimate the uncertainty using quadratic equation.

       This accounts for both the normal quadratic dependence,
       and also the linear dependence that you get when 
       the signal is near zero.
       
       deltaLogLike  : log-likelihood value we are solving for
       pos_errs      : Positive side errors on the scan
       neg_errs      : Negative side errors on the scan
    */
    int signalUncertainty_quad(double deltaLogLike,
			       double& posErr,
			       double& negErr);

    
    /* get the index (int the master list) corresponding to a source by name */
    int getTemplateIndex(const std::string& srcName) const;

    /* update the model template for a single source by name */
    int updateTemplateForSource(const std::string& srcName);


    // access --------------------------------------------------------

    // Information about the baseline model
    inline const std::string& testSourceName() const { return m_testSourceName; }
    inline size_t npix() const { return m_npix; }
    inline size_t nebins() const { return m_nebins; }
    inline size_t nBkgModel() const { return m_allModels.size(); }
    inline bool useReduced() const { return m_useReduced; }
    inline bool useWeights() const { return m_useWeights; }
    inline double tol() const { return m_tol; }
    inline int maxIter() const { return m_maxIter; }
    inline double initLambda() const { return m_initLambda; }

    inline const double& loglike_ref() const { return m_loglike_ref; }    

    // Information about the currently cached model
    inline size_t nFreeCurrent() const { return m_currentModels.size(); }
    inline int testSourceIndex() const { return m_currentTestSourceIndex; }
    inline const std::vector<bool>& currentFree() const { return m_currentFreeSources; }

    // Information about the current fit
    inline const CLHEP::HepVector& currentPars() const { return m_currentPars; }
    inline const CLHEP::HepSymMatrix& currentCov() const { return m_currentCov; }
    inline const std::vector<float>& currentModel() const { return m_currentBestModel; }
    inline double currentLogLike() const { return m_currentLogLike; }
    inline double currentEDM() const { return m_currentEDM; }    
    inline int firstEnergyBin() const { return m_firstEnergyBin; }
    inline int lastEnergyBin() const { return m_lastEnergyBin; }

    // parts of the models
    inline const std::vector<float>& refValues() const { return m_refValues; }
    inline const std::vector<std::string>& templateSourceNames() const { return m_templateSourceNames; }
    inline const std::vector<std::vector<float> >& allModels() const { return m_allModels; }
    inline const std::vector<float>& allFixed() const { return m_allFixed; }   
    inline const std::vector<float>& weights() const { return m_weights; }   
    inline const std::vector<float>& targetModel() const { return m_targetModel; }

    // info about the iteration
    inline size_t firstBin() const { return m_firstBin; }
    inline size_t lastBin() const { return m_lastBin; }
    
    inline const std::vector<float>& currentFixed() const { return m_currentFixed; }
    inline const std::vector<float>& targetRedModel() const { return m_targetRedModel; }    
    inline const std::vector<float>& weightsRed() const { return m_weightsRed; }   

    // access to the priors
    const FitScanMVPrior* getPrior(Prior_Version whichPrior=No_Prior, 
				   bool include_test_source=false) const;

    // access to the snapshot of the baseline model when it was built
    inline const Snapshot* snapshot() const { return m_snapshot; }

    // access to the SourceModel
    inline const BinnedLikelihood& sourceModel() const { return m_modelWrapper.getMasterComponent(); }

  protected:

    void reduceModels();

    void build_from_model();
    
    void update_fixed_from_model();

    void update_free_from_model(const std::vector<std::string>& changed_sources);

    void refactor_from_model();

    void get_status_and_scales(std::vector<bool>& freeSources,
			       std::vector<float>& pars_scales);

    void extract_init_priors_from_model();

    void cleanup();


  private:

    // The wrapper around the model object
    FitScanModelWrapper& m_modelWrapper;

    // A snapshot of the reference model
    Snapshot* m_snapshot;

    // The name of the test source
    const std::string m_testSourceName;
    
    // Fit tolerance
    double m_tol;    
    // maximum number of interactions
    int m_maxIter;
    // Initial damping parameter for step calculation
    double m_initLambda;
       
    // number of energy bins in the counts map
    size_t m_nebins;
    // number of pixels in the counts map
    size_t m_npix;

    // reference to counts data
    const std::vector<float>& m_data;
    // master list of all the models, except for the test source
    std::vector<std::vector<float> > m_allModels;
    // master list of the names of the sources
    std::vector<std::string> m_templateSourceNames;
    // master version of the sum of the fixed model components
    std::vector<float> m_allFixed;
    // master version of the initial model
    std::vector<float> m_initModel;
    // master version of the weights
    std::vector<float> m_weights;
    // set of refrence values of all the paramters
    std::vector<float> m_refValues;
    // current test source model map
    std::vector<float> m_targetModel;    

    // specific for sparse model fitting
					       
    // Use the reduced vectors?
    bool m_useReduced;
    // Use the weights
    bool m_useWeights;

    // Reduced data vector
    std::vector<float> m_dataRed;
    // Indices of non-zero bins
    std::vector<int> m_nonZeroBins;
    // Indices of last bin in each energy bin
    std::vector<int> m_energyBinStopIdxs;
    // master list of all the reduced models, except for the test source  
    std::vector<std::vector<float> > m_allRedModels;
    // master version of the reduced models sum of the fixed model components
    std::vector<float> m_allRedFixed;
    // master version of the reduced weights
    std::vector<float> m_weightsRed;     
    // current test source reduced model map
    std::vector<float> m_targetRedModel;


    // log-likelihood for the reference fit
    double m_loglike_ref;

    // this is the list of models for the current fit
    // it may or may not include the test source model
    std::vector<const std::vector<float>* > m_currentModels;
    // Vector of sources that are currently free
    std::vector<bool> m_currentFreeSources;
    // this is the fixed component for the current fit
    std::vector<float> m_currentFixed;
    // this is the set of refrence values for the paremeters in the current fit
    std::vector<float> m_currentRefValues;
    // these are the indices of the components in the master list
    std::vector<int> m_currentSourceIndices;
    // this is the index of the test source in the current fit
    int m_currentTestSourceIndex;
    // these are the initial parameters for the current fit
    CLHEP::HepVector m_initPars;
    // these are the current parameters for the current fit
    CLHEP::HepVector m_currentPars;
    // this is the current covarience martix for the current fit
    CLHEP::HepSymMatrix m_currentCov;
    // this is the current gradient for the current fit
    CLHEP::HepVector m_currentGrad;
    // these are the prior to use for the current fit 
    FitScanMVPrior* m_prior_test;
    FitScanMVPrior* m_prior_bkg;

    // these are the global priors to use for all the broad band fits
    FitScanMVPrior* m_global_prior_test;
    FitScanMVPrior* m_global_prior_bkg;    

    // these are the initital priors, inherited from the input model
    FitScanMVPrior* m_init_prior_test;
    FitScanMVPrior* m_init_prior_bkg;        

    // this is the best-fit model for the current fit
    std::vector<float> m_currentBestModel;
    // this is the log-likelihood for the current fit
    double m_currentLogLike;
    // this is the estimated distance to minimum for the current fit
    double m_currentEDM;

    // this is the energy bin for the current fit
    // -1 means fit all ranges
    size_t m_firstEnergyBin;
    size_t m_lastEnergyBin;

    // this are the overall bin ranges for the current fit
    size_t m_firstBin;
    size_t m_lastBin;
    
  };


  /* A class to perform a series of related fits in a single ROI with a test source 

     This class can scan over:
     
     1) The test source location (a so-called TSMap)
     2) Energy bins (i.e., to make an SED)
     3) The normalization of the test source (i.e., to make to so-called "Castro" plot)

  */

  class FitScanner {

  public:

    /* Utility function to build an evtbin::Binner object with the energy bins */
    static evtbin::Binner* buildEnergyBinner(const std::vector<double>& energies);


    static astro::SkyProj* buildSkyProj(const std::string &projName,
					const astro::SkyDir& dir,
					double pixSize,
					int nPix, bool galactic=false);
    

    typedef enum  { 
      PRIMARY_HDU = 0x1,
      FITS_IMAGE = 0x2,
      FITDATA_TABLE = 0x4,
      SCANDATA_TABLE = 0x8
    } Out_Location;

  public:
    
    // C'tor from WCS grid of directions
    FitScanner(BinnedLikelihood& binnedLike,
	       optimizers::Optimizer& optimizer,
	       const astro::SkyProj& proj,
	       int nx, int ny);

    // C'tor from WCS grid of directions
    FitScanner(SummedLikelihood& summedLike,
	       optimizers::Optimizer& optimizer,
	       const astro::SkyProj& proj,
	       int nx, int ny);

    // C'tor from HEALPix region set of directions
    FitScanner(BinnedLikelihood& binnedLike,
	       optimizers::Optimizer& optimizer,
	       const astro::HealpixProj& proj,
	       const std::string& region);

    // C'tor from HEALPix region set of directions
    FitScanner(SummedLikelihood& summedLike,
	       optimizers::Optimizer& optimizer,
	       const astro::HealpixProj& proj,
	       const std::string& region);

    

    // D'tor, does cleanup
    virtual ~FitScanner() throw();
    
    /* Build a TS map.
       This scans over the directions and calculates the Test Statistics w.r.t. the null 
       hypothesis for each direction 

       This actually just calls run_tscube with doSED and doNorm set to false.
 
       covScale_bb   : Scale factor to apply to global fitting cov. matrix in broadband fits ( < 0 -> fixed )
       tol           : Critetia for fit convergence (estimated vertical distance to min < tol )
       maxIter       : Maximum number of iterations for the Newton's method fitter
       tolType       : Absoulte (0) or relative (1) criteria for convergence
       remakeTestSource : If true, recomputes the test source image (otherwise just shifts it)
       ST_scan_level : Level to which to do ST-based fitting (for testing)
       src_model_out : Name of XML file to write re-fit source model
       initLambda    : Initial damping parameter for step size calculation. (0 disables damping)
       useWeights    : Use the likelihood weights 

       returns 0 for success, or a error code
     */
    int run_tsmap(double covScale_bb = -1.0, 
		  double tol=1e-3, 
		  int tolType = 0, 
		  int maxIter = 30, 
		  bool remakeTestSource = false,
		  int ST_scan_level = 0,
		  std::string src_model_out = "",
		  double initLambda = 0.0,
		  bool useWeights = false);

     /* Build an SED.
	This calculates the spectrum as a function of energy
       and can also scan over the normalization.

       This actually just calls run_tscube with doTSMap set to false.

       nNorm         : Number of points in the likelihood v. normalization scan
       covScale_bb   : Scale factor to apply to global fitting cov. matrix in broadband fits ( < 0 -> fixed )
       covScale      : Scale factor to apply to broadband fitting cov. matrix in bin-by-bin fits ( < 0 -> fixed )
       normSigma     : Number of sigma to use for the scan range 
       tol           : Critetia for fit convergence (estimated vertical distance to min < tol )
       maxIter       : Maximum number of iterations for the Newton's method fitter       
       tolType       : Absoulte (0) or relative (1) criteria for convergence
       remakeTestSource : If true, recomputes the test source image (otherwise just shifts it)
       ST_scan_level : Level to which to do ST-based fitting (for testing)
       src_model_out : Name of XML file to write re-fit source model
       initLambda    : Initial damping parameter for step size calculation. (0 disables damping)
       useWeights    : Use the likelihood weights 

       returns 0 for success, or a error code
     */
    int run_SED(int nNorm = 10, double normSigma = 5.0, 
		double covScale_bb = -1.0, double covScale = -1.0, 
		double tol = 1e-3, int maxIter = 30, 
		int tolType = 0, 
		bool remakeTestSource = false,
		int ST_scan_level = 0,
		std::string src_model_out = "", 
		double initLambda = 0.0,
		bool useWeights = false);


    /* Build a TS cube.
       For each point in a TS Map this also calculate the spectrum as a function of energy
       and can also scan over the normalization

       doTSMap       : Scan over the grid of test directions
       doSED         : Compute the energy bin-by-bin fits
       nNorm         : Number of points in the likelihood v. normalization scan
       covScale_bb   : Scale factor to apply to global fitting cov. matrix in broadband fits ( < 0 -> fixed )
       covScale      : Scale factor to apply to broadband fitting cov. matrix in bin-by-bin fits ( < 0 -> fixed )
       normSigma     : Number of sigma to use for the scan range 
       tol           : Critetia for fit convergence (estimated vertical distance to min < tol )
       maxIter       : Maximum number of iterations for the Newton's method fitter
       tolType       : Absoulte (0) or relative (1) criteria for convergence
       remakeTestSource : If true, recomputes the test source image (otherwise just shifts it)
       ST_scan_level : Level to which to do ST-based fitting (for testing)
       src_model_out : Name of XML file to write re-fit source model
       initLambda    : Initial damping parameter for step size calculation. (0 disables damping)
       useWeights    : Use the likelihood weights 

       returns 0 for success, or a error code

       throws exceptions if some inconsistency in the models or data is detected.
     */
    int run_tscube(bool doTSMap=true, bool doSED=true, int nNorm = 10, double normSigma = 5.0, 
		   double covScale_bb = -1.0, double covScale = -1.0, 
		   double tol = 1e-3, int maxIter = 30, int tolType = 0, 
		   bool remakeTestSource = false,
		   int ST_scan_level = 0,
		   std::string src_model_out = "",
		   double initLambda = 0.0,
		   bool useWeights = false);

    /* Write the stored data to a FITS file */
    int writeFitsFile(const std::string& fitsFile,
		      const std::string& creator,
		      std::string fits_template = "",
		      bool copyGTIs=false) const;
      

    // Access Functions -----------

    // Likelihood object
    inline const FitScanModelWrapper* modelWrapper() const { return m_modelWrapper; }

    /* Optimizer used to do full fitting 

       In practice this will not be used unless it is specifcially requested.
     */
    inline const optimizers::Optimizer* optimizer() const { return m_opt; }

    // The projection used to bin the data
    inline const astro::ProjBase* proj() const { return m_proj; }

    // The currect direction for the test source
    inline const astro::SkyDir& testSourceDir() const { return m_testSourceDir; }

    /* The test source

       In practice this will not be moved during the scanning unless 
       full-fitting is requested.  

       Instead, the test source will be used to calculate a reference image, 
       and the image will then be moved by offsetting it by a the correct number of pixels
     */
    inline const Source* testSource() const { return m_testSource; }    

    // Name of the test source (for bookkeeping)
    inline const std::string& testSourceName() const { return m_testSourceName; }

    // The number of pixels we are scanning over
    int nPixels() const;

    // The number of energy bins we are scanning over
    int nEBins() const;
    
    // The number of normalization points we are scanning over
    int nNorms() const;


    // Modifiers (use with extreme caution) -----------

    /* Use a powerlaw point source 
       
       Note: this builds the initial image of the test source in the center of the pixel rounding down from the
       center of the counts map.  For maps with even numbers of pixels the ScienceTools convention is to put
       the reference direction on a pixel edge.  In that case this would build the test source offset
       by 1/2 pixel down.  For maps with odd numbers of pixels this would just build the test source 
       in the center of the central pixel.

       This doesn't matter if the remakesrc option is true, in that case the test source will be remade in the center of 
       each pixel in the output scan.  
       
       On the other hand, if the remakesrc option is false, the initial image will be simply be shifted, and then it does matter
       where it was built.
    */    
    int setPowerlawPointTestSource(optimizers::FunctionFactory& factory, 
				   double index=2.0);


    /* Use a source with a premade spectrum 
       
       Note: this builds the initial image of the test source in the center of the pixel rounding down from the
       center of the counts map.  For maps with even numbers of pixels the ScienceTools convention is to put
       the reference direction on a pixel edge.  In that case this would build the test source offset
       by 1/2 pixel down.  For maps with odd numbers of pixels this would just build the test source 
       in the center of the central pixel.

       This doesn't matter if the remakesrc option is true, in that case the test source will be remade in the center of 
       each pixel in the output scan.  
       
       On the other hand, if the remakesrc option is false, the initial image will be simply be shifted, and then it does matter
       where it was built.
    */ 
    int setPointTestSource(optimizers::Function& spectrum);
				   
    /* Use a source from the model 

       Note: this take the image of the test source as is.  To get the same results as
       with setPowerlawPointTestSource you either need to place the test source at the center of a pixel or
       set shiftToPixelCenter to true.

       For maps with even numbers of pixels the ScienceTools convention is to put
       the reference direction on a pixel edge.  In that case this would build the test source offset
       by 1/2 pixel down. 

       This doesn't matter if the remakesrc option is true, in that case the test source will be remade in the center of 
       each pixel in the output scan.  
       
       On the other hand, if the remakesrc option is false, the initial image will be simply be shifted, and then it does matter
       where it was built. 
    */
    int setTestSourceByName(const std::string& sourceName, bool shiftToPixelCenter=true);

    /* Use a generic source 

       Note: this take the image of the test source as is.  To get the same results as
       with setPowerlawPointTestSource you need to place the test source at the center of a pixel.
       For maps with even numbers of pixels the ScienceTools convention is to put
       the reference direction on a pixel edge.  In that case this would build the test source offset
       by 1/2 pixel down. 

       This doesn't matter if the remakesrc option is true, in that case the test source will be remade in the center of 
       each pixel in the output scan.  
       
       On the other hand, if the remakesrc option is false, the initial image will be simply be shifted, and then it does matter
       where it was built. 
    */
    int setTestSource(Source& source, bool owned=true);


    // Debbugging stuff
    inline int verbose_null() const { return m_verbose_null; }
    inline int verbose_bb() const { return m_verbose_bb; }
    inline int verbose_scan() const { return m_verbose_scan; }
    inline bool writeTestImages() const { return m_writeTestImages; }
    inline bool useReduced() const { return m_useReduced; }

    inline void set_verbose_null(int val) { m_verbose_null = val; }
    inline void set_verbose_bb(int val) { m_verbose_bb = val; }
    inline void set_verbose_scan(int val) { m_verbose_scan = val; }
    inline void set_writeTestImages(bool val) { m_writeTestImages = val; }
    inline void set_useReduced(bool val) { m_useReduced = val; }

  protected:

    /* This adds the test source to the source model */
    int addTestSourceToModel();

    /* This removes the test source from the source model */
    void removeTestSourceFromModel();

    /* Set the direction of the test source, based on the loop parameters */
    int setTestSourceDir(int ix, int iy);

    /* This does the baseline fit
       i.e., the fit without the test source */
    int baselineFit(double tol = 1e-3, int tolType = 0);

    /* This does the baseline fit with Newton's Method,
       for the normalization parameters only */
    int baselineFit_Newton(double tol = 1e-3, int maxIter = 30, double initLambda = 0.0);

    /* This does the broadband fit
       i.e., the fit with the source across the entire energy range */
    int fitTestSourceBroadband(double tol = 1e-3, int tolType = 0);   

    /* This does the sed fitting with Newton's method,                                                       
       for the normalization paramters only */
    int sed_binned_newton(int nnorm, double normSigma,
			  double constrainScale,
			  std::vector<double>& norm_mles,
			  std::vector<double>& pos_errs,
			  std::vector<double>& neg_errs,
			  std::vector<double>& logLike_mles,
			  std::vector<double>& uls,
			  std::vector<int>& sed_fit_status,
			  std::vector<std::vector<double> >& norms,
			  std::vector<std::vector<double> >& logLikes);

    /* Build and cache an image of the test source */
    int buildTestModelCache();
    
    /* Build an n-dimensional histogram based on the loop parameters */
    HistND* buildHist(const std::string& name, 
		      bool do_pix = true,
		      bool do_energy = false,
		      bool do_norm = false);        

    /* Write a histogram as a FITS image */
    int writeFitsImage(const std::string& fitsFile,
		       const std::string& extName,
		       const HistND& hist) const;
    
    /* Write a histogram as a FITS table */
    int writeFitsTable_byPixel(const std::string& fitsFile,
			       const std::string& extName,
			       const std::vector<std::pair<std::string,std::pair<HistND*,std::string> > >& colData) const;

    /* write the energy bins */
    int writeFits_EnergyBins(const std::string& fitsFile) const;

    /* write the Good time intervals */
    int writeFits_GTIs(const std::string& fitsFile) const;    
    
    /* write the baseline fit info */
    int writeFits_Baseline(const std::string& fitsFile) const;    


    /* Convert the dimension string to the format expected by FITS */
    bool convertDimString(const std::string& inString,
			  std::string& outString, 
			  bool do_pix = false,
			  bool do_energy = true,
			  bool do_norm = true) const;

    /* set the TDIM keyword */
    void setDimKeyword(tip::Header& header,
		       int icol,
		       const std::string& dimString) const;

    /* delete the test model caches */
    void deleteTestModelCaches();

  private:

    // The log-likelihood object (also the source model)
    FitScanModelWrapper* m_modelWrapper;

    // The optimizer
    optimizers::Optimizer* m_opt;    

    // The projection (only needed if we are looping over directions)
    const astro::ProjBase* m_proj;

    // The direction of the test source (changes if we are looping over directions)
    mutable astro::SkyDir m_testSourceDir;  
 
    // Binners for the various loop parameters
    const evtbin::Binner* m_dir1_binner;
    const evtbin::Binner* m_dir2_binner;
    const evtbin::Binner* m_energy_binner;
    const evtbin::Binner* m_norm_binner;

    // The test source 
    Source* m_testSource;

    // Do we own the test source
    bool m_testSourceOwned;

    // The name of the test source (useful for bookkeeping)
    std::string m_testSourceName;

    // These are the output data from the scanning
    std::vector< std::pair< std::string,std::pair<HistND*,std::string> > > m_scanData;

    // These tell us where to write the outputs
    std::map< std::string, unsigned int > m_outLocs;

    bool m_scanHasTSMap;
    // these are the parameters from the baseline fit
    CLHEP::HepVector m_baselinePars;
    // this is the covarience martix from the baseline fit
    CLHEP::HepSymMatrix m_baselineCovs;

    // This is where we cache all the info about the model and the fits
    FitScanCache* m_cache;

    // This is what we use to move around the image of the test source
    std::vector<TestSourceModelCache*> m_testSourceCaches;    

    // For debugging
    int m_verbose_null;
    int m_verbose_bb;
    int m_verbose_scan;
    bool m_writeTestImages;
    bool m_useReduced;

  };

}

#endif // Likelihood_CountsMapHealpix_h
