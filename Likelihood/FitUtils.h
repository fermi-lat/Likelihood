/**
 * @file FitUtils.h
 * @brief Functions to perform likelihood scans
 * @author E. Charles
 *
 *  This file contains a number of functions useful for likelihood scans.
 *
 *  These also include a number of vector manipulations.  In principle these could be 
 *  replaced with builtin C++ functionality, or functions in libraries like boost or thrust.
 *
 *  The second point is that many of these functions assume a specify type of vector (std::vector<float>).
 *  Although this could be replaced with templated code, the use of floats here is intentional.   For the 
 *  specific computations that I am doing, using doubles would be wasteful both computationally and in terms 
 *  of data volume.  
 *
 *  The third point is that in many cases I've chosen to used iterator, which makes the code somewhat less readable.
 *  This is purely for speed of execution.  I've tried to document what the actual function does in each case.
 *  
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/FitUtils.h,v 1.5 2016/02/05 22:31:12 echarles Exp $
 */

#ifndef Likelihood_FitUtils_h
#define Likelihood_Fitutils_h

#include <vector>
#include <map>
#include <string>
#include <list>

namespace CLHEP {
  class HepVector;
  class HepSymMatrix;
}

namespace Likelihood {
  
  class BinnedLikelihood;
  class FitScanMVPrior;
  class Source;
  class SourceMap;
  class SourceModel;

  namespace FitUtils {
    
    /* Fill an STL vector with values from a CLHEP vector */
    void Vector_Hep_to_Stl(const CLHEP::HepVector& hep,
			   std::vector<float>& stl);

    /* Fill an STL vector with values from a CLHEP Matrix */
    void Matrix_Hep_to_Stl(const CLHEP::HepSymMatrix& hep,
			   std::vector<float>& stl);
    
    /* Fill a CLHEP vector with values from an STL vector */
    void Vector_Stl_to_Hep(const std::vector<float>& stl,
			   CLHEP::HepVector& hep);
    
    /* Sum from start to stop and put them into value */
    void sumVector(std::vector<float>::const_iterator start,
		   std::vector<float>::const_iterator stop,
		   float& value);

    /* Set everything from start to stop equal to value */
    void setVectorValue(const float& val,
			std::vector<float>::iterator start,
			std::vector<float>::iterator stop);

    /* Add two vectors with optional factors.
       
       out = fact1 * v1 + fact2 * v2
       
       This checks that the distances: 
       ( stop1 - start1 ),  
       ( stop2 - start2 ),
       ( out_stop - out_start ) 
       are all the same and throws and exception if they are not.
    */
    void vectorAdd(std::vector<float>::const_iterator start1, 
		   std::vector<float>::const_iterator stop1, 
		   std::vector<float>::const_iterator start2, 
		   std::vector<float>::const_iterator stop2, 
		   std::vector<float>::iterator out_start,
		   std::vector<float>::iterator out_stop,
		   float fact1 = 1., float fact2 = 1.);

    /* Multiply two vectors.
       
       out = v1 * v2
       
       This checks that the distances: 
       ( stop1 - start1 ),  
       ( stop2 - start2 ),
       ( out_stop - out_start ) 
       are all the same and throws and exception if they are not.
    */
    void vectorMultiply(std::vector<float>::const_iterator start1, 
			std::vector<float>::const_iterator stop1, 
			std::vector<float>::const_iterator start2, 
			std::vector<float>::const_iterator stop2, 
			std::vector<float>::iterator out_start,
			std::vector<float>::iterator out_stop);
    
    /* Multiply in place vector by a        
	v = scalar * v       
    */  
    void multipyByScalar(std::vector<float>::iterator vstart,
			 std::vector<float>::iterator vstop,
			 double scalar);
    
    /* Subtract a vector from another vector.
       
       out = v1 - v2
       
       This checks that the distances: 
       ( stop1 - start1 ),  
       ( stop2 - start2 ),
       ( out_stop - out_start ) 
       are all the same and throws and exception if they are not.
    */   
    void vectorSubtract(std::vector<float>::const_iterator start1, 
			std::vector<float>::const_iterator stop1, 
			std::vector<float>::const_iterator start2, 
			std::vector<float>::const_iterator stop2, 
			std::vector<float>::iterator out_start,
			std::vector<float>::iterator out_stop);
      
    /* Calculate the inner product of two vectors and put them into value.

       value = Sum_i v1_i v2_i

       This checks that the distances: 
       ( stop1 - start1 ),  
       ( stop2 - start2 )
       are the same and throws and exception if they are not.
     */  
    void innerProduct(std::vector<float>::const_iterator start1, 
		      std::vector<float>::const_iterator stop1, 
		      std::vector<float>::const_iterator start2, 
		      std::vector<float>::const_iterator stop2, 
		      float& value);
 
    /* Calculate the inner product of three vectors and put them into value.

       value = Sum_i v1_i v2_i v3_i 

       This checks that the distances: 
       ( stop1 - start1 ),  
       ( stop2 - start2 ),
       ( stop3 - start3 )
       are the same and throws an exception if they are not.
    */  
    void innerProduct(std::vector<float>::const_iterator start1, 
		      std::vector<float>::const_iterator stop1, 
		      std::vector<float>::const_iterator start2, 
		      std::vector<float>::const_iterator stop2, 
		      std::vector<float>::const_iterator start3, 
		      std::vector<float>::const_iterator stop3, 
		      float& value);

    /* Calculate the fractional difference between data and model 

       out = ( data - model ) / data

       This checks that  data, model and output are all the same size 
       and throws an exception if they are not.

       This will throw an exception if model <= 0 in a bin where data != 0
    */  
    void fracDiff(std::vector<float>::const_iterator data_start, 
		  std::vector<float>::const_iterator data_stop,
		  std::vector<float>::const_iterator model_start,
		  std::vector<float>::const_iterator model_stop,
		  std::vector<float>::iterator out_start,
		  std::vector<float>::iterator out_stop);
 
    /* Calculate the data over the model squared  

       out = data / (model*model)

       This checks that  data, model and output are all the same size 
       and throws an exception if they are not.

       This will throw an exception if model <= 0 in a bin where data != 0
    */  
    void data_over_model2(std::vector<float>::const_iterator data_start, 
			  std::vector<float>::const_iterator data_stop,
			  std::vector<float>::const_iterator model_start,
			  std::vector<float>::const_iterator model_stop,
			  std::vector<float>::iterator out_start,
			  std::vector<float>::iterator out_stop);
        
    /* reshape a 1D vector in to 2D vector 

       This will throw an exception if n_i * n_j is not equal to the size of invect
     */
    void reshapeVector(const std::vector<float>& invect,
		       size_t n_i, size_t n_j,
		       std::vector<std::vector<float> >& outvect);


    /* Get the symetric error from the positive and negative errors */
    double symmetricError(double pos_err, double neg_err);
    
    /* Sum the model components * normalization to get the total model counts per bin
       
       total_i = fixed_i + templates_ij * norms_j

       Note that this gives you the option of only summing over part of the vector, using
       firstBin and lastBin.
    */    
    void sumModel(const CLHEP::HepVector& norms,
		  const std::vector<const std::vector<float>* >& templates,
		  const std::vector<float>& fixed,
		  std::vector<float>& total,
		  size_t firstBin = 0,
		  size_t lastBin = 0);
    
    /* Cacluate the gradiant and Hessian for a fit in which only the normalization of the
       components are allowed to vary.

       The elements are:

       g_a = Sum ( ( data - model ) / model ) * templates_a
       h_ab = Sum ( ( data / ( model*model ) ) * templates_a * templates_b

       Note that this gives you the option of only summing over part of the vector, using
       firstBin and lastBin.
    */      
    void getGradientAndHessian(const std::vector<float>& data,
			       const CLHEP::HepVector& norms,
			       const std::vector<const std::vector<float>* >& templates,
			       const std::vector<float>& fixed,
			       const FitScanMVPrior* prior,
			       std::vector<float>& model,
			       CLHEP::HepVector& gradient,
			       CLHEP::HepSymMatrix& hessian,
			       size_t firstBin = 0,
			       size_t lastBin = 0,
			       int verbose = 0);


    /* Invert the Hessian to get the covariance matrix.
       Use the gradient vector and the covariance matrix to the delta for the next iteration.
       Calculate the estimated vertical distance to the minimum using the delta and the gradient.
       
       V = h^-1
       delta_a = Sum V_ab * g_b
       edm = Sum delta_a * g_a
    */       
    int getDeltaAndCovarAndEDM(const CLHEP::HepSymMatrix& hessian,
			       const CLHEP::HepVector& gradient,
			       const CLHEP::HepVector& norms,
			       CLHEP::HepSymMatrix& covar,
			       CLHEP::HepVector& delta,
			       double& edm);
    
    /* Extract the set of bins that are non-zero into parallel vectors of indices
       this is done to make it faster to iterate over sparse data

       dataVect:      Input data map
       npix:          Number of pixels (used to insert seperator bins)
       nonZeroBins:   The indices of the non-zero bins
       dataRed:       Reduced data vector, include bins with 0 as separators between energies
       energyBinStartIdxs:  Indicies in the reduced vector marking the start of the energy layers

    */
    void extractNonZeroBins(const std::vector<float>& dataVect,
			    int npix,
			    std::vector<int>& nonZeroBins,
			    std::vector<float>& dataRed,
			    std::vector<int>& energyBinStopIdxs);
      



    /* Extract the bins that have non-zero counts, and reduce the correspond model maps

       dataVect:      Input data map
       model:         Input model template
       modelRed:      Reduced model template
       
       Note that modelRed vectors has a several bins, 
       corresponding to the sum of all of the zero counts bins in each energy layer.

       This allows use to re-use the other fitting functions in this file
     */
    void sparsifyModel(const std::vector<int>& nonZeroBins,
		       const std::vector<float>& model,
		       std::vector<float>& modelRed);


    /* Calculate and return the Poission negative log likelihood
       
       retVal = Sum  data_i - data_i * log(model_i) 

       For speed, the log is not executed for bins where data == 0
     */
    double negativeLogLikePoisson(std::vector<float>::const_iterator data_start, 
				  std::vector<float>::const_iterator data_stop,
				  std::vector<float>::const_iterator model_start,
				  std::vector<float>::const_iterator model_stop);
    

    /* Fit the normalization using Newton's method
       
       data:          The observed data
       initNorms:     Initial values of the normalizations
       templates:     Templates of the free sources
       fixed:         Template with the sum of all the fixed sources
       prior:         If provided, a multivariate prior on the fit
       tol:           Tolerance for estimating covergence.  Done when edm < tol
       maxIter:       Maximum number of iterations before failing the fit
       norms:         Filled with the fit results for the normalizations
       covar:         Filled with the covariance matrix from the fit
       gradient:      Filled with the gradient at the best fit point
       model:         Filled with the best-fit model
       edm:           Filled with the estimated (veritcal) distance to minimum
       logLikeVal:    Filled with the log likelihood at the best-fit point
       firstBin:      First bin to use
       lastBin:       Last bin to use ( 0 -> end )
       verbose:       0 : none; 1 : start & converged; 2 : params;  3 : matrices
    */
    int fitNorms_newton(const std::vector<float>& data,
			const CLHEP::HepVector& initNorms,
			const std::vector<const std::vector<float>* >& templates,
			const std::vector<float>& fixed,
			const FitScanMVPrior* prior,
			double tol, int maxIter,
			CLHEP::HepVector& norms,
			CLHEP::HepSymMatrix& covar,
			CLHEP::HepVector& gradient,
			std::vector<float>& model,
			double& edm,
			double& logLikeVal,
			size_t firstBin = 0, 
			size_t lastBin = 0,
			int verbose = 0);
    

    /* Fit the log of the normalization using Newton's method
       
       data:          The observed data
       initNorms:     Initial values of the normalizations
       templates:     Templates of the free sources
       fixed:         Template with the sum of all the fixed sources
       prior:         If provided, a multivariate prior on the fit
       tol:           Tolerance for estimating covergence.  Done when edm < tol
       maxIter:       Maximum number of iterations before failing the fit
       norms:         Filled with the fit results for the normalizations
       covar:         Filled with the covariance matrix from the fit
       gradient:      Filled with the gradient at the best fit point
       model:         Filled with the best-fit model
       edm:           Filled with the estimated (veritcal) distance to minimum
       logLikeVal:    Filled with the log likelihood at the best-fit point
       firstBin:      First bin to use
       lastBin:       Last bin to use ( 0 -> end )
       verbose:       0 : none; 1 : start & converged; 2 : params;  3 : matrices
    */
    int fitLogNorms_newton(const std::vector<float>& data,
			   const CLHEP::HepVector& initNorms,
			   const std::vector<const std::vector<float>* >& templates,
			   const std::vector<float>& fixed,
			   const FitScanMVPrior* prior,
			   double tol, int maxIter,
			   CLHEP::HepVector& norms,
			   CLHEP::HepSymMatrix& covar,
			   CLHEP::HepVector& gradient,
			   std::vector<float>& model,
			   double& edm,
			   double& logLikeVal,
			   size_t firstBin = 0, 
			   size_t lastBin = 0,
			   int verbose = 0);
    

    /* Extract a vector of spectral normalization values from a Source object

       source:    The source object
       energies:  The energies at which to evalute the spectrum
       specVals:  Filled with the normalization values at the input energies
     */
    void extractSpectralVals(const Source& source,
			     const std::vector<double>& energies,
			     std::vector<double>& specVals);

    
    /* Extract a vector of spectral normalization values from a Source object

       source:    The source object
       energies:  The energies at which to evalute the spectrum
       nPreds:  Filled with the normalization values at the input energies
     */
    void extractNPreds(const Source& source,
		       const std::vector<double>& energies,
		       std::vector<double>& nPreds);


    /* Extract the predicted counts model from a Source object

       source:    The source object
       logLike:   The BinnedLikelihood object with the template counts map and energy binning
       model:     Predicted counts model
       rescaleToNormOne:  If true the model is rescaled to the value it would have is the normalization parameter was 1.0
     */
    void extractModelFromSource(const Source& source,
				const BinnedLikelihood& logLike,
				std::vector<float>& model,
				bool rescaleToNormOne = false);

    /* Extract a vector of spectral normalization values from a SourceMap object

       sourceMap:   The source object
       energies:    The energies at which to evalute the spectrum
       specVals:    The normalazation values at the input energies
       modelCounts: Filled with the predicted model count
       
       Note, the output model counts vector has one less energy plane that the 
       energies and specVals vectors.   This is because we have integrated 
       over the energy bins in making the modelCounts vector.
     */ 
    void extractModelCounts(const SourceMap& sourceMap,
			    const std::vector<double>& energies,
			    const std::vector<double>& specVals,
			    std::vector<float>& modelCounts);

    /* Extract all of the predicted counts models for free sources from a 
       BinnedLikelihood object

       logLike:     The BinnedLikelihood object
       test_name:   Name of the test source
       templates:   Models for all of the free sources except the test source
       fixed:       Summed model for all of the fixed sources
       test_source_model: Model for all the test source
       refPars:     Values of the normalizations of the free sources.       
     */ 
    void extractModels(const BinnedLikelihood& logLike,
		       const std::string& test_name,
		       std::vector<std::vector<float> >& templates,		       
		       std::vector<float>& fixed,
		       std::vector<float>& test_source_model,
		       std::vector<float>& refPars);

    /* Refactors the free and fixed model components

       templates_in  :  Master list of input model templates
       fixed_in      :  Fixed component in input model
       scales_in     :  Scale factors to apply to input model templates
       freePars      :  Vector marking free sources
       test_source_model :  Test source model, null ptr -> do not include test source
       templates_out :  Vector of pointers to model templates for free sources
       fixed_out     :  Sum of fixed model components
       scales_out    :  Scale factors for free components
     */
    void refactorModels(const std::vector<std::vector<float> >& templates_in,
			const std::vector<float>& fixed_in,
			const std::vector<float>& scales_in,
			const std::vector<bool>& freePars,
			const std::vector<float>* test_source_model,
			std::vector<const std::vector<float>* >& templates_out,
			std::vector<float>& fixed_out,
			std::vector<float>& scales_out);
 
    /* Fit the normalization using Newton's method

       logLike:       The BinnedLikelihood object
       test_name:     Name of the test source
       prior:         If provided, a multivariate prior on the fit
       tol:           Tolerance for estimating covergence.  Done when edm < tol
       maxIter:       Maximum number of iterations before failing the fit
       norms:         Filled with the fit results for the normalizations
       covar:         Filled with the covariance matrix from the fit
       gradient:      Filled with the gradient at the best fit point
       model:         Filled with the best-fit model
       edm:           Filled with the estimated (vertical) distance to minimum
       logLikeVal:    Filled with the log likelihood at the best-fit point
       firstBin:      First bin to use
       lastBin:       Last bin to use ( 0 -> end )
    */
    int fitModelNorms_newton(const BinnedLikelihood& logLike,
			     const std::string& test_name,
			     double tol, int maxIter,
			     CLHEP::HepVector& norms,
			     CLHEP::HepSymMatrix& covar,
			     CLHEP::HepVector& gradient,
			     std::vector<float>& model,
			     double& edm,
			     double& logLikeVal,
			     size_t firstBin = 0, size_t lastBin = 0);    
  
    /* Print a CLHEP vector to std::cout */
    void printVector(const std::string& name,
		     const CLHEP::HepVector& vect);

    /* Print a CLHEP SymMatrix to std::cout */
    void printSymMatrix(const std::string& name,
			const CLHEP::HepSymMatrix& mat);

  };

} // namespace Likelihood

#endif // Likelihood_FitUtils_h
