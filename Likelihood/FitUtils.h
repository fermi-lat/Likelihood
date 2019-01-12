/**
 * @file FitUtils.h
 * @brief Functions to perform likelihood scans
 * @author E. Charles
 *
 *  This file contains a number of functions useful for likelihood fitting and likelihood scans.
 *
 *  These also include a number of vector manipulations.  In principle these could be 
 *  replaced with builtin C++ functionality, or functions in libraries like boost or thrust.
 *
 *  The second point is that many of these functions assume a specify type of vector (std::vector<float>).
 *  Although this could be replaced with templated code, the use of floats here is intentional.   For the 
 *  specific computations that I am doing, using doubles would be wasteful both computationally and in terms 
 *  of data volume.  
 *
 *  The third point is that in many cases I've chosen to use iterators, which makes the code somewhat less readable.
 *  This is purely for speed of execution.  I've tried to document what the actual function does in each case.
 *  
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/FitUtils.h,v 1.14 2016/07/11 23:44:04 mdwood Exp $
 */

#ifndef Likelihood_FitUtils_h
#define Likelihood_Fitutils_h

#include <vector>
#include <map>
#include <string>
#include <list>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "Likelihood/Accumulator.h"

namespace CLHEP {
  class HepVector;
  class HepMatrix;
  class HepSymMatrix;
}

namespace optimizers {
  class Parameter;
}

namespace Likelihood {
  
  class BinnedLikelihood;
  class BinnedCountsCache;
  class FitScanMVPrior;
  class Source;
  class SourceMap;
  class SourceModel;
  class Drm;
  class Drm_Cache;

  namespace FitUtils {
    

    /// Integrates weights over a pixel to get the counts
    double pixelCounts_linearQuad(double emin, double emax, double y1, double y2);

    /// Integrates weights over a pixel to get the counts
    double pixelCounts_loglogQuad(double emin, double emax, double y1, double y2, double log_ratio);
  
    /// Expands a vector of energy bin edges by taking equal size steps in log space
    void expand_energies(std::vector<double> & energies, int edisp_bins);

    /* Compute the pseudo-inverse of a matrix from its singular value
       decomposition.  For a matrix A with SVD,

       A = U * S * V^t

       its psuedo-inverse is given by

       A^-1 = V * S^-1 * U^t

       where singular values in S^-1 are set to zero.

     */
    void svd_inverse(const CLHEP::HepMatrix& u,
		     const CLHEP::HepMatrix& v,
		     const CLHEP::HepVector& s,
		     CLHEP::HepMatrix& inverse);
    
    /* Solve a system of equations using singular value decomposition.  

       h * d = g

       where h is the hessian matrix and g is the gradient vector.
       Returns a solution vector delta and the SVD decomposition of
       the hessian.  The parameter eps sets the fractional threshold
       for setting singular values to zero.

     */
    int svd_solve(const CLHEP::HepSymMatrix& hessian,
		  const CLHEP::HepVector& gradient,
		  CLHEP::HepMatrix& u,
		  CLHEP::HepMatrix& v,
		  CLHEP::HepVector& s,
		  CLHEP::HepVector& delta,
		  double eps = 1E-16);
    
    gsl_vector * Vector_Hep_to_Gsl(const CLHEP::HepVector& hep);

    gsl_matrix * Matrix_Hep_to_Gsl(const CLHEP::HepMatrix& hep);

    /* Fill a CLHEP vector with values from a GSL vector */
    void Vector_Gsl_to_Hep(const gsl_vector * gsl,
			   CLHEP::HepVector& hep);

    /* Fill a CLHEP matrix with values from a GSL matrix */
    void Matrix_Gsl_to_Hep(const gsl_matrix * gsl,
			   CLHEP::HepMatrix& hep);

    /* Fill an STL vector with values from a CLHEP vector */
    void Vector_Hep_to_Stl(const CLHEP::HepVector& hep,
			   std::vector<float>& stl);

    /* Fill an STL vector with values from a CLHEP Matrix */
    void Matrix_Hep_to_Stl(const CLHEP::HepSymMatrix& hep,
			   std::vector<float>& stl);
    
    /* Fill a CLHEP vector with values from an STL vector */
    void Vector_Stl_to_Hep(const std::vector<float>& stl,
			   CLHEP::HepVector& hep);
    
    /* Fill a CLHEP matrix with values from an STL vector */
    void Matrix_Stl_to_Hep(const std::vector<float>& stl,
			   CLHEP::HepSymMatrix& hep);

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
    void multiplyByScalar(std::vector<float>::iterator vstart,
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
 

    /* Calculate the inner product of four vectors and put them into value.

       value = Sum_i v1_i v2_i v3_i v4_i

       This checks that the distances: 
       ( stop1 - start1 ),  
       ( stop2 - start2 ),
       ( stop3 - start3 )
       ( stop4 - start4 )
       are the same and throws an exception if they are not.
    */  
    void innerProduct(std::vector<float>::const_iterator start1, 
		      std::vector<float>::const_iterator stop1, 
		      std::vector<float>::const_iterator start2, 
		      std::vector<float>::const_iterator stop2, 
		      std::vector<float>::const_iterator start3, 
		      std::vector<float>::const_iterator stop3, 
		      std::vector<float>::const_iterator start4, 
		      std::vector<float>::const_iterator stop4, 
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
    
    
    /* Sum a set of vectorn
       
       total_i = fixed_i + Sum_j templates_ij 

       Note that this gives you the option of only summing over part of the vector, using
       firstBin and lastBin.
    */    
    void sumModel_Init(const std::vector<std::vector<float> >& templates,
		       const std::vector<float>& fixed,
		       std::vector<float>& total);
    
    

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
    
    /* Evaluate the total log-likelihood for the sum of data and prior terms. */
    double getLogLikelihood(const std::vector<float>& data,
			    const CLHEP::HepVector& norms,
			    const std::vector<const std::vector<float>* >& templates,
			    const std::vector<float>& fixed,
			    const FitScanMVPrior* prior,			       
			    const std::vector<float>* weights,    
			    std::vector<float>& model,
			    size_t firstBin = 0,
			    size_t lastBin = 0,
			    int verbose = 0);

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
			       const std::vector<float>* weights,
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

    int getDeltaAndCovarAndEDM(const CLHEP::HepSymMatrix& hessian,
			       const CLHEP::HepVector& gradient,
			       const CLHEP::HepVector& norms,
			       CLHEP::HepSymMatrix& covar,
			       CLHEP::HepVector& delta,
			       double& edm,
			       double lambda);

    int getDeltaAndCovarAndEDM_STL(const std::vector<float>& hessian,
				   const std::vector<float>& gradient,
				   const std::vector<float>& norms,
				   std::vector<float>& covar,
				   std::vector<float>& delta,
				   std::vector<double>& edm,
				   double lambda=-1.0);
    
    
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

    /* Extract the bins that have non-zero counts, and reduce the correspond model maps

       dataVect:      Input data map
       weights:       Input weights
       model:         Input summed model
       weightsRed:    Reduced weights
       
       Note that modelRed vectors has a several bins, 
       corresponding to the sum of all of the zero counts bins in each energy layer.

       This allows use to re-use the other fitting functions in this file
     */
    void sparsifyWeights(const std::vector<int>& nonZeroBins,
			 const std::vector<float>& weights,
			 const std::vector<float>& model,
			 std::vector<float>& weightsRed);


    /* Calculate and return the Poission negative log likelihood
       
       retVal = Sum  data_i - data_i * log(model_i) 

       For speed, the log is not executed for bins where data == 0
     */
    double logLikePoisson(std::vector<float>::const_iterator data_start, 
			  std::vector<float>::const_iterator data_stop,
			  std::vector<float>::const_iterator model_start,
			  std::vector<float>::const_iterator model_stop);
    
    /* Calculate and return the Poission negative log likelihood
       
       retVal = Sum  w_i * (data_i - data_i * log(model_i) )

       For speed, the log is not executed for bins where data == 0
     */
    double logLikePoisson(std::vector<float>::const_iterator data_start, 
			  std::vector<float>::const_iterator data_stop,
			  std::vector<float>::const_iterator model_start,
			  std::vector<float>::const_iterator model_stop,			 
			  std::vector<float>::const_iterator w_start,
			  std::vector<float>::const_iterator w_stop);


    /* Fit the normalization using Newton's method
       
       data:          The observed data
       initNorms:     Initial values of the normalizations
       templates:     Templates of the free sources
       fixed:         Template with the sum of all the fixed sources
       prior:         If provided, a multivariate prior on the fit
       tol:           Tolerance for estimating covergence.  Done when edm < tol
       maxIter:       Maximum number of iterations before failing the fit
       lambda:        Initial damping parameter for step size calculation. (0 disables damping)
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
			const std::vector<float>* weights,
			double tol, int maxIter, double lambda,
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
       lambda:        Initial damping parameter for step size calculation. (0 disables damping)
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
			   const std::vector<float>* weights,
			   double tol, int maxIter, double lambda,
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

    /* Extract an array of derivatives of the spectram from a Source object

       source:     The source object
       energies:   The energies at which to evalute the spectrum
       paramNames: The names of the params w.r.t. which to evaluate the derivaties
       deriveVals: Filled with the normalization values at the input energies
     */
    void extractSpectralDerivs(const Source& source,
			       const std::vector<double>& energies,
			       const std::vector<std::string>& paramNames,
			       std::vector<std::vector<double> >& derivVals);


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
    void extractModelCounts(SourceMap& sourceMap,
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
       weights:     Likelihood weights (optional)
       useUnitRefVals: Set the reference normalizations to one. 
     */ 
    void extractModels(const BinnedLikelihood& logLike,
		       const std::string& test_name,
		       std::vector<std::string>& freeSrcNames,
		       std::vector<std::vector<float> >& templates,		       
		       std::vector<float>& fixed,
		       std::vector<float>& test_source_model,
		       std::vector<float>& refPars,
		       std::vector<float>* weights = 0,
		       bool useUnitRefVals = false);

    
    /* Extract all of the predicted counts for all of the fixed sources
       BinnedLikelihood object

       logLike:     The BinnedLikelihood object
       test_name:   The name of the test source ( not added to model)
       fixed:       Summed model for all of the fixed sources
       latched:     An option vector of latched source ( not added to model )
     */ 
    void extractFixedModel(const BinnedLikelihood& logLike,
			   const std::string& test_name,
			   std::vector<float>& fixed,
			   const std::vector<std::string>* latched=0);
    

    /* Extract the prior on a parameter
       par               : The parameter in question
       centralVal        : Filled with the central value
       uncertainty       : Filled with the uncertainty 
       returns trun if Parameter has a prior
    */ 
    bool extractPrior(const optimizers::Parameter& par,
		      double& centralVal,
		      double& uncertainty);

    /* Extract the priors on the free source normalizations from the 
       BinnedLikelihood object

       logLike           : The BinnedLikelihood object
       freeSrcNames      : Names of all the free sources except the test source
       centralVals       : Filled with the central values from the priors
       uncertainties     : Filled with the untertainty values from the priors
       parHasPrior       : Filled with the flags showing which parameters have priors

       return true if any priors were extract, false otherwise
    */ 
    bool extractPriors(const BinnedLikelihood& logLike,
		       const std::vector<std::string>& freeSrcNames,
		       CLHEP::HepVector& centralVals,
		       CLHEP::HepVector& uncertainties,
		       std::vector<bool>& parHasPrior);

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
       lambda:        Initial damping parameter for step size calculation. (0 disables damping)
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
			     double tol, int maxIter, double lambda,
			     CLHEP::HepVector& norms,
			     CLHEP::HepSymMatrix& covar,
			     CLHEP::HepVector& gradient,
			     std::vector<float>& model,
			     double& edm,
			     double& logLikeVal,
			     size_t firstBin = 0, size_t lastBin = 0);    

    /* Get the quanities that appears in the log-log quadrature formula for the spectral term

       For a given energy bin, i, these are:
       spec_wts[i].first = spec[i] * energy[i] * log_energy_ratio[i] / 2.
       spec_wts[i].second = spec[i+1] * energy[i+1] * log_energy_ratio[i] / 2.

       spec       : The specturm (or spectral derivative) for the source in question
       dataCache  : Object with info about the binning
       spec_wts   : The specturm (or spectral derivative) weights for the source in question

    */  
    void get_spectral_weights(const std::vector<double>& spec,
			      const BinnedCountsCache& dataCache,
			      std::vector<std::pair<double, double> >& spec_weights);


    /* Get the range of energy bins to consider when computing the energy dispersion

	dataCache  : Object with info about the binning
        edisp_val  : Flag saying how to apply energy dispersion
	k          : The index of the current energy bin
	kmin       : Index of the lowest energy bin to consider
	kmax      : Index of the lowest energy bin to consider

    */  
    void get_edisp_range(const BinnedCountsCache& dataCache, int edisp_val,
			 size_t k,  
			 size_t& kmin, size_t& kmax);
	
		 
    /* Get the constants we need for the energy dispersion

	srcMap     : The source map for the source in question
 	dataCache  : Object with info about the binning
        edisp_val  : Flag saying how to apply energy dispersion
	k          : The current energy bin
	kmin       : Filled with the index of the lowest energy bin to loop over
	kmax       : Filled with the index of the highest energy bin to loop over
	edisp_col  : Filled with the energy disperson factors

	Note that edisp_col.size() == kmax - kmin
	I.e., only the factors for the bins we are looping over are extracted
    */  
    void get_edisp_constants(SourceMap& srcMap, const BinnedCountsCache& dataCache, 
			     int edisp_val, 
			     size_t k, size_t& kmin, size_t& kmax,
			     std::vector<double>& edisp_col);    

    
     /* Get the model counts contribution to a particular pixel 

	srcMap     : The SourceMap for the source in question
	spec_wts   : The specturm (or spectral derivative) weights for the source in question
	xi         : The energy disperion correction factor
	npix       : Number of pixes per energy layer
	kref       : Index of the energy bin in question
	ipix       : Pixel index (within the energy layer)

	returns the contribution
     */
    double model_counts_contribution(SourceMap& srcMap,
				     const std::vector<std::pair<double, double> > & spec_wts,
				     const double& xi, 
				     size_t npix, size_t kref, size_t ipix);
				      
     /* Get the model counts contribution including the energy dispersion

	srcMap     : The SourceMap for the source in question
	spec_wts   : The specturm (or spectral derivative) weights for the source in question
	edisp_col  : Energy dispersion factors
	ipix       : Index of the pixel in question
	npix       : Number of pixels per energy layer
	kmin       : Index of the first energy layer to consider true counts from
	kmax       : Index of the last energy layer to consider true counts from

	returns the contribution
     */
    double model_counts_edisp(SourceMap& srcMap,
			      const std::vector<std::pair<double, double> > & spec_wts,
			      const std::vector<double> & edisp_col,
			      size_t ipix, size_t npix, size_t kmin, size_t kmax);

    /* Get the total model counts contribution to a energy layer

	npred_vals     : The total npred for that energy layer
	npred_weights  : The weight factors for that energy layer
	spec_wts   : The specturm (or spectral derivative) weights for the source in question
	xi         : The energy dispersion correction factor
	kref       : The energy bin
	counts     : Filled with the counts contribution
	counts_wt  : Filled with the weighted counts contribution
     */
    void npred_contribution(const std::vector<double>& npred_vals,
			    const std::pair<double,double>& npred_weights,
			    const std::vector<std::pair<double, double> > & spec_wts,
			    const double& xi, 
			    size_t kref,
			    double& counts,
			    double& counts_wt);

    /* Get the model counts contribution to a energy layer

	npred_vals     : The total npred for that energy layer
	npred_weights  : The weight factors for that energy layer
	spec_wts   : The specturm (or spectral derivative) weights for the source in question
	edisp_col  : Energy dispersion factors
	kmin       : Index of the first energy layer to consider true counts from
	kmax       : Index of the last energy layer to consider true counts from
	counts     : Filled with the total model count contribution
	counts_wt  : Filled with the total weighted model count contribution
    */
    void npred_edisp(const std::vector<double> & npred_vals,
		     const std::vector<std::pair<double,double> > & npred_weights,
		     const std::vector<std::pair<double,double> > & spec_wts,
		     const std::vector<double> & edisp_col,
		     size_t kmin,
		     size_t kmax,
		     double& counts,
		     double& counts_wt);

     /* Add (or subtract) the counts for a source onto a vector 	
	This is used by several functions.

	modelCounts: The vector being added to.
	srcMap     : The SourceMap for the source in question
	dataCache  : Object with info about the binning
	edisp_val : How to apply the energy dispersion
	subtract   : If true, subtract from the vector.  	
     */     

    void addSourceCounts(std::vector<double> & modelCounts,
			 SourceMap& srcMap,
			 const BinnedCountsCache& dataCache,
			 int edisp_val,
			 bool subtract);
 
     /* Add (or subtract) the counts for a fxied source onto the vectors that 
	collect that info for fixed sources.
	This is used by several functions.

	fixed_counts_spec : The vector being added to.
	fixed_counts_spec_wt : The vector being added to.
	fixed_counts_spec_edisp : The vector being added to.
	fixed_counts_spec_edisp_wt : The vector being added to.
	srcMap     : The SourceMap for the source in question
	dataCache  : Object with info about the binning
	edisp_val : How to apply the energy dispersion
	subtract   : If true, subtract from the vector.  	
     */     
    void addFixedNpreds(std::vector<double>& fixed_counts_spec,
			std::vector<double>& fixed_counts_spec_wt,
			std::vector<double>& fixed_counts_spec_edisp,
			std::vector<double>& fixed_counts_spec_edisp_wt,
			SourceMap& srcMap,
			const BinnedCountsCache& dataCache,
			int edisp_val,
			bool subtract);

#ifndef SWIG
     /* Add (or subtract) the contributions to the derivatives from a single source.

	posDerivs  : The vector being added to.
	negDerivs  : The vector being added to.
	freeIndex  : The overall index of the first parameter for this source
	srcMap     : The SourceMap for the source in question
	data_over_model : A vector with the ratio of data to model in each filled pixel
	dataCache  : Object with info about the binning
	edisp_val  : How to apply the energy dispersion
	kmin       : Index of energy bin to start summation
	kmax       : Index of energy bin to stop summation
     */     
    void addFreeDerivs(std::vector<Kahan_Accumulator>& posDerivs,
		       std::vector<Kahan_Accumulator>& negDerivs,
		       long freeIndex,
		       SourceMap& srcMap,
		       const std::vector<double>& data_over_model, 
		       const BinnedCountsCache& dataCache,
		       int edisp_val,
		       size_t kmin, size_t kmax);
#endif //SWIG

          
    /* Add (or subtract) the counts for a source onto a vector.
       This version loops over all the pixels, not just the filled ones.

	modelMap   : The vector being added to.
	srcMap     : The SourceMap for the source in question
	dataCache  : Object with info about the binning
	use_mask   : Use the Weights mask
	edisp_val : How to apply the energy dispersion
     */     
    void updateModelMap(std::vector<float> & modelMap,
			SourceMap& srcMap,
			const BinnedCountsCache& dataCache,				       
			bool use_mask,
			int edisp_val);

    /* Print a CLHEP vector to std::cout */
    void printVector(const std::string& name,
		     const CLHEP::HepVector& vect);

    /* Print a CLHEP SymMatrix to std::cout */
    void printMatrix(const std::string& name,
		     const CLHEP::HepSymMatrix& mat);

    /* Print a CLHEP Matrix to std::cout */
    void printMatrix(const std::string& name,
		     const CLHEP::HepMatrix& mat);

  };

} // namespace Likelihood

#endif // Likelihood_FitUtils_h
