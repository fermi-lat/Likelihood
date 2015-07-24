/**
 * @file ScanUtils.h
 * @brief Functions to perform likelihood scans
 * @author E. Charles
 *
 *
 *  This file contains a number of functions useful for likelihood scans.
 *
 *  These are basically utility function to do very specific things with the
 *  BinnedLikelihood and SummedLikelihood interfaces
 *
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/ScanUtils.h,v 1.1 2015/07/17 18:41:54 echarles Exp $
 */

#ifndef Likelihood_ScanUtils_h
#define Likelihood_Scanutils_h

#include <vector>
#include <string>
#include <list>

namespace optimizers {
  class Optimizer;
}

namespace Likelihood {

  class Source;
  class SourceModel;
  class FitScanModelWrapper;
  
  namespace ScanUtils {
    
    /* Freeze all the parameters of a given source */
    void freezeSourceParams(Source& source);
    
    /* Freeze all the parameters for all the sources in a source model */
    void freezeAllParams(SourceModel& srcModel);

    /* Free or Freeze all the normalizations of a list of source in a source model */
    void freeNormalizations(SourceModel& srcModel,
			    const std::list<std::string>& freeSources,
			    bool free=true);
    
    /* Apply Gaussian Prior constraints on the normalization in a source model 
       using a set of input values and covariances 

       srcModel         : The SourceModel object
       srcNameList      : List of source to constrain
       vals             : Central values for constraints
       covs             : Covariance matrix between constraints    
    */
    void constrainSourceNorms_valuesAndCovs(SourceModel& srcModel,
					    const std::list<std::string>& srcNameList,
					    const std::vector<double>& vals,
					    const std::vector<std::vector<double> >& covs);
    
    /* Apply Gaussian Prior constraints on the normalization in a source model 
       using the current fit covariance matrix

       srcModel         : The SourceModel object
       srcNameList      : List of source to constrain
       covScale         : An overall scale factor to apply to the widths of the priors
    */
    void constrainSourceNorms_current(SourceModel& srcModel,
				      const std::list<std::string>& srcNameList,
				      const double& covScale);

    /* Scan the likelihood as a function of the normalization of the test source
     */
    void scan_norm_binned(FitScanModelWrapper& modelWrapper,
			  const std::string& signal_name,			  
			  optimizers::Optimizer& optimizer,
			  double tol, int tolType,
			  size_t kmin, size_t kmax,
			  size_t nnorms,
			  std::vector<double>& norms,
			  double& norm_mle,
			  double& logLike_mle,	
			  std::vector<double>& logLikes);
  
    /* Scan the likelihood as a function of the normalization of the test source
     */
    void sed_binned(FitScanModelWrapper& modelWrapper,
		    const std::string& signal_name,			  		    
		    optimizers::Optimizer& optimizer,
		    double tol,int tolType,
		    size_t nnorms,
		    std::vector<std::vector<double> >& norms,
		    std::vector<double>& norm_mles,
		    std::vector<double>& logLike_mles,
		    std::vector<std::vector<double> >& logLikes);

  };

} // namespace Likelihood

#endif // Likelihood_ScanUtils_h
