/**
 * @file Convolve.cxx
 * @brief Functions to perform convolutions of HEALPix maps
 * @author E. Charles
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/Likelihood/src/ConvolveHealpix.cxx,v 1.2 2015/03/03 06:00:00 echarles Exp $
 */


#include "Likelihood/ScanUtils.h"

#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/Source.h"
#include "optimizers/OptimizerFactory.h"
#include "optimizers/Optimizer.h"
#include "optimizers/Parameter.h"


namespace Likelihood {

  namespace ScanUtils {

    void freezeSourceParams(Source& source) {
      std::vector<std::string> parNames;
      optimizers::Function& spec = source.spectrum();
      spec.getParamNames(parNames);
      for ( std::vector<std::string>::const_iterator itr = parNames.begin();
	    itr != parNames.end(); itr++ ) {
	spec.parameter(*itr).setFree(false);
      }      
    }

    void freezeAllParams(SourceModel& srcModel) {
      const std::map<std::string, Source *>& sources = srcModel.sources();
      for ( std::map<std::string, Source *>::const_iterator itr = sources.begin();
	    itr != sources.end(); itr++ ) {
	Source* aSrc = itr->second;
	freezeSourceParams(*aSrc);
      }
    }

    void freeNormalizations(SourceModel& srcModel,
			    const std::list<std::string>& freeSources,
			    bool free) {
      for ( std::list<std::string>::const_iterator itr = freeSources.begin();
	    itr != freeSources.end(); itr++ ) {
	Source* aSrc = srcModel.getSource(*itr);
	aSrc->spectrum().normPar().setFree(free);
      }
    }  
    
    void constrainSourceNorms_valuesAndCovs(SourceModel& srcModel,
					    const std::list<std::string>& srcNameList,
					    const std::vector<double>& vals,
					    const std::vector<std::vector<double> >& covs) {
      return;
    }
	
										    
    void constrainSourceNorms_current(SourceModel& srcModel,
				      const std::list<std::string>& srcNameList,
				      const double& covScale) {
      return;
    }

    void scan_norm_binned(BinnedLikelihood& logLike,
			  const std::string& signal_name,
			  optimizers::Optimizer& optimizer,
			  double tol, int tolType,
			  size_t kmin,size_t kmax,
			  size_t nnorms,
			  std::vector<double>& norms,
			  double& norm_mle,
			  double& logLike_mle,	
			  std::vector<double>& logLikes) {

      // Set the energy bounds
      std::pair<int, int> klims = logLike.klims();
      if ( klims.first != kmin ||
	   klims.second != kmax ) {
	logLike.set_klims(kmin,kmax);
      }

      double emin = logLike.energies()[kmin];
      double emax = logLike.energies()[kmax];

      // get the source and the normalization parameter
      Source* signal_source = logLike.getSource(signal_name);
      optimizers::Parameter& normPar = signal_source->spectrum().normPar();
      normPar.setFree(true);
      logLike.syncParams();

      // find the global min
      int status = optimizer.find_min(false,tol,tolType);
      
      // Move the errors into the parameters.   
      // FIXME, i'm not sure this is done correctly
      const std::vector<double>& uncertainies = optimizer.getUncertainty(false);
      const std::vector<std::vector<double> > cov = optimizer.covarianceMatrix();
      std::vector<optimizers::Parameter>& pars = logLike.parameters();
      size_t npar = pars.size();
      // count the number of free parameters
      size_t nfree(0);
      for ( size_t ipar(0); ipar < npar; ipar++ ) {
	optimizers::Parameter& par = pars[ipar];
	if ( par.isFree() ) {
	  par.setError(uncertainies[nfree]);
	  double check_1 = std::sqrt(cov[nfree][nfree]);
	  nfree++;
	}
      }
      if ( nfree != uncertainies.size() ) {
	std::cout << "Wrong number of free parameters " << nfree << ' ' << uncertainies.size() << std::endl;
      }
      
      // grab the log-likelihood, and various normalization factors
      // FIXME
      logLike_mle = logLike.value();
      norm_mle = normPar.getValue();
      double norm_err = uncertainies[0];

      double scan_min = std::max(1.0e-10,norm_mle - (10.*norm_err));
      double scan_max = std::max(2.5e-0, norm_mle + (10.*norm_err));

      double lin_step = (scan_max - scan_min) / float(nnorms-1);
      norms.clear();
      for ( int is(0); is < nnorms; is++ ) {
	norms.push_back(scan_min);
	scan_min += lin_step;	
      }      

      // clear out the result vector
      logLikes.clear();   

      std::pair<double, double> bounds = normPar.getBounds();
      normPar.setFree(false);

      // Loop over the normalization values
      for ( std::vector<double>::const_iterator itrNorm = norms.begin();
	    itrNorm != norms.end(); itrNorm++ ) {
	
	// Set the normalization parameter
	normPar.setValue(*itrNorm);
	logLike.syncParams();

	// if there are any other free parameters re-minimize 
	if ( nfree > 1 ) {
	  status = optimizer.find_min_only(true,tol,tolType);
	}
	
	// grap the delta log-likelihood value
	double logLike_val = logLike.value();
	logLikes.push_back(logLike_val);
      }
      // Clean up
      normPar.setFree(true);
      logLike.syncParams();
    }
    
    void sed_binned(BinnedLikelihood& logLike,
		    const std::string& signal_name,			  
		    optimizers::Optimizer& optimizer,
		    double tol, int tolType,
		    size_t nnorms,
		    std::vector< std::vector<double> >& norms,
		    std::vector<double>& norm_mles,
		    std::vector<double>& logLike_mles,
		    std::vector<std::vector<double> >& logLikes) {
      // get the number of energy bins
      int nebins = logLike.energies().size() - 1;

      norms.resize(nebins);
      norm_mles.resize(nebins);
      logLike_mles.resize(nebins);
      logLikes.resize(nebins);

      for ( int ie(0); ie < nebins; ie++ ) {
	scan_norm_binned(logLike,signal_name,
			 optimizer,tol,tolType,
			 ie,ie+1,nnorms,
			 norms[ie],
			 norm_mles[ie],logLike_mles[ie],
			 logLikes[ie]);
      }

      // Reset the klims
      logLike.set_klims(0,nebins+1);

    }

  } // namespace ConvolveHealpix
 
} // namespace Likelihood
