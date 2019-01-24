/**
 * @file Convolve.cxx
 * @brief Functions to perform convolutions of HEALPix maps
 * @author E. Charles
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/FitUtils.cxx,v 1.25 2016/09/14 20:11:02 echarles Exp $
 */


#include "Likelihood/FitUtils.h"

#include <cmath>

#include <stdexcept>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_linalg.h"

#include "Likelihood/Drm.h"
#include "Likelihood/Source.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/SourceModel.h"
// #include "Likelihood/CountsMapBase.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/FitScanner.h"
#include "Likelihood/WeightMap.h"

#include <vector>
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "optimizers/Parameter.h"
#include "optimizers/Function.h"
#include "optimizers/dArg.h"

namespace Likelihood {

  namespace FitUtils {

    /// Integrates weights over a pixel to get the counts
    double pixelCounts_linearQuad(double emin, double emax, double y1, double y2) {
      return (y1 + y2)*(emax - emin)/2.;
    }

    /// Integrates weights over a pixel to get the counts
    double pixelCounts_loglogQuad(double emin, double emax, double y1, double y2, double log_ratio) {
      return (y1*emin + y2*emax)/2.*log_ratio;
    }
  
    /// Expands a vector of energy bin edges by taking equal size steps in log space
    void expand_energies(std::vector<double> & energies, int edisp_bins) {      
      size_t nee = energies.size();

      if ( nee < 2 ) return;
      
      // Make the output vector and copy the input vector to the middle of it
      std::vector<double> out_vect(nee+2*edisp_bins, 0);
      std::copy(energies.begin(), energies.end(), out_vect.begin() + edisp_bins);
      
      // Get the high and low ednges in low space
      double log_val_lo = std::log10(energies[0]);
      double log_val_hi = std::log10(energies[nee-1]);

      double log_step_lo = log_val_lo - std::log10(energies[1]);  // This will be negative
      double log_step_hi = log_val_hi - std::log10(energies[nee-2]);
      
      for ( size_t i(1); i <= edisp_bins; i++ ) {
	log_val_lo += log_step_lo;
	log_val_hi += log_step_hi;
	out_vect[edisp_bins-i] = std::pow(10., log_val_lo);
	out_vect[nee+edisp_bins+i-1] = std::pow(10., log_val_hi);
      }

      // Swap the out_vect and the input vector
      energies.swap(out_vect);
    }

    void log_energy_ratios(const std::vector<double>& energies,
			   std::vector<double>& log_ratios) {
      log_ratios.resize(energies.size() -1);
      for ( size_t i(0); i < log_ratios.size(); i++ ) {
	log_ratios[i] = std::log(energies[i+1] / energies[i]);
      }
    }

    void svd_inverse(const CLHEP::HepMatrix& u,
		     const CLHEP::HepMatrix& v,
		     const CLHEP::HepVector& s,
		     CLHEP::HepMatrix& inverse) {

      CLHEP::HepMatrix w(s.num_row(),s.num_row());

      for ( size_t i(0); i < s.num_row(); i++ ) {
	if(s[i] <= 0) {
	  w[i][i] = 0;
	} else {
	  w[i][i] = 1./s[i];
	}
      }
      inverse = v*w*u.T();
    }

    int svd_solve(const CLHEP::HepSymMatrix& hessian,
		  const CLHEP::HepVector& gradient,
		  CLHEP::HepMatrix& u,
		  CLHEP::HepMatrix& v,
		  CLHEP::HepVector& s,
		  CLHEP::HepVector& delta,
		  double eps) {

      CLHEP::HepMatrix h(hessian);
      const int nrow = hessian.num_row();
      const int ncol = hessian.num_col();

      delta = CLHEP::HepVector(nrow);
      
      gsl_matrix * U = Matrix_Hep_to_Gsl(h);
      gsl_vector * G = Vector_Hep_to_Gsl(gradient);
      gsl_matrix * V = gsl_matrix_alloc(nrow, ncol);
      gsl_vector * S = gsl_vector_alloc(nrow);
      gsl_vector * work = gsl_vector_alloc(nrow);
      gsl_vector * x = gsl_vector_alloc(nrow);
      
      int ierr = gsl_linalg_SV_decomp(U,V,S,work);
      if(ierr) {
	return ierr;
      }

      double threshold = 0.5*sqrt(double(nrow+ncol)+1.0)*gsl_vector_get(S,0)*eps;
      for ( size_t i(0); i < nrow; i++ ) {
	if( fabs(gsl_vector_get(S,i)) < threshold ) {
	  gsl_vector_set(S,i,0);
	}
      }

      gsl_linalg_SV_solve(U,V,S,G,x);

      Vector_Gsl_to_Hep(x,delta);

      Matrix_Gsl_to_Hep(U,u);
      Matrix_Gsl_to_Hep(V,v);
      Vector_Gsl_to_Hep(S,s);

      gsl_matrix_free(U);
      gsl_vector_free(G);
      gsl_matrix_free(V);
      gsl_vector_free(S);
      gsl_vector_free(work);
      gsl_vector_free(x);

      return 0;
    }
    
    gsl_vector * Vector_Hep_to_Gsl(const CLHEP::HepVector& hep) {

      const int nrow = hep.num_row();
      gsl_vector * v = gsl_vector_alloc(nrow);
      for ( size_t i(0); i < nrow; i++ ) {
	gsl_vector_set(v,i,hep[i]);
      }
      return v;
    }

    gsl_matrix * Matrix_Hep_to_Gsl(const CLHEP::HepMatrix& hep) {

      const int nrow = hep.num_row();
      const int ncol = hep.num_col();
      gsl_matrix * m = gsl_matrix_alloc(nrow,ncol);
      for ( size_t i(0); i < nrow; i++ ) {
	for ( size_t j(0); j < ncol; j++ ) {
	  gsl_matrix_set(m,i,j,hep[i][j]);
	}
      }
      return m;
    }

    void Vector_Gsl_to_Hep(const gsl_vector * gsl,
			   CLHEP::HepVector& hep) {

      const int nrow = gsl->size;
      hep = CLHEP::HepVector(nrow);
      for ( size_t i(0); i < nrow; i++ ) {
	hep[i] = gsl_vector_get(gsl,i);
      }
    }

    void Matrix_Gsl_to_Hep(const gsl_matrix * gsl,
			   CLHEP::HepMatrix& hep) {

      const int nrow = gsl->size1;
      const int ncol = gsl->size2;
      hep = CLHEP::HepMatrix(nrow,ncol);
      for ( size_t i(0); i < nrow; i++ ) {
	for ( size_t j(0); j < ncol; j++ ) {
	  hep[i][j] = gsl_matrix_get(gsl,i,j);
	}
      }
    }

    void Vector_Hep_to_Stl(const CLHEP::HepVector& hep,
			   std::vector<float>& stl) {
      stl.resize(hep.num_row());
      for ( size_t i(0); i < stl.size(); i++ ) {
	stl[i] = hep[i];
      }
    }

    void Matrix_Hep_to_Stl(const CLHEP::HepSymMatrix& hep,
			   std::vector<float>& stl) {
      stl.resize(hep.num_row()*hep.num_col());
      size_t idx(0);
      for ( size_t i(0); i < hep.num_row(); i++ ) {
	for ( size_t j(0); j < hep.num_col(); j++ ) {
	  stl[idx] = hep[i][j];
	  idx += 1;
	}
      }
    }
    
    void Vector_Stl_to_Hep(const std::vector<float>& stl,
			   CLHEP::HepVector& hep) {
      hep = CLHEP::HepVector(stl.size());
      for ( size_t i(0); i < stl.size(); i++ ) {
	hep[i] = stl[i];
      }
    }

    void Matrix_Stl_to_Hep(const std::vector<float>& stl,
			   CLHEP::HepSymMatrix& hep) {
      size_t nrow = ::round(std::sqrt(double(stl.size())));
      if(nrow*nrow < stl.size()) nrow += 1;
      hep = CLHEP::HepSymMatrix(nrow);
      size_t idx(0);
      for ( size_t i(0); i < hep.num_row(); i++ ) {
	for ( size_t j(0); j < hep.num_col(); j++ ) {
	  hep[i][j] = stl[idx];
	  idx += 1;
	}
      } 
    }

    void sumVector(std::vector<float>::const_iterator start,
		   std::vector<float>::const_iterator stop,
		   float& value) {
      // Use a double to store the running sum, since we might be 
      // adding together a lot of things
      double sum(0.);
      for ( std::vector<float>::const_iterator itr = start;
	    itr != stop; itr++ ) {
	sum += *itr;
      }
      value = sum;
    }

    void setVectorValue(const float& val, 
			std::vector<float>::iterator start,
			std::vector<float>::iterator stop ){
      for ( std::vector<float>::iterator itr = start;
	    itr != stop; itr++ ) {
	*itr = val;
      }
    }

    void vectorAdd(std::vector<float>::const_iterator start1, 
		   std::vector<float>::const_iterator stop1, 
		   std::vector<float>::const_iterator start2, 
		   std::vector<float>::const_iterator stop2, 
		   std::vector<float>::iterator out_start,
		   std::vector<float>::iterator out_stop,
		   float fact1, float fact2) {
      size_t s1 = stop1 - start1;
      size_t s2 = stop2 - start2;
      size_t s3 = out_stop - out_start;

      if ( s1 != s2 ||
	   s1 != s3 ) {
	throw std::runtime_error("Size of vectors in FitUtils::vectorAdd does not match.");
	return;
      }
      std::vector<float>::const_iterator itr1 = start1;
      std::vector<float>::const_iterator itr2 = start2;
      std::vector<float>::iterator itr_out = out_start;

      for ( ; itr1 != stop1; itr1++, itr2++, itr_out++ ) {
	*itr_out = ((fact1*(*itr1)) + (fact2*(*itr2)));
      }

   }

    void vectorMultiply(std::vector<float>::const_iterator start1, 
			std::vector<float>::const_iterator stop1, 
			std::vector<float>::const_iterator start2, 
			std::vector<float>::const_iterator stop2, 
			std::vector<float>::iterator out_start,
			std::vector<float>::iterator out_stop) {
      size_t s1 = stop1 - start1;
      size_t s2 = stop2 - start2;
      size_t s3 = out_stop - out_start;

      if ( s1 != s2 ||
	   s1 != s3 ) {
	throw std::runtime_error("Size of vectors in FitUtils::vectorMultiply does not match.");
	return;
      }
      std::vector<float>::const_iterator itr1 = start1;
      std::vector<float>::const_iterator itr2 = start2;
      std::vector<float>::iterator itr_out = out_start;

      for ( ; itr1 != stop1; itr1++, itr2++, itr_out++  ) {
	*itr_out = ((*itr1) * (*itr2));
      }
    }

    void vectorSubtract(std::vector<float>::const_iterator start1, 
			std::vector<float>::const_iterator stop1, 
			std::vector<float>::const_iterator start2, 
			std::vector<float>::const_iterator stop2, 
			std::vector<float>::iterator out_start,
			std::vector<float>::iterator out_stop) {
      size_t s1 = stop1 - start1;
      size_t s2 = stop2 - start2;
      size_t s3 = out_stop - out_start;

      if ( s1 != s2 ||
	   s1 != s3 ) {
	throw std::runtime_error("Size of vectors in FitUtils::vectorSubtract does not match.");
	return;
      }
      std::vector<float>::const_iterator itr1 = start1;
      std::vector<float>::const_iterator itr2 = start2;
      std::vector<float>::iterator itr_out = out_start;

      for ( ; itr1 != stop1; itr1++, itr2++, itr_out++ ) {
	*itr_out = (*itr1) - (*itr2);
      }
    }
   
    void multiplyByScalar(std::vector<float>::iterator vstart,
			 std::vector<float>::iterator vstop,
			 double scalar) {
      for ( std::vector<float>::iterator itr = vstart ; itr != vstop; itr++ ) {
	*itr *= scalar;
      }
    }    

    void innerProduct(std::vector<float>::const_iterator start1, 
		      std::vector<float>::const_iterator stop1, 
		      std::vector<float>::const_iterator start2, 
		      std::vector<float>::const_iterator stop2, 
		      float& value) {
      size_t s1 = stop1 - start1;
      size_t s2 = stop2 - start2;
      if ( s1 != s2 ) {
	throw std::runtime_error("Size of vectors in FitUtils::innerProduct does not match.");
	return;
      }
      // use a double to store the sum, snice we might be adding a lot of things
      double sum(0.);
      std::vector<float>::const_iterator itr1 = start1;
      std::vector<float>::const_iterator itr2 = start2;
      for ( ; itr1 != stop1; itr1++, itr2++ ) {
	sum += ((*itr1) * (*itr2));
      }
      value = sum;
   }
    
    void innerProduct(std::vector<float>::const_iterator start1, 
		      std::vector<float>::const_iterator stop1, 
		      std::vector<float>::const_iterator start2, 
		      std::vector<float>::const_iterator stop2, 
		      std::vector<float>::const_iterator start3, 
		      std::vector<float>::const_iterator stop3, 
		      float& value) {
      size_t s1 = stop1 - start1;
      size_t s2 = stop2 - start2;
      size_t s3 = stop3 - start3;
      if ( s1 != s2 ||
	   s1 != s3 ) {
	throw std::runtime_error("Size of vectors in FitUtils::innerProduct does not match.");
	return;
      }   
      // use a double to store the sum, snice we might be adding a lot of things
      double sum(0);
      std::vector<float>::const_iterator itr1 = start1;
      std::vector<float>::const_iterator itr2 = start2;
      std::vector<float>::const_iterator itr3 = start3;
      for ( ; itr1 != stop1; itr1++, itr2++, itr3++ ) {
	sum += ((*itr1) * (*itr2) * (*itr3));
      }
      value = sum;
    }

    void innerProduct(std::vector<float>::const_iterator start1, 
		      std::vector<float>::const_iterator stop1, 
		      std::vector<float>::const_iterator start2, 
		      std::vector<float>::const_iterator stop2, 
		      std::vector<float>::const_iterator start3, 
		      std::vector<float>::const_iterator stop3, 
		      std::vector<float>::const_iterator start4, 
		      std::vector<float>::const_iterator stop4, 
		      float& value) {
      size_t s1 = stop1 - start1;
      size_t s2 = stop2 - start2;
      size_t s3 = stop3 - start3;
      size_t s4 = stop4 - start4;
      
      if ( s1 != s2 ||
	   s1 != s3 || 
	   s1 != s4 ) {
	throw std::runtime_error("Size of vectors in FitUtils::innerProduct does not match.");
	return;
      }   
      // use a double to store the sum, snice we might be adding a lot of things
      double sum(0);
      std::vector<float>::const_iterator itr1 = start1;
      std::vector<float>::const_iterator itr2 = start2;
      std::vector<float>::const_iterator itr3 = start3;
      std::vector<float>::const_iterator itr4 = start4;
      for ( ; itr1 != stop1; itr1++, itr2++, itr3++, itr4++ ) {
	sum += ((*itr1) * (*itr2) * (*itr3) * (*itr4));
      }
      value = sum;
    }

    void fracDiff(std::vector<float>::const_iterator data_start, 
		  std::vector<float>::const_iterator data_stop,
		  std::vector<float>::const_iterator model_start,
		  std::vector<float>::const_iterator model_stop,
		  std::vector<float>::iterator out_start,
		  std::vector<float>::iterator out_stop) {
      size_t s1 = data_stop - data_start;
      size_t s2 = model_stop - model_start;
      size_t s3 = out_stop - out_start;
      if ( s1 != s2 ||
	   s1 != s3 ) {
	throw std::runtime_error("Size of vectors in FitUtils::fracDiff does not match.");
	return;
      }
      std::vector<float>::const_iterator itr_data = data_start;
      std::vector<float>::const_iterator itr_model = model_start;
      std::vector<float>::iterator itr_out = out_start;
      for ( ; itr_data != data_stop; itr_data++, itr_model++, itr_out++ ) {
	if ( *itr_data <= 0 ) {
	  *itr_out = 1.;
	  continue;
	}
	if ( *itr_model <= 0. ) {
	  throw std::runtime_error("Negative model counts in FitUtils::fracDiff.");
	  return;
	}
	*itr_out = ( 1. - ( *itr_data / *itr_model ) );
      }
    }

    void data_over_model2(std::vector<float>::const_iterator data_start, 
			  std::vector<float>::const_iterator data_stop,
			  std::vector<float>::const_iterator model_start,
			  std::vector<float>::const_iterator model_stop,
			  std::vector<float>::iterator out_start,
			  std::vector<float>::iterator out_stop) {
      size_t s1 = data_stop - data_start;
      size_t s2 = model_stop - model_start;
      size_t s3 = out_stop - out_start;
      if ( s1 != s2 ||
	   s1 != s3 ) {
	throw std::runtime_error("Size of vectors in FitUtils::data_over_model2 does not match.");
	return;
      }
      std::vector<float>::const_iterator itr_data = data_start;
      std::vector<float>::const_iterator itr_model = model_start;
      std::vector<float>::iterator itr_out = out_start;
      for ( ; itr_data != data_stop; itr_data++, itr_model++, itr_out++ ) {
	if ( *itr_data <= 0 ) {
	  *itr_out = 0.;
	  continue;
	}
	if ( *itr_model <= 0. ) {
	  throw std::runtime_error("Negative model counts in FitUtils::data_over_model2.");
	  return;
	}
	*itr_out = *itr_data / ( (*itr_model) * (*itr_model) );
      }
    }    
   
    void reshapeVector(const std::vector<float>& invect,
		       size_t n_i, size_t n_j,
		       std::vector<std::vector<float> >& outvect) {
      if ( invect.size() != (n_i*n_j) ) {
	throw std::runtime_error("ize of vector in FitUtils::reshapeVector does not match output shape.");
	return;
      }
      // setup the output
      outvect.resize(n_i);
      size_t idx(0);
      for ( size_t i(0); i < n_i; i++ ) {
	outvect[i].resize(n_j);
	for ( size_t j(0); j < n_j; j++, idx++ ) {
	  outvect[i][j] = invect[idx];
	}
      }
    }

    double symmetricError(double pos_err, double neg_err) {
      return neg_err < 0 ? pos_err : (pos_err + neg_err) / 2.;
    }

    void sumModel_Init(const std::vector<std::vector<float> >& templates,
		       const std::vector<float>& fixed,
		       std::vector<float>& total) {
 
      std::copy( fixed.begin(), fixed.end(), total.begin() );
      for ( size_t i(0); i < templates.size(); i++ ) {
	const std::vector<float>& tmpl = templates[i];
	std::vector<float>::const_iterator tmpl_itr = tmpl.begin();
	std::vector<float>::iterator outItr = total.begin();
	for ( std::vector<float>::iterator outItr = total.begin(); outItr != total.end(); outItr++, tmpl_itr++ ) {
	  *outItr += *tmpl_itr;
	}
      } 
    }    
    
    void sumModel(const CLHEP::HepVector& norms,
		  const std::vector<const std::vector<float>* >& templates,
		  const std::vector<float>& fixed,
		  std::vector<float>& total,
		  size_t firstBin,
		  size_t lastBin) {      
      
      size_t npar = norms.num_row();
      size_t nbin = total.size();
      if ( fixed.size() != nbin ) {
	throw std::runtime_error("Size of fixed component does not match size of total model in FitUtils::sumModel");
	return;
      }
      if ( templates.size() != npar ) {
	throw std::runtime_error("Number of templates does not match number of free parameters in FitUtils::sumModel");
	return;
      }
      // copy over the fixed components
      std::vector<float>::const_iterator copy_start = fixed.begin() + firstBin;
      std::vector<float>::const_iterator copy_stop = lastBin == 0 ? fixed.end() : fixed.begin() + lastBin;
      std::vector<float>::iterator out_start = total.begin() + firstBin;
      std::vector<float>::iterator out_stop = lastBin == 0 ? total.end() : total.begin() + lastBin;
      std::vector<float>::iterator outItr = out_start;

      std::copy( copy_start, copy_stop, outItr );
      for ( size_t i(0); i < npar; i++ ) {
	float norm = norms[i];
	const std::vector<float>& tmpl = *(templates[i]);
	std::vector<float>::const_iterator tmpl_itr = tmpl.begin() + firstBin;
	for ( outItr = out_start; outItr != out_stop; outItr++, tmpl_itr++ ) {
	  *outItr += norm * ( *tmpl_itr );
	}
      }
    }
    
    double getLogLikelihood(const std::vector<float>& data,
			    const CLHEP::HepVector& norms,
			    const std::vector<const std::vector<float>* >& templates,
			    const std::vector<float>& fixed,
			    const FitScanMVPrior* prior,
			    const std::vector<float>* weights,   
			    std::vector<float>& model,
			    size_t firstBin,
			    size_t lastBin,
			    int verbose) {

      std::vector<float>::const_iterator data_start = data.begin() + firstBin;
      std::vector<float>::const_iterator data_end = lastBin == 0 ? data.end() : data.begin() + lastBin;
      std::vector<float>::const_iterator model_start = model.begin() + firstBin;
      std::vector<float>::const_iterator model_end = lastBin == 0 ? model.end() : model.begin() + lastBin;

      sumModel(norms,templates,fixed,model,firstBin,lastBin);
      double logLikeVal = weights == 0 ? 
	logLikePoisson(data_start,data_end,model_start,model_end) :
	logLikePoisson(data_start,data_end,model_start,model_end, 
		       weights->begin() + firstBin,
		       lastBin == 0 ? weights->end() :  weights->begin() + lastBin);
      if ( prior ) {
	double logLikePrior(0.);
	prior->logLikelihood(norms,logLikePrior);
	logLikeVal += logLikePrior;
      }
      return logLikeVal;
    }

    void getGradientAndHessian(const std::vector<float>& data,
			       const CLHEP::HepVector& norms,
			       const std::vector<const std::vector<float>* >& templates,
			       const std::vector<float>& fixed,
			       const FitScanMVPrior* prior,			       
			       const std::vector<float>* weights,
			       std::vector<float>& model,
			       CLHEP::HepVector& gradient,
			       CLHEP::HepSymMatrix& hessian,
			       size_t firstBin,
			       size_t lastBin,
			       int verbose) {      

      size_t start = firstBin;
      size_t stop = lastBin == 0 ? data.size() : lastBin;

      size_t sz = stop - start;
      std::vector<float> fdiff(sz);
      std::vector<float> w2(sz);
      
      sumModel(norms,templates,fixed,model,firstBin,lastBin);

      fracDiff(data.begin() + start, data.begin() + stop,
	       model.begin() + start, model.begin() + stop,
	       fdiff.begin(),fdiff.end());
	       
      data_over_model2(data.begin() + start, data.begin() + stop,
		       model.begin() + start, model.begin() + stop,
		       w2.begin(),w2.end());

      size_t npar = norms.num_row();
      for ( size_t i(0); i < npar; i++ ) {
	const std::vector<float>& tmpl_i = *(templates[i]);
	float gval(0.);
	// this is g_i = Sum fdiff * templates[i]
	if ( weights != 0 ) {
	  innerProduct( weights->begin() + start, weights->begin() + stop,
		        fdiff.begin(), fdiff.end(),
			tmpl_i.begin() + start, tmpl_i.begin() + stop, 
			gval );
	} else {
	  innerProduct( fdiff.begin(), fdiff.end(),
			tmpl_i.begin() + start, tmpl_i.begin() + stop, 
			gval );
	}
	gradient[i] = gval;
	// note that the inner loop starts at i.
	// we only need elements at or above the diagonal
	for ( size_t j(i); j < npar; j++ ) {
	  const std::vector<float>& tmpl_j = *(templates[j]);
	  float hval(0.);
	  // this is h_ij = Sum w2 * templates[i] * templates[j] 
	  if ( weights != 0 ) {
	    innerProduct( weights->begin() + start, weights->begin() + stop,
			  w2.begin(), w2.end(),
			  tmpl_i.begin() + start, tmpl_i.begin() + stop, 
			  tmpl_j.begin() + start, tmpl_j.begin() + stop, 
			  hval );
	  } else {
	    innerProduct( w2.begin(), w2.end(),
			  tmpl_i.begin() + start, tmpl_i.begin() + stop, 
			  tmpl_j.begin() + start, tmpl_j.begin() + stop, 
			  hval );
	  }
	  hessian[i][j] = hval;
	}
      }
      
      // Add in the terms from the prior, if requested
      if ( prior ) {
	CLHEP::HepVector gPrior;
	prior->gradient(norms,gPrior);
	gradient += gPrior;
	hessian += prior->hessian();
      }
      
      return;
    }

    int getDeltaAndCovarAndEDM(const CLHEP::HepSymMatrix& hessian,
			       const CLHEP::HepVector& gradient,
			       const CLHEP::HepVector& norms,
			       CLHEP::HepSymMatrix& covar,
			       CLHEP::HepVector& delta,
			       double& edm) {      
      covar.assign(hessian);
      int ifail(0);
      covar.invert(ifail);
      if ( ifail ) {
	return ifail;
      }
      
      delta = covar * gradient;   
      int npar = gradient.num_row();
      edm = 0.;

      // Protect against components going negative.
      // The retry=true block below recalculates
      // the delta vector after setting the 
      // offending gradient element to zero.
      // That is important if there are other components
      // that are highly anti-correlated with the
      // component that is going negative.
      int retry(0);
      CLHEP::HepVector g2(gradient);
      for ( int i(0); i < npar; i++ ) {
	if ( delta[i] > norms[i] ) {
	  delta[i] = 0.9999*norms[i];
	  g2[i] = 0.;
	  retry += 1;
	}
	edm += std::fabs(delta[i] * g2[i]);
      }      
      
      if ( retry > 0 ) {
        delta = covar*g2;
        for ( int i(0); i < npar; i++ ) {
          if ( delta[i] > norms[i] ) {
            delta[i] = 0.9999*norms[i];
            g2[i] = 0.;
	    delta[i] = 0.;
            retry += 1;
          }
        }
      }

      return 0;
    }

    int getDeltaAndCovarAndEDM(const CLHEP::HepSymMatrix& hessian,
			       const CLHEP::HepVector& gradient,
			       const CLHEP::HepVector& norms,
			       CLHEP::HepSymMatrix& covar,
			       CLHEP::HepVector& delta,
			       double& edm,
			       double lambda) {      

      const int npar = gradient.num_row();
      CLHEP::HepSymMatrix alpha = hessian;
      for( int i(0); i < npar; i++)
	alpha[i][i] *= (1+lambda);

      CLHEP::HepMatrix u, v;
      CLHEP::HepVector s;
      int ifail = svd_solve(alpha, gradient, u, v, s, delta);
      if(ifail) {
	return ifail;
      }

      covar.assign(hessian);
      covar.invert(ifail);
      for ( int i(0); i < npar; i++ ) {
	if(covar[i][i] < 0)
	  ifail = 1;
      }

      if ( ifail ) {
	CLHEP::HepMatrix inverse;
	svd_inverse(u,v,s,inverse);
	covar.assign(inverse);
      }

      const double lo_bound = 0;
      int retry(0);
      CLHEP::HepVector g2(gradient);
      for ( int i(0); i < npar; i++ ) {

	if ( norms[i] - delta[i] < lo_bound ) {
	  delta[i] = std::min(0.9999*norms[i], norms[i]-lo_bound);
	  g2[i] = 0.;
	  retry += 1;
	}
      }
      
      edm = 0.;
      for ( int i(0); i < npar; i++ ) {
	edm += std::fabs(delta[i] * gradient[i]);
      }

      if ( retry > 0 ) {

	int ifail = svd_solve(alpha, g2, u, v, s, delta);
	if(ifail) {
	  return ifail;
	}

	for ( int i(0); i < npar; i++ ) {
	  if ( norms[i] - delta[i] < lo_bound ) {
	    delta[i] = std::min(0.9999*norms[i], norms[i]-lo_bound);
	    g2[i] = 0.;
	    retry += 1;
	  }
	}
      }

      return 0;
    }


    int getDeltaAndCovarAndEDM_STL(const std::vector<float>& hessian,
				   const std::vector<float>& gradient,
				   const std::vector<float>& norms,
				   std::vector<float>& covar,
				   std::vector<float>& delta,
				   std::vector<double>& edm,
				   double lambda) {
      CLHEP::HepSymMatrix hessian_hep;
      CLHEP::HepVector gradient_hep;
      CLHEP::HepVector norms_hep;
      FitUtils::Matrix_Stl_to_Hep(hessian, hessian_hep);
      FitUtils::Vector_Stl_to_Hep(gradient, gradient_hep);
      FitUtils::Vector_Stl_to_Hep(norms, norms_hep);

      CLHEP::HepSymMatrix covar_hep;
      CLHEP::HepVector delta_hep;

      int status(0);
      if ( lambda > 0. ) {
	status = getDeltaAndCovarAndEDM(hessian_hep, gradient_hep, norms_hep,
					covar_hep, delta_hep, edm[0], lambda);
      } else {
	status = getDeltaAndCovarAndEDM(hessian_hep, gradient_hep, norms_hep,
					covar_hep, delta_hep, edm[0]);
      }
      FitUtils::Matrix_Hep_to_Stl(covar_hep, covar);
      FitUtils::Vector_Hep_to_Stl(delta_hep, delta);
      return status;
    }
    


    int fitNorms_newton(const std::vector<float>& data,
			const CLHEP::HepVector& initNorms,
			const std::vector<const std::vector<float>* >& templates,
			const std::vector<float>& fixed,
			const FitScanMVPrior* prior,
			const std::vector<float>* weights,
			double tol, int maxIter, double initLambda,
			CLHEP::HepVector& norms,
			CLHEP::HepSymMatrix& covar,
			CLHEP::HepVector& gradient,
			std::vector<float>& model,
			double& edm,
			double& logLikeVal,
			size_t firstBin, 
			size_t lastBin,
			int verbose) {


      // copy over the initial parameters
      norms = initNorms;
      
      // local stuff
      size_t npar = initNorms.num_row();
      double logLikeInit(0.);
      double logLikeIter(0.);
      double deltaLogLike(0.);

      std::vector<float>::const_iterator data_start = data.begin() + firstBin;
      std::vector<float>::const_iterator data_end = lastBin == 0 ? data.end() : data.begin() + lastBin;
      std::vector<float>::const_iterator model_start = model.begin() + firstBin;
      std::vector<float>::const_iterator model_end = lastBin == 0 ? model.end() : model.begin() + lastBin;
 
      if ( verbose > 0 ) {
	std::cout << "Init Newton's method fit.  NPar = " <<  npar << std::endl;
	printVector("Init Pars: ",initNorms);
      }
            
      logLikeVal = getLogLikelihood(data,norms,templates,fixed,prior,weights,
				    model,firstBin,lastBin,verbose);
      logLikeInit = logLikeVal;

      if ( npar == 0 ) {
	std::cout << "warning, no free parameters " << std::endl;
	return 0;
      }
      
      gradient = CLHEP::HepVector(npar);
      CLHEP::HepVector delta(npar);
      CLHEP::HepSymMatrix hessian(npar);      

      // Set the EDM larger than the tolerance
      edm = 100*tol;
      int iter(0);
      double lambda = initLambda;
      bool stepOk = true;

      // Check to see if there are any counts.  
      float data_total(0.);
      sumVector(data_start,data_end,data_total);
      if ( data_total < 1e-9 ) {
	// no counts, set some defaults and bail out
	norms[npar-1] = 0.;
	// the gradient can be useful
	getGradientAndHessian(data,norms,templates,fixed,prior,weights,model,
			      gradient,hessian,firstBin,lastBin,verbose);
	covar = CLHEP::HepSymMatrix(npar);
	// this will skip the loop below
	edm = 0.;
	// make sure we don't fail the sanity check below
	logLikeInit = -1.0e99;
      }
      
      // loop until convergence or max iterations
      while ( iter < maxIter && std::fabs(edm) > tol ) {
	
	// do the derivative stuff
	getGradientAndHessian(data,norms,templates,fixed,prior,weights,model,
			      gradient,hessian,firstBin,lastBin,verbose);
	if ( verbose > 2 ) {
	  printMatrix("Hesse: ",hessian);
	  printVector("Grad: ",gradient);
	  if ( prior ) {
	    CLHEP::HepVector gPrior;
	    prior->gradient(norms,gPrior);
	    printVector("Prior Grad: ",gPrior);
	    printMatrix("Prior Hesse: ",prior->hessian());
	  }
	}
	// invert the hessian and get the delta for the next iteration
	// this also compute the estimated distance to the minimum
	int covOk = 0;

	if(lambda > 0)
	  covOk = getDeltaAndCovarAndEDM(hessian,gradient,norms,
					 covar,delta,edm,lambda);
	else
	  covOk = getDeltaAndCovarAndEDM(hessian,gradient,norms,
					 covar,delta,edm);

	if ( verbose > 2 ) {
	  printMatrix("Cov: ",covar);
	}

	if ( covOk !=0 ) {
	  // non-zero value means it failed to invert the cov. matrix
	  if ( verbose > 0 ) {
	    std::cout << "Matrix inversion failed with code " << covOk << std::endl;
	  }
	  return covOk;
	}

	if(lambda > 0) {

	  CLHEP::HepVector norms_tmp = norms;
	  norms_tmp -= delta;

	  // Check whether the fit has improved
	  logLikeIter = getLogLikelihood(data,norms_tmp,templates,fixed,prior,weights,
					 model,firstBin,lastBin,verbose);
	  deltaLogLike = logLikeIter-logLikeVal;

	  if (logLikeIter > logLikeVal) {
	    norms -= delta;
	    lambda = std::max(1E-8,lambda/5.);
	    logLikeVal = logLikeIter;
	    stepOk = true;
	  } else {
	    lambda *= 3;
	    stepOk = false;
	  }

	} else {
	  // update the parameter values
	  norms -= delta;
	}

	// catch fits that are diverging
	for ( size_t iCheck(0); iCheck < norms.num_row(); iCheck++ ) {
	  if ( norms[iCheck] > 1e9 ) {
	    if ( verbose > 0 ) {
	      std::cout << "A fit using Newton's method seems to have diverged." << std::endl;
	      printVector("Norms:",norms);
	    }
	    return -4;
	  }
	}

	if ( verbose > 1 ) {
	  std::cout << iter << ' ' << std::endl;
	  printVector("Delta:",delta);
	  printVector("Norms:",norms);
	  std::cout << "EDM " << edm << std::endl;
	  std::cout << "lambda " << lambda << std::endl;
	  std::cout << "logLike " << logLikeIter << " " << deltaLogLike
		    << std::endl;
	}

	// Break if tolerance criteria are met
	if((lambda <= 0 && std::fabs(edm) < tol) || 
	   (lambda > 0 && std::fabs(edm) < tol && 
	    deltaLogLike < tol && stepOk)) {
	  break;
	}

	iter++;
      }
      
      // get the final log-likelihood
      if(lambda == 0) {
	logLikeVal = getLogLikelihood(data,norms,templates,fixed,prior,weights,model,
				      firstBin,lastBin,verbose);
      } else {
	getGradientAndHessian(data,norms,templates,fixed,prior,weights,model,
			      gradient,hessian,firstBin,lastBin,verbose);
	getDeltaAndCovarAndEDM(hessian,gradient,norms,
			       covar,delta,edm,0.0);
      }

      // check to see if we reach the max iterations
      if ( iter >= maxIter ) {
	if ( verbose > 0 ) {
	  std::cout << "Reached max iterations." << std::endl;
	  printVector("Init:",initNorms);
	  printVector("Norms:",norms);
	  std::cout << "Log-likelihood Init " << logLikeInit << std::endl;
	  std::cout << "Log-likelihood      " << logLikeVal << std::endl 
		    << std::endl;	  
	}
	return -2;
      }

      if ( verbose > 0 ) {
	std::cout << "Converged, niter = " << iter << std::endl;	
	printVector("Init:",initNorms);
	printVector("Norms:",norms);
	std::cout << "Log-likelihood Init " << logLikeInit << std::endl;
	std::cout << "Log-likelihood      " << logLikeVal << std::endl 
		  << std::endl;
      }

      // Sanity check
      if (  logLikeVal - logLikeInit < -1. ) {
	// The fit got worse, reset things
	if ( verbose > 0 ) {
	  std::cout << "Fit seems to have made things worse " << logLikeVal << ' ' << logLikeInit << std::endl
		    << "  This is usually caused by degeneracy between the test source and a poorly constrained source." << std::endl;
	}
	logLikeVal = logLikeInit;
	norms = CLHEP::HepVector(npar);
	sumModel(norms,templates,fixed,model,firstBin,lastBin);
	edm = 0.;
	return -8;
      }

      return 0;
    }

    
    int fitLogNorms_newton(const std::vector<float>& data,
			   const CLHEP::HepVector& initNorms,
			   const std::vector<const std::vector<float>* >& templates,
			   const std::vector<float>& fixed,
			   const FitScanMVPrior* prior,
			   const std::vector<float>* weights,
			   double tol, int maxIter, double initLambda,
			   CLHEP::HepVector& norms,
			   CLHEP::HepSymMatrix& covar,
			   CLHEP::HepVector& gradient,
			   std::vector<float>& model,
			   double& edm,
			   double& logLikeVal,
			   size_t firstBin, 
			   size_t lastBin,
			   int verbose) {

      logLikeVal = 0.;

      // local stuff
      size_t npar = initNorms.num_row();

      if ( npar == 0 ) {
	std::cout << "warning, no free parameters " << std::endl;
	norms = initNorms;
	sumModel(norms,templates,fixed,model,firstBin,lastBin);
      
	std::vector<float>::const_iterator data_start = data.begin() + firstBin;
	std::vector<float>::const_iterator data_end = lastBin == 0 ? data.end() : data.begin() + lastBin;
	std::vector<float>::const_iterator model_start = model.begin() + firstBin;
	std::vector<float>::const_iterator model_end = lastBin == 0 ? model.end() : model.begin() + lastBin;
	
	logLikeVal = weights == 0 ? 
	  logLikePoisson(data_start,data_end,model_start,model_end) :
	  logLikePoisson(data_start,data_end,model_start,model_end,
			 weights->begin() + firstBin,
			 lastBin == 0 ? weights->end() :  weights->begin() + lastBin);
	

	if ( prior ) {
	  double logLikePrior(0.);
	  prior->logLikelihood(norms,logLikePrior);
	  logLikeVal += logLikePrior;
	}
	
	if ( verbose > 0 ) {
	  std::cout << "Log-likelihood " << logLikeVal << std::endl << std::endl;
	}

	return 0;
      }

      if ( verbose > 0 ) {
	std::cout << "Init Newton's method fit.  NPar = " <<  npar << std::endl;
	printVector("Init Pars: ",initNorms);
      }

      // copy over the initial parameters      
      norms = initNorms;
      // fix the last parameter ( we don't want zeros)
      if ( norms[npar-1] < 1e-6 ) {
	norms[npar-1] = 1e-6;
      }

      CLHEP::HepVector logNorms(npar);
      
      gradient = CLHEP::HepVector(npar);
      CLHEP::HepVector delta(npar);
      CLHEP::HepSymMatrix hessian(npar);      


      // Set the EDM larger than the tolerance
      edm = 100*tol;
      int iter(0);
      
      // loop until convergence or max iterations
      while ( std::fabs(edm) > tol &&
	      iter < maxIter ) {
	
	// do the derivative stuff
	getGradientAndHessian(data,norms,templates,fixed,prior,weights,model,gradient,hessian,firstBin,lastBin,verbose);
	
	// add in the extra factors from the chain rule
	for ( size_t iChain(0); iChain < npar; iChain++ ) {
	  gradient[iChain] *= norms[iChain];
	  logNorms[iChain] = log(norms[iChain]);
	  for ( size_t jChain(iChain); jChain < npar; jChain++ ) {
	    hessian[iChain][jChain] *= (norms[iChain]*norms[jChain]);
	  }	  
	}
	
	if ( verbose > 2 ) {
	  printVector("Log Norms: ",logNorms);
	  printMatrix("Hesse: ",hessian);
	  printVector("Grad: ",gradient);
	  if ( prior ) {
	    CLHEP::HepVector gPrior;
	    prior->gradient(norms,gPrior);
	    printVector("Prior Grad: ",gPrior);
	    printMatrix("Prior Hesse: ",prior->hessian());
	  }
	}

	// invert the hessian and get the delta for the next iteration
	// this also compute the estimated distance to the minimum
	int covOk = getDeltaAndCovarAndEDM(hessian,gradient,norms,covar,delta,edm);
	if ( verbose > 2 ) {
	  printMatrix("Cov: ",covar);
	}
		
	if ( covOk !=0 ) {
	  // non-zero value means it failed to invert the cov. matrix
	  if ( verbose > 0 ) {
	    std::cout << "Matrix inversion failed." << std::endl;
	  }
	  return covOk;
	}
	
	// update the parameter values
	logNorms -= delta;

	if ( verbose > 1 ) {
	  std::cout << iter << ' ' << std::endl;
	  printVector("Delta:",delta);
	  printVector("LogNorms:",logNorms);
	  std::cout << "EDM " << edm << std::endl;
	}

	// catch fits that are diverging
	for ( size_t iCheck(0); iCheck < norms.num_row(); iCheck++ ) {
	  if ( std::fabs(logNorms[iCheck]) > 20 ) {
	    if ( verbose > 0 ) {
	      std::cout << "A fit using Newton's method seems to have diverged." << std::endl;
	      printVector("LogNorms:",logNorms);
	    }
	    return -4;
	  }
	  norms[iCheck] = exp(logNorms[iCheck]);
	  // here we only loop up to iCheck, this makes sure we are using 
	  // the up-to-date values of iNorms
	  for ( size_t jCheck(0); jCheck < iCheck; jCheck ++ ) {
	    covar[iCheck][jCheck] /= ( norms[iCheck] * norms[jCheck] );
	  }
	}
	iter++;
      }
      // check to see if we reach the max iterations
      if ( iter >= maxIter ) {
	if ( verbose > 0 ) {
	  std::cout << "Reached max iterations." << std::endl;
	}
	return -2;
      }
      
      if ( verbose > 0 ) {
	std::cout << "Converged, niter = " << iter << std::endl;	
	printVector("Init:",initNorms);
	printVector("Norms:",norms);
      }

      // get the final log-likelihood
      sumModel(norms,templates,fixed,model,firstBin,lastBin);
      
      std::vector<float>::const_iterator data_start = data.begin() + firstBin;
      std::vector<float>::const_iterator data_end = lastBin == 0 ? data.end() : data.begin() + lastBin;
      std::vector<float>::const_iterator model_start = model.begin() + firstBin;
      std::vector<float>::const_iterator model_end = lastBin == 0 ? model.end() : model.begin() + lastBin;

      logLikeVal = weights == 0 ? 
	logLikePoisson(data_start,data_end,model_start,model_end) :
	logLikePoisson(data_start,data_end,model_start,model_end,
		       weights->begin() + firstBin,
		       lastBin == 0 ? weights->end() : weights->begin() + lastBin);
      
      if ( prior ) {
	double logLikePrior(0.);
	prior->logLikelihood(norms,logLikePrior);
	logLikeVal += logLikePrior;
      }

      if ( verbose > 0 ) {
	std::cout << "Log-likelihood " << logLikeVal << std::endl << std::endl;
      }
      return 0;
    }

    void extractSpectralVals(const Source& source,
			     const std::vector<double>& energies,
			     std::vector<double>& specVals) {

      const optimizers::Function& spec = source.spectrum();
      specVals.clear();
      for ( std::vector<double>::const_iterator itr = energies.begin();
	    itr != energies.end(); itr++ ) {
	optimizers::dArg darg(*itr);
	specVals.push_back( spec(darg) );
      }
      return;
    }

    void extractSpectralDerivs(const Source& source,
			       const std::vector<double>& energies,
			       const std::vector<std::string>& paramNames,
			       std::vector<std::vector<double> >& derivVals) {

      const optimizers::Function& spec = source.spectrum();
      derivVals.resize(paramNames.size());
      size_t idx(0);
      for ( std::vector<std::string>::const_iterator itr = paramNames.begin();
	    itr != paramNames.end(); itr++, idx++ ) {
	derivVals[idx].resize(energies.size());
	size_t idx2(0);
	for ( std::vector<double>::const_iterator itr2 = energies.begin();
	      itr2 != energies.end(); itr2++, idx2++ ) {
	  optimizers::dArg eArg(*itr2);
	  derivVals[idx][idx2] = spec.derivByParam(eArg, *itr);
	}
      }
    }


    void extractNPreds(const Source& source,
		       const std::vector<double>& energies,
		       std::vector<double>& nPreds) {
      size_t nE = energies.size();
      for ( size_t iE(0); iE < nE-1; iE++ ) {
	nPreds.push_back(  source.Npred(energies[iE],energies[iE+1]) );
      }
      return;
    }


    void extractModelFromSource(const Source& source,
				const BinnedLikelihood& logLike,
				std::vector<float>& model,
				bool rescaleToNormOne) {
      
      model.resize(logLike.data_map_size(), 0);
      setVectorValue(0., model.begin(), model.end());

      if(logLike.hasSrcNamed(source.getName())) {
	bool hasSourceMap = logLike.hasSourceMap(source.getName());
	BinnedLikelihood& nclike = const_cast<BinnedLikelihood&>(logLike);
	nclike.computeModelMap(source.getName(),model);
	if ( !hasSourceMap ) {
	  nclike.eraseSourceMap(source.getName());
	}

      } else {	

	SourceMap* theMap(0);
	
	// Otherwise, build a map ourselves
	Source* nc_source = const_cast<Source*>(&source);
	BinnedLikelihood& nclike = const_cast<BinnedLikelihood&>(logLike);
	theMap = new SourceMap(*nc_source,
			       &logLike.dataCache(),
			       logLike.observation(),
			       logLike.config(),
			       nclike.drm());
	
	nclike.updateModelMap_fromSrcMap(model, source, theMap);	
	//std::vector<double> specVals;
	//FitUtils::extractSpectralVals(source,logLike.energies(),specVals);
	//FitUtils::extractModelCounts(*theMap,logLike.energies(),specVals,model);
	
	delete theMap;
      }

      if ( rescaleToNormOne ) {
	double refValue = source.spectrum().normPar().getValue();
	double scaleFactor = 1. / refValue;
	FitUtils::multiplyByScalar(model.begin(),model.end(),scaleFactor);
      }
    }
  

    void extractModelCounts(SourceMap& sourceMap,
			    const std::vector<double>& energies,
			    const std::vector<double>& specVals,			    
			    std::vector<float>& modelCounts) {
      
      // This functionality is copied from Source::pixelCounts
      // It integtrates the spectrum for each pixel at each energy bin by
      // doing quadrature in log-log space, suggested by J. Ballet.         

      // The spectral values are defined at the bin edges, so there
      // is one less energy bin than spectral value
      size_t nebins = specVals.size() - 1;
      const std::vector<float>& model = sourceMap.model();
      // This is the stride from one energy bin to the next in the model vector
      size_t npix = model.size() / specVals.size();
      
      size_t idx(0);
      // Loop over energy bins
      for ( size_t i(0); i < nebins; i++ ) {
	// Latch spectral info
	double x1 = energies[i]*specVals[i];
	double x2 = energies[i+1]*specVals[i+1];
	double pf = std::log(energies[i+1]/energies[i]) / 2.;
	// Loop over pixels
	for ( size_t j(0); j < npix; j++, idx++ ) {
	  // Integrate the spectrum over the energy bin for this pixel
	  // using the weights from the model
	  float val = pf*(x1*model[idx] + x2*model[idx+npix]);
	  if ( val <= 0 ) {
	    std::cout << "Negative model component " << sourceMap.name() << ' ' << idx << ' ' << val << ' ' 
		      << pf << ' ' << x1 << ' ' << x2 << ' ' 
		      << model[idx] << ' ' << model[idx+npix] << std::endl;
	  }
	  if ( val > 1e10 ) {
	    std::cout << "Overlarge model component " << sourceMap.name() << ' ' << idx << ' ' << val << ' ' 
		      << pf << ' ' << x1 << ' ' << x2 << ' ' 
		      << model[idx] << ' ' << model[idx+npix] << std::endl;
	  }	  
	  modelCounts[idx] = val;
	}
      }
    }  

    void extractModels(const BinnedLikelihood& logLike,
		       const std::string& test_name,	
		       std::vector<std::string>& freeSrcNames,	       
		       std::vector<std::vector<float> >& templates,		       
		       std::vector<float>& fixed,
		       std::vector<float>& test_source_model,
		       std::vector<float>& refPars,
		       std::vector<float>* weights,
		       bool useUnitRefVals) {

      BinnedLikelihood& nc_logLike = const_cast<BinnedLikelihood&>(logLike);

      // These are the energy bin edges
      const std::vector<double>& energies = logLike.energies();

      // Get the size of the models, and allocate the space 
      // in the output vectors
      size_t nbins = logLike.countsMap().data().size();      
      fixed.resize( nbins, 0. );
      test_source_model.resize( nbins, 0. );
      refPars.clear();
      
      // This is just for debugging
      float modelTotal(0.);

      // Use names of all the sources to loop over the sources
      std::vector<std::string> srcNames;
      logLike.getSrcNames(srcNames);

      // Get the weights, if requested
      if ( weights != 0 ) {
	const WeightMap* wts_map = logLike.weightMap();
	if ( wts_map == 0 ) {
	  throw std::runtime_error("Requested to use likelihood weights, but no weights are present in BinnedLikelihood");
	}
	weights->resize( wts_map->model().size() );
	std::copy(wts_map->model().begin(),wts_map->model().end(),weights->begin());
	float sumW(0.);
	sumVector(weights->begin(),weights->end(),sumW);
	std::cout << "Sum of weights " << sumW << " for " << weights->size() << std::endl;
      }

      // Loop over sources
      for ( std::vector<std::string>::const_iterator itr = srcNames.begin();
	    itr != srcNames.end(); itr++ ) {

	// We treat the test source specially
	// To make it faster to turn it on and off
	bool isTest = (*itr) == test_name;

	// We are actually fitting a scale factor w.r.t. the baseline fit
	// so we want to latch the normalization parameter values
	// from the baseline fit
	Source* aSrc = nc_logLike.getSource(*itr);
	double refValue = aSrc->spectrum().normPar().getValue();
	bool hasSourceMap = logLike.hasSourceMap(*itr);

	// Here we extract the model counts
	// Note that we have to it differently for fixed
	// Sources b/c we don't want to have to keep
	// source maps for all the fixed sources
	std::vector<float> modelCounts(nbins, 0.);
	nc_logLike.computeModelMap(*itr,modelCounts);

	if ( !hasSourceMap ) {
	  nc_logLike.eraseSourceMap(*itr);
	} 

	// Here we cache the model
	if ( isTest ) {
	  // We rescale the model for the test source 
	  // to correspond to a value of 1.0 for the normalization 
	  // parameter
	  double scaleFactor = 1./refValue;
	  vectorAdd(test_source_model.begin(),test_source_model.end(),
		    modelCounts.begin(),modelCounts.end(),
		    test_source_model.begin(),test_source_model.end(),
		    1.0,scaleFactor);	  
	} else if ( aSrc->spectrum().normPar().isFree() ) {
	  // The source is free, put it on the vector of templates
	  // and latch the reference value

	  if(useUnitRefVals && refValue > 0) {
	    multiplyByScalar(modelCounts.begin(),modelCounts.end(),1./refValue);
	    refValue = 1.0;
	  }
	  freeSrcNames.push_back(aSrc->getName());
	  templates.resize(templates.size() + 1);
	  std::vector<float>& mt = templates.back();
	  mt.resize(nbins);
	  std::copy(modelCounts.begin(), modelCounts.end(), mt.begin());
	  refPars.push_back(refValue);
	} else {	  
	  // The source is fixed, add it to the total fixed model
          vectorAdd(fixed.begin(),fixed.end(),
                    modelCounts.begin(),modelCounts.end(),
                    fixed.begin(),fixed.end());
	}
	if ( false ) {
	  sumVector(modelCounts.begin(),modelCounts.end(),modelTotal);
	}
      }      
      if ( false ) {
	sumVector(fixed.begin(),fixed.end(),modelTotal);
      }
      return;
    }        


    void extractFixedModel(const BinnedLikelihood& logLike,
			   const std::string& test_name,
			   std::vector<float>& fixed,
			   const std::vector<std::string>* latched){

      BinnedLikelihood& nc_logLike = const_cast<BinnedLikelihood&>(logLike);

      // Get the size of the models, and allocate the space 
      // in the output vectors
      size_t nbins = logLike.countsMap().data().size();      
      fixed.resize( nbins, 0. );
      std::fill(fixed.begin(),fixed.end(),0);

      // This is just for debugging
      float modelTotal(0.);

      // Use names of all the sources to loop over the sources
      std::vector<std::string> srcNames;
      logLike.getSrcNames(srcNames);

      // Loop over sources
      for ( std::vector<std::string>::const_iterator itr = srcNames.begin();
	    itr != srcNames.end(); itr++ ) {

	// We treat the test source specially
	// To make it faster to turn it on and off
	bool isTest = (*itr) == test_name;
	if ( isTest ) continue;

	bool is_latched = false;
	if ( latched ) {
	  std::vector<std::string>::const_iterator itrfind = 
	    std::find(latched->begin(),latched->end(),*itr);
	  if(itrfind != latched->end()) {
	    is_latched = true;
	  }
	}
	if( is_latched ) continue;

	// We are actually fitting a scale factor w.r.t. the baseline fit
	// so we want to latch the normalization parameter values
	// from the baseline fit
	Source* aSrc = nc_logLike.getSource(*itr);
	if ( aSrc->spectrum().normPar().isFree() ) continue;

	bool hasSourceMap = logLike.hasSourceMap(*itr);

	// Here we extract the model counts
	// Note that we have to it differently for fixed
	// Sources b/c we don't want to have to keep
	// source maps for all the fixed sources
	std::vector<float> modelCounts(nbins, 0.);
	nc_logLike.computeModelMap(*itr,modelCounts);

	if( !hasSourceMap ) {
	  nc_logLike.eraseSourceMap(*itr);
	} 

	// The source is fixed, add it to the total fixed model
	vectorAdd(fixed.begin(),fixed.end(),
		  modelCounts.begin(),modelCounts.end(),
		  fixed.begin(),fixed.end());
      }      
      return;
     
    }


    bool extractPrior(const optimizers::Parameter& par,
		      double& centralVal,
		      double& uncertainty) {
      if ( ! par.has_prior() ) return false;
      optimizers::Parameter& ncpar = const_cast<optimizers::Parameter&>(par);
      const optimizers::Function& prior = ncpar.log_prior();
      if ( prior.genericName() != "LogGaussian" && 
	   prior.genericName() != "GaussianError") {
	throw std::runtime_error("Prior is not a LogGaussian or GaussianError");	
      }
      centralVal = prior.getParam("Mean").getTrueValue();
      uncertainty = prior.getParam("Sigma").getTrueValue();
      return true;
    }

    bool extractPriors(const BinnedLikelihood& logLike,
		       const std::vector<std::string>& freeSrcNames,
		       CLHEP::HepVector& centralVals,
		       CLHEP::HepVector& uncertainties,
		       std::vector<bool>& parHasPrior) {
      size_t npar = freeSrcNames.size();
      centralVals = CLHEP::HepVector(npar,0);
      uncertainties = CLHEP::HepVector(npar,0);
      parHasPrior.resize(npar);
      bool retVal(false);
      for ( size_t i(0); i < freeSrcNames.size(); i++ ) {
	const std::string& srcName = freeSrcNames[i];
	const Source& aSrc = logLike.source(srcName);
	const optimizers::Parameter& np = aSrc.spectrum().normPar();
	parHasPrior[i] = extractPrior(np,centralVals[i],uncertainties[i]);
	retVal |= parHasPrior[i];
      }      
      return retVal;
    }

    


    void refactorModels(const std::vector<std::vector<float> >& templates_in,
			const std::vector<float>& fixed_in,
			const std::vector<float>& scales_in,
			const std::vector<bool>& freePars,
			const std::vector<float>* test_source_model,
			std::vector<const std::vector<float>* >& templates_out,
			std::vector<float>& fixed_out,
			std::vector<float>& scales_out) {

      size_t s1 = templates_in.size();
      size_t s2 = freePars.size();
      if ( s1 != s2 ) {
	throw std::runtime_error("Number of templates does not match number of free parameters in FitsUtils::refactorModels");
	return;
      }
      size_t nbins = fixed_in.size();
      fixed_out.resize(nbins,0.);
      std::copy(fixed_in.begin(),fixed_in.end(),fixed_out.begin());

      templates_out.clear();
      scales_out.clear();
      
      for ( size_t i(0); i < s1; i++ ) {
	const std::vector<float>& tmpl = templates_in[i];
	if ( freePars[i] ) {
	  scales_out.push_back(scales_in[i]);
	  templates_out.push_back(&tmpl);	  
	} else {
	  vectorAdd(fixed_out.begin(),fixed_out.end(),
		    tmpl.begin(),tmpl.end(),
		    fixed_out.begin(),fixed_out.end(),
		    1.,scales_in[i]);
	}
      }

      if ( test_source_model )  {
	// Add the test source with a normalization of zero
	templates_out.push_back(test_source_model);
	scales_out.push_back(0.);
      }
    }

    void extractNonZeroBins(const std::vector<float>& dataVect,
                            int npix,
                            std::vector<int>& nonZeroBins,
                            std::vector<float>& dataRed,
                            std::vector<int>& energyBinStopIdxs) {

      size_t n = dataVect.size();
      nonZeroBins.clear();
      dataRed.clear();
      energyBinStopIdxs.clear();
      
      static const float epsilon(1e-15);
      for ( size_t i(0); i < n; i++ ) {
	if ( dataVect[i] > epsilon ) { 
	  nonZeroBins.push_back(i);
	  dataRed.push_back(dataVect[i]);
	}
	// put in placeholders between the energy layers	
	// and latch the index in the reduced vector
        if ( (i+1) % npix == 0 ) {
	  dataRed.push_back(0.);
	  nonZeroBins.push_back(-1*(i+1));
	  energyBinStopIdxs.push_back( dataRed.size() );
	}
      }     
    }

    void sparsifyModel(const std::vector<int>& nonZeroBins,
		       const std::vector<float>& model,
		       std::vector<float>& modelRed) {
      modelRed.clear();
      int nBins(0);
      int i(0);
      int nZero(0);
      float zeroCountSum(0.);
     
      // This is a bit tricky.
      // The real loop is over the bin index in the full model (i)
      // This loop here is just telling us where the next non-zero bin is
      // Then we loop up to that bin.
      // We do things this way to keep track of the contribution
      // from the zero bins in each energy layer
      for ( std::vector<int>::const_iterator itr = nonZeroBins.begin();
	    itr != nonZeroBins.end(); itr++ ) {
	int nextNonZero = *itr;
	if ( *itr < 0 ) {
	  // This is a placeholder between energy layers
	  nextNonZero *= -1;
	} 
	// Loop to the next bin
	for ( ; i < nextNonZero; i++ ) {	    
	  zeroCountSum += model[i];
	  nZero++;
	} 
	if ( *itr < 0 ) {
	  // This is a placeholder between energy layers
	  // push back and reset the sum of the model counts in the zero count bins
	  modelRed.push_back(zeroCountSum);
	  zeroCountSum = 0;	  
	  nBins++;	  
	} else {
	  // This is the next non-zero bin, push it back on the model
	  modelRed.push_back(model[i]);
	  i++;
	}
      }
      // This is just for debugging
      if ( false ) {
	float fullSum(0.);
	float redSum(0.);
	sumVector(model.begin(),model.end(),fullSum);
	sumVector(modelRed.begin(),modelRed.end(),redSum);
	std::cout << "Compare: " << fullSum << ' ' << redSum << std::endl;
	std::cout << "NBins:   " << nBins << ' ' << i << ' ' << nZero << std::endl;
      }

    }

    void sparsifyWeights(const std::vector<int>& nonZeroBins,
			 const std::vector<float>& weights,
			 const std::vector<float>& model,
			 std::vector<float>& weightsRed){
      weightsRed.clear();
      int nBins(0);
      int i(0);
      int nZero(0);
      float zeroCountSum_weight(0.);
      float zeroCountSum_model(0.);
     
      // This is a bit tricky.
      // The real loop is over the bin index in the full model (i)
      // This loop here is just telling us where the next non-zero bin is
      // Then we loop up to that bin.
      // We do things this way to keep track of the contribution
      // from the zero bins in each energy layer
      for ( std::vector<int>::const_iterator itr = nonZeroBins.begin();
	    itr != nonZeroBins.end(); itr++ ) {
	int nextNonZero = *itr;
	if ( *itr < 0 ) {
	  // This is a placeholder between energy layers
	  nextNonZero *= -1;
	} 
	// Loop to the next bin
	for ( ; i < nextNonZero; i++ ) {	    
	  zeroCountSum_weight += (weights[i] * model[i]);
	  zeroCountSum_model += model[i];	  
	  nZero++;
	} 
	if ( *itr < 0 ) {
	  // This is a placeholder between energy layers
	  // push back and reset the sum of the model counts in the zero count bins
	  float zeroCountSum = zeroCountSum_model > 0 ? zeroCountSum_weight / zeroCountSum_model : 1.;
	  weightsRed.push_back(zeroCountSum);
	  zeroCountSum_weight = 0;	  
	  zeroCountSum_model = 0;
	  nBins++;	  
	} else {
	  // This is the next non-zero bin, push it back on the model
	  weightsRed.push_back(weights[i]);
	  i++;
	}
      }
      // This is just for debugging
      if ( false ) {
	float fullSum(0.);
	float redSum(0.);
	sumVector(weights.begin(),weights.end(),fullSum);
	sumVector(weightsRed.begin(),weightsRed.end(),redSum);
	std::vector<float> modelRed;
	sparsifyModel(nonZeroBins,model,modelRed);

	float ip_0(0.); 
	innerProduct(model.begin(),model.end(),weights.begin(),weights.end(),ip_0);
	float ip_1(0.); 
	innerProduct(modelRed.begin(),modelRed.end(),weightsRed.begin(),weightsRed.end(),ip_1);
	  

	std::cout << "Compare:  " << fullSum << ' ' << redSum << std::endl;
	std::cout << "Compare2: " << ip_0 << ' ' << ip_1 << std::endl;
	std::cout << "NBins:    " << nBins << ' ' << i << ' ' << nZero << std::endl;
	
      }

    }


    double logLikePoisson(std::vector<float>::const_iterator data_start, 
			  std::vector<float>::const_iterator data_stop,
			  std::vector<float>::const_iterator model_start,
			  std::vector<float>::const_iterator model_stop) {
      size_t s1 = data_stop - data_start;
      size_t s2 = model_stop - model_start;
      if ( s1 != s2 ) {      
	throw std::runtime_error("Size of model does not match size of data in FitUtils::logLikePoisson.");
	return 0.;
      }
      std::vector<float>::const_iterator itr_data = data_start;
      std::vector<float>::const_iterator itr_model = model_start;
      // We do the running sum in two parts, to avoid numerical issues
      // with adding something very small to something very large 
      // Note also that we are using doubles for the running sums
      double nPred(0.);
      double logTerm(0.);
      for ( ; itr_data != data_stop; itr_data++, itr_model++ ) {
	if ( *itr_model < 0. ) {
	  throw std::runtime_error("Negative model counts in FitUtils::logLikePoisson.");
	  return 0.;
	}
	nPred += *itr_model;
	// logs are expensive, don't do this unless the number of data counts is > 0.
	if ( *itr_data > 1e-9 ) {
	  if (  *itr_model <= 0. ) { 
	    throw std::runtime_error("Negative or zero model counts for pixel with data counts in  FitUtils::logLikePoisson.");
	  }
	  logTerm += ( *itr_data * std::log(*itr_model) );
	}
      }
      return (logTerm-nPred);
    }

    double logLikePoisson(std::vector<float>::const_iterator data_start, 
			  std::vector<float>::const_iterator data_stop,
			  std::vector<float>::const_iterator model_start,
			  std::vector<float>::const_iterator model_stop,
			  std::vector<float>::const_iterator w_start,
			  std::vector<float>::const_iterator w_stop) {
      size_t s1 = data_stop - data_start;
      size_t s2 = model_stop - model_start;
      size_t s3 = w_stop - w_start;
      if ( s1 != s2 ) {      
	throw std::runtime_error("Size of model does not match size of data in FitUtils::logLikePoisson.");
	return 0.;
      }
      if ( s1 != s3 ) {      
	throw std::runtime_error("Size of weights does not match size of data in FitUtils::logLikePoisson.");
	return 0.;
      }
      std::vector<float>::const_iterator itr_data = data_start;
      std::vector<float>::const_iterator itr_model = model_start;
      std::vector<float>::const_iterator itr_w = w_start;
      // We do the running sum in two parts, to avoid numerical issues
      // with adding something very small to something very large 
      // Note also that we are using doubles for the running sums
      double nPred(0.);
      double logTerm(0.);
      for ( ; itr_data != data_stop; itr_data++, itr_model++, itr_w++ ) {
	if ( *itr_model < 0. ) {
	  throw std::runtime_error("Negative model counts in FitUtils::logLikePoisson.");
	  return 0.;
	}
	nPred += (*itr_model) * (*itr_w);
	double dw = (*itr_data) * (*itr_w);
	// logs are expensive, don't do this unless the weighted number of data counts is > 0.
	if ( dw > 1e-9 ) {
	  if (  *itr_model <= 0. ) { 
	    throw std::runtime_error("Negative or zero model counts for pixel with weighted data counts in FitUtils::logLikePoisson.");
	  }
	  logTerm += ( dw * std::log(*itr_model) );
	}
      }
      return (logTerm-nPred);
    }

    int fitModelNorms_newton(const BinnedLikelihood& logLike,
			     const std::string& test_name,
			     double tol, int maxIter, double initLambda,
			     CLHEP::HepVector& norms,
			     CLHEP::HepSymMatrix& covar,
			     CLHEP::HepVector& gradient,
			     std::vector<float>& model,
			     double& edm,
			     double& logLikeVal,
			     size_t firstBin, size_t lastBin) {

      
      std::vector<std::vector<float> > template_master;
      std::vector<std::string> template_names;
      std::vector<float> fixed_master;
      std::vector<float> test_source_master;
      std::vector<float> refPars_master;      

      extractModels(logLike,test_name,template_names,
		    template_master,fixed_master,
		    test_source_master,refPars_master);
      size_t nPar = template_master.size();
      
      std::vector<float> scales_in(nPar,1.);
      std::vector<bool> freePars(nPar,true);
      std::vector<const std::vector<float>* > templates;
      std::vector<float> fixed;
      std::vector<float> scales;
      // EAC FIXME
      std::vector<float>* weights(0);


      refactorModels(template_master,fixed_master,scales_in,freePars,
		     &test_source_master,templates,fixed,scales);
	    
      CLHEP::HepVector initNorms;
      FitScanMVPrior* prior(0);
      covar = CLHEP::HepSymMatrix(scales.size());      
      Vector_Stl_to_Hep(scales,initNorms);
      Vector_Stl_to_Hep(scales,norms);

      model.resize(fixed.size());
      int retVal = fitNorms_newton(logLike.countsMap().data(),
				   initNorms,templates,fixed,prior,weights,
				   tol,maxIter,initLambda,norms,covar,gradient,
				   model,edm,logLikeVal,
				   firstBin,lastBin);
      return retVal;
    }


    void get_spectral_weights(const std::vector<double>& spec,
			      const std::vector<double>& energies,
			      const std::vector<double>& log_energy_ratios,
			      std::vector<std::pair<double, double> >& spec_weights) {
      size_t ne = log_energy_ratios.size();
      spec_weights.clear();
      spec_weights.resize(ne);
      for ( size_t i(0); i < ne; i++ ) {
	spec_weights[i].first = spec[i]*energies[i]*log_energy_ratios[i] / 2.;
	spec_weights[i].second = spec[i+1]*energies[i+1]*log_energy_ratios[i] / 2.;
      }
    }
    

    void get_edisp_range(SourceMap& srcMap,
			 size_t k, 
			 size_t& kmin, size_t& kmax) {

      kmin = k;
      kmax = kmin + 1;
      if ( srcMap.edisp_val() > 0 ) {
	// We are going to scan over energy, so expand the range
	kmin += srcMap.edisp_bins() - srcMap.edisp_val();
	kmax += srcMap.edisp_bins() + srcMap.edisp_val();
      }
    }

    
    void get_edisp_constants(SourceMap& srcMap, 
			     size_t k, size_t& kmin, size_t& kmax,
			     std::vector<double>& edisp_col){
      
      get_edisp_range(srcMap, k, kmin, kmax);
      edisp_col.clear();
      
      if ( srcMap.edisp_val() == 0 ) { 
	// Don't apply energy dispersion
	edisp_col.push_back(1.);
      } else if ( srcMap.edisp_val() < 0 ) { // Apply 'old' version of energy dispersion
	const Drm_Cache* drm_cache = srcMap.drm_cache();
	int kref(-1);
	double xi = drm_cache->get_correction(k, kref);
	edisp_col.push_back(xi);
	if ( kref >= 0 ) { 
	  // Doesn't have counts in this bin, so use the reference bin instead.
	  kmin = kref + srcMap.edisp_bins();
	  kmax = kmin + 1;
	} else {
	  kmin += srcMap.edisp_bins();
	  kmax += srcMap.edisp_bins();
	}
      } else if ( srcMap.edisp_val() > 0 ) {
	const std::vector<double>& drm_col = srcMap.drm()->col(k);
	edisp_col.resize(kmax-kmin);
	// Copy out the relevant energy bins
	std::copy(drm_col.begin() + kmin, drm_col.begin() + kmax, edisp_col.begin());
	// Now adjust kmin and kmax so that they are in the space of the source map
	kmin += srcMap.edisp_offset();
	kmax += srcMap.edisp_offset();
      }
    }


    double model_counts_contribution(SourceMap& srcMap,
				     const std::vector<std::pair<double, double> > & spec_wts,
				     const double& xi, 
				     size_t npix, size_t kref, size_t ipix) {    
      size_t jmin = kref*npix + ipix;
      size_t jmax = jmin + npix;
      double y1 = srcMap[jmin]*spec_wts[kref].first;
      double y2 = srcMap[jmax]*spec_wts[kref].second;
      return xi*(y1 + y2);
    }
  
    void npred_contribution(const std::vector<double>& npred_vals,
			    const std::pair<double,double>& weighted_npreds,
			    const std::vector<std::pair<double, double> > & spec_wts,
			    const double& xi, 
			    size_t kref,
			    double& counts,
			    double& counts_wt) {

      double npred_0 = npred_vals[kref] * spec_wts[kref].first;
      double npred_1 = npred_vals[kref+1] * spec_wts[kref].second;      
      double npred_wt_0 = weighted_npreds.first * spec_wts[kref].first;
      double npred_wt_1 = weighted_npreds.second * spec_wts[kref].second;
      counts = xi*(npred_0 + npred_1);
      counts_wt = xi*(npred_wt_0 + npred_wt_1);
    }


    double model_counts_edisp(SourceMap& srcMap,
			      const std::vector<std::pair<double, double> > & spec_wts,
			      const std::vector<double> & edisp_col,
			      size_t ipix,
			      size_t npix, 
			      size_t kmin,
			      size_t kmax) {
       
      size_t jmin = kmin*npix + ipix;
      size_t jmax = jmin + npix;      
      size_t idx(0);
      double ret_val(0.);
      for ( size_t k(kmin); k < kmax; k++, idx++, jmin+=npix, jmax+=npix ) {
	double y1 = srcMap[jmin] * spec_wts[k].first;
	double y2 = srcMap[jmax] * spec_wts[k].second;
	ret_val += edisp_col[idx] * (y1 + y2);
      }

      return ret_val;
    }			      


    void npred_edisp(const std::vector<double>& npred_vals,
		     const std::vector<std::pair<double,double> >& weighted_npreds,
		     const std::vector<std::pair<double, double> > & spec_wts,
		     const std::vector<double> & edisp_col,
		     size_t kmin,
		     size_t kmax,
		     double& counts,
		     double& counts_wt) {

      counts = 0.;
      counts_wt = 0.;
      size_t idx(0);	
    
      for ( size_t k(kmin); k < kmax; k++, idx++) {
	double y1 = npred_vals[k] * spec_wts[k].first;
	double y2 = npred_vals[k+1] * spec_wts[k].second;
	counts += edisp_col[idx] * (y1 + y2);
	double y1_wt = weighted_npreds[idx].first * spec_wts[k].first;
	double y2_wt = weighted_npreds[idx].second * spec_wts[k].second;		
	counts_wt += edisp_col[idx] * (y1_wt + y2_wt);
      }
    }    
    
    void addSourceCounts(std::vector<double> & modelCounts,
			 SourceMap& srcMap,
			 const BinnedCountsCache& dataCache,
			 int edisp_val,
			 bool subtract) {
      
      double my_sign = subtract ? -1.0 : 1.0;

      const std::vector<std::pair<double, double> >& spec_wts = srcMap.specWts();

      size_t nebins = dataCache.num_ebins();
      size_t npix = dataCache.num_pixels();
      const std::vector<size_t>& pix_ranges = dataCache.firstPixels();
      std::vector<double> edisp_col;

      for (size_t k(0); k < nebins; k++ ) {

	size_t kmin_edisp(0);
	size_t kmax_edisp(0);
	get_edisp_constants(srcMap, k, kmin_edisp, kmax_edisp, edisp_col);	
	size_t j_start = pix_ranges[k];
	size_t j_stop = pix_ranges[k+1];      
	for (size_t j(j_start); j < j_stop; j++) {
	  size_t ipix = dataCache.filledPixels()[j];
	  double counts = edisp_val > 0 ? 
	    model_counts_edisp(srcMap, spec_wts, edisp_col, ipix, npix, kmin_edisp, kmax_edisp) :
	    model_counts_contribution(srcMap, spec_wts, edisp_col[0], npix, kmin_edisp, ipix);
	  double addend = my_sign*counts;
	  modelCounts[j] += addend;
	}
      }
    }

    void addFixedNpreds(std::vector<double>& fixed_counts_spec,
			std::vector<double>& fixed_counts_spec_wt,
			std::vector<double>& fixed_counts_spec_edisp,
			std::vector<double>& fixed_counts_spec_edisp_wt,
			SourceMap& srcMap,
			const BinnedCountsCache& dataCache,
			int edisp_val,
			bool subtract) {
      
      const std::vector<double>& npred_vals = srcMap.npreds();
      const std::vector<std::vector<std::pair<double,double> > >& weighted_npreds = srcMap.weighted_npreds();
      const std::vector<std::pair<double, double> >& spec_wts = srcMap.specWts();

      size_t ne(dataCache.num_ebins());
      std::vector<double> edisp_col;

      std::vector<double> ones(2*srcMap.edisp_bins() +1, 1.);

      for (size_t k(0); k < ne; k++) {

	size_t kmin_edisp(0);
	size_t kmax_edisp(0);
	get_edisp_constants(srcMap, k, kmin_edisp, kmax_edisp, edisp_col);	
	double counts(0.);
	double counts_wt(0.);
	if ( edisp_val > 0 ) {
	  npred_edisp(npred_vals, weighted_npreds.at(k), spec_wts, ones, kmin_edisp, kmax_edisp, counts, counts_wt);
	} else {
	  npred_contribution(npred_vals, weighted_npreds.at(k).at(0), spec_wts, ones[0], kmin_edisp, counts, counts_wt);
	}

	double counts_edisp(0.);
	double counts_edisp_wt(0.);
	if ( edisp_val > 0 ) {
	  npred_edisp(npred_vals, weighted_npreds.at(k), spec_wts, edisp_col, kmin_edisp, kmax_edisp, counts_edisp, counts_edisp_wt);
	} else {
	  npred_contribution(npred_vals, weighted_npreds.at(k).at(0), spec_wts, edisp_col[0], kmin_edisp, counts_edisp, counts_edisp_wt);
	}	

	if ( subtract ) {
	  counts *= -1.;
	  counts_wt *= -1.;
	  counts_edisp *= -1.;
	  counts_edisp_wt *= -1.;
	}

	fixed_counts_spec[k] += counts;
	fixed_counts_spec_wt[k] += counts_wt;
	fixed_counts_spec_edisp[k] += counts_edisp;
	fixed_counts_spec_edisp_wt[k] += counts_edisp_wt;
      }

    }

    void addFreeDerivs(std::vector<Kahan_Accumulator>& posDerivs,
		       std::vector<Kahan_Accumulator>& negDerivs,
		       long freeIndex,
		       SourceMap& srcMap,
		       const std::vector<double>& data_over_model,
		       const BinnedCountsCache& dataCache,
		       int edisp_val,
		       size_t kmin, size_t kmax) {
      
      const std::vector< std::vector<double> > & specDerivs = srcMap.cached_specDerivs();
      const std::vector<double>& energies = srcMap.energies();
      const std::vector<double>& log_energy_ratios = srcMap.log_energy_ratios();


      // We need this stuff for the second term...
      const std::vector<double> & npreds = srcMap.npreds();
      const std::vector<std::vector<std::pair<double,double> > >& weighted_npreds = srcMap.weighted_npreds();
      

      size_t npix = dataCache.num_pixels();
      const std::vector<size_t>& pix_ranges = dataCache.firstPixels();    

      long iparam(freeIndex);
      for (size_t i(0); i < specDerivs.size(); i++, iparam++) {

	std::vector<std::pair<double, double> > spec_wts;
	get_spectral_weights(specDerivs[i], energies, log_energy_ratios, spec_wts);
	std::vector<double> edisp_col;

	for (size_t k(kmin); k < kmax; k++ ) {
	  size_t kmin_edisp(0);
	  size_t kmax_edisp(0);	  
	  get_edisp_constants(srcMap, k, kmin_edisp, kmax_edisp, edisp_col);

	  // First term, derivate of n_obs log n_model = ( n_obs / n_model ) * ( d model / d param ) 
	  // loop over all the filled pixels
	  size_t j_start = pix_ranges[k];
	  size_t j_stop = pix_ranges[k+1];      
	  for (size_t j(j_start); j < j_stop; j++) {
	    if ( data_over_model[j] <= 0. ) {
	      continue;
	    }
	    size_t ipix = dataCache.filledPixels()[j];
	    double counts_deriv = edisp_val > 0 ? 
	      model_counts_edisp(srcMap, spec_wts, edisp_col, ipix, npix, kmin_edisp, kmax_edisp) :
	      model_counts_contribution(srcMap, spec_wts, edisp_col[0], npix, kmin_edisp, ipix);
	    double addend = data_over_model[j]*counts_deriv;
	    if (addend > 0) {
	      posDerivs[iparam].add(addend);
	    } else {
	      negDerivs[iparam].add(addend);
	    }
	  }	

	  // Second term, the derivatives of the nPreds. 
	  // Loop over the energy layers
	  double counts_deriv(0.);
	  double counts_deriv_wt(0.);
	  if ( edisp_val > 0 ) {
	    npred_edisp(npreds, weighted_npreds.at(k), spec_wts, edisp_col, kmin_edisp, kmax_edisp, counts_deriv, counts_deriv_wt);
	  } else {
	    npred_contribution(npreds, weighted_npreds.at(k).at(0), spec_wts, edisp_col[0], kmin_edisp, counts_deriv, counts_deriv_wt);
	  }
	  if (-counts_deriv_wt > 0) {
	    posDerivs[iparam].add(-counts_deriv_wt);
	  } else {
	    negDerivs[iparam].add(-counts_deriv_wt);
	  }
	} // Loop on energy bins
      } // Loop on parameters
    }


    void updateModelMap(std::vector<float> & modelMap,
			SourceMap& srcMap,
			const BinnedCountsCache& dataCache,				       
			bool use_mask) {
      
      const WeightMap* mask = use_mask ? srcMap.weights() : 0;
      size_t npix = dataCache.num_pixels();
      
      const std::vector<std::pair<double, double> >& spec_wts = srcMap.specWts();
      std::vector<double> edisp_col;
 
      for (size_t k(0); k < dataCache.num_ebins(); k++) {	

	size_t kmin_edisp(0);
	size_t kmax_edisp(0);
	get_edisp_constants(srcMap, k, kmin_edisp, kmax_edisp, edisp_col);

	// This is the index in the output model map, which is also the 
	// index in the mask (if there is a mask)
	size_t jmin = k * npix;

	for (size_t ipix(0); ipix < npix; ipix++, jmin++) {
	  // EAC, skip masked pixels
	  if ( mask && 
	       ( mask->model()[jmin] <= 0. ) ) {
	    continue;
	  }
	  double counts = srcMap.edisp_val() >= 0 ?
	    model_counts_edisp(srcMap, spec_wts, edisp_col, ipix, npix, kmin_edisp, kmax_edisp) :
	    model_counts_contribution(srcMap, spec_wts, edisp_col[0], npix, kmin_edisp, ipix);
	  modelMap[jmin] += counts;
	}
      } 
    }  

    void printVector(const std::string& name,
		     const CLHEP::HepVector& vect) {
      std::cout << name << ' ';
      for ( int i(0); i < vect.num_row(); i++ ) {
	std::cout << vect[i] << ' '; 
      }
      std::cout << std::endl;
    }

    void printMatrix(const std::string& name,
		     const CLHEP::HepSymMatrix& mat) {
      std::cout << name << std::endl;
      for ( int i(0); i < mat.num_row(); i++ ) {
	for ( int j(0); j < mat.num_col(); j++ ) {	  
	  std::cout << mat[i][j] << ' '; 
	}
	std::cout << std::endl;
      }
    }

    void printMatrix(const std::string& name,
		     const CLHEP::HepMatrix& mat) {
      std::cout << name << std::endl;
      for ( int i(0); i < mat.num_row(); i++ ) {
	for ( int j(0); j < mat.num_col(); j++ ) {	  
	  std::cout << mat[i][j] << ' '; 
	}
	std::cout << std::endl;
      }
    }
  


  } // namespace FitUtils
 
} // namespace Likelihood
