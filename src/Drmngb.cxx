/** 
 * @file Drmngb.cxx
 * @brief Drmngb class implementation
 * @author P. Nolan
 *
 */

#include "Likelihood/Drmngb.h"
#include "Likelihood/Parameter.h"
#include "Likelihood/Exception.h"
#include <vector>
#include <algorithm>
#include <iostream>

namespace Likelihood {

  typedef std::vector<double>::iterator dptr;
  typedef std::vector<Parameter>::iterator pptr;

  std::vector<double> & Drmngb::getUncertainty(void) {
    return m_uncertainty;
  }

  int Drmngb::getRetCode(void) const {
    return m_retCode;
  }
  
  void Drmngb::find_min(int verbose, double tol) {

    // Unpack model parameters into the arrays needed by Drmngb
    
    std::vector<Parameter> params;
    m_stat->getFreeParams(params);
    const int nparams = params.size();
    
    std::vector<double> paramVals;
    std::vector<double> paramBounds;
    for (pptr p = params.begin(); p != params.end(); p++) {
      paramVals.push_back(p->getValue());
      paramBounds.push_back(p->getBounds().first);
      paramBounds.push_back(p->getBounds().second);
    }
    
    // Create the variables and arrays used by DRMNGB.
    // These serve as storage between calls to drmngb
    // so they must be declared outside the loop.
    // Most of them don't need to be initialized.

    double funcVal;
    const int liv = 59 + nparams;
    const int lv = 71 + nparams*(nparams+19)/2;
    std::vector<double> gradient(nparams), scale(nparams, 1.);
    std::vector<int> iv(liv);
    std::vector<double> v(lv);

    // Set default values for internal settings
    const int kind = 2;
    divset_(&kind, &iv[0], &liv, &lv, &v[0]);
    if (verbose == 0) {iv[20] = 0;}
    v[31] = tol;

    // Call the optimizing function in an infinite loop.
    for (;;) {
      drmngb_(&paramBounds[0], &scale[0], &funcVal, &gradient[0], &iv[0], &liv,
	     &lv, &nparams, &v[0], &paramVals[0]);
      int value = iv[0];
      if (value == 1) { // request for a function value
	funcVal = -m_stat->value(paramVals);
      }
      else if (value == 2) { // request for the gradient
	m_stat->getFreeDerivs(gradient);
	for (dptr p = gradient.begin(); p != gradient.end(); p++) {
	  *p = -*p;
	}
      }
      else {  // Finished.  Exit loop.
	m_retCode = value;
	if (value > 6) {throw Exception("DRMNGB error", value);}
	break;
      }
    }

    // Get parameter values and put them back into the Statistic
    int j = 0;
    for (pptr p = params.begin(); p != params.end(); p++, j++) {
      p->setValue(paramVals[j]);
    }
    (*m_stat)(paramVals);

    // Get the Cholesky factor of the Hessian.  It's a lower triangular
    // matrix stored in compact fashion, so we treat it as 1-dimensional.
    const dptr vp = v.begin() + iv[41] - 1;
    std::vector<double> hess(vp, vp + nparams*(nparams+1)/2);

    // Invert the Hessian to produce the covariance matrix.  The result 
    //is also stored as a compact lower triangular matrix.
    int info;
    const char uplo = 'L';
    dpptri_(&uplo, &nparams, &hess[0], &info, 1);

    // Extract the diagonal elements of the covariance matrix.
    // Their square roots are the parameter uncertainties.
    m_uncertainty.clear();
    for (int i = 0; i < nparams; i++) {
      m_uncertainty.push_back(sqrt(hess[i*(i+3)/2]));
    }

  } // End of find_min
}


