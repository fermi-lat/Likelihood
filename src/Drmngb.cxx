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

    for (;;) {  // Call drmngb_ in an infinite loop
      drmngb_(&paramBounds[0], &scale[0], &funcVal, &gradient[0], &iv[0], &liv,
	     &lv, &nparams, &v[0], &paramVals[0]);
      int value = iv[0];
      if (value == 1) { // request for a function value
	funcVal = -m_stat->value(paramVals);
      }
      else if (value == 2) { // request for the gradient
	m_stat->getFreeDerivs(gradient);
	for (dptr p = gradient.begin();
	     p != gradient.end(); p++) {
	  *p = -*p;
	}
      }
      else {  // Finished.  Exit loop.
	m_retCode = value;
	if (value > 6) {throw Exception("DRMNGB error", value);}
	break;
      }
    }

    // Get parameter values
    int j = 0;
    for (pptr p = params.begin();
	 p != params.end(); p++, j++) {
      p->setValue(paramVals[j]);
    }
    // Put parameter values back into the Statistic
    (*m_stat)(paramVals);

    // Get the Cholesky factor of the Hessian.  It's a triangular
    // matrix stored in compact fashion, so we treat it as 1-dimensional
    dptr vp = v.begin() + iv[41] - 1;
    std::vector<double> hess(vp, vp + nparams*(nparams+1)/2);

    // Undo the Cholesky factorization for the diagonal elements to
    // find the parameter uncertainties.
    m_uncertainty.clear();
    int rowLength = 1;  // Each row is longer than the one before
    for (dptr p = hess.begin(); p != hess.end(); /* */ ) { // loop over rows
      double diag = 0.;
      dptr endRow = p + rowLength;
      while (p != endRow) { // Diagonal = sum of factored row^2
	diag += (*p) * (*p);
	p++;
      }
      rowLength++;
      // Since this is a likelihood function, 1-sigma uncertainty
      // is 1/sqrt(Hessian diagonals)
      m_uncertainty.push_back(1./sqrt(diag));
    }
    
  } // End of find_min
}


