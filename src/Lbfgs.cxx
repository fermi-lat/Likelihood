/** @file Lbfgs.cxx
 * @brief Lbfgs class implementation
 *
 */

#include "Likelihood/Lbfgs.h"
#include "Likelihood/Parameter.h"
#include <vector>
#include <utility>
#include <cmath>
#include <stdlib.h>
#include <string.h>

namespace Likelihood {

  void Lbfgs::setMaxVarMetCorr(const int m)
  {m_maxVarMetCorr = m;}

  void Lbfgs::setPgtol(const double pgtol)
  {m_pgtol = pgtol;}

  void Lbfgs::setMaxIterations(const int mx)
  {m_maxIterations = mx;}

  int Lbfgs::getRetCode(void) const
  {return m_retCode;}

  void Lbfgs::find_min(int verbose, double tol) {
    
    m_verbose = verbose;
    
    std::vector<Parameter> Params;
    m_stat->getFreeParams(Params);
    
    const int nparams = Params.size();
    double *paramVals = new double[nparams];
    double *paramMins = new double[nparams];
    double *paramMaxs = new double[nparams];
    int iprint = verbose - 2;
    m_numEvals = 0;
    
    for (int i = 0; i < nparams; i++) {
      paramVals[i] = Params[i].getValue();
      paramMins[i] = Params[i].getBounds().first;
      paramMaxs[i] = Params[i].getBounds().second;
    }
    
    lbfgsMin(nparams, paramVals, paramMins, paramMaxs,
			  iprint, tol);
    
    delete [] paramVals;
    delete [] paramMins;
    delete [] paramMaxs;
  }
  
  void Lbfgs::lbfgsMin(const int n, double *x, const double *xmin, 
		      const double *xmax, const int iprint, 
		      const double tol)
  {
    double factr = tol * 1.0e+16;
    
    int* nbd = new int[n];    
    for (int i=0; i < n; i++) {
      nbd[i] = 2;  //  All parameters are bounded, for now.
    }
    
    char task[60], csave[60];
    memset(task,  ' ', sizeof(task));
    memset(csave, ' ', sizeof(csave));
    logical* lsave = new logical[4];
    int* intWorkArray = new int[3*n];
    int* isave = new int[44];
    double funcVal, fmin;
    double* gradient = new double[n];
    double* xopt = new double[n];
    double* dsave = new double[29];
    double* doubleWorkArray = new 
      double[(2*m_maxVarMetCorr+4)*n+12*m_maxVarMetCorr*(m_maxVarMetCorr+12)];
    strncpy(task, "START", strlen("START"));
    
    for (;;) {
      setulb_(&n, &m_maxVarMetCorr, x, xmin, xmax, nbd, &funcVal, gradient, 
	      &factr, &m_pgtol, doubleWorkArray, intWorkArray, 
	      task, &iprint, csave, lsave, isave, dsave, 
	      sizeof(task), sizeof(csave));
      if (strncmp(task, "FG", strlen("FG")) == 0) {
	std::vector<double>paramValues(n); // evaluate at trial point
	for (int i = 0; i < n; i++) {
	  paramValues[i] = x[i];
	}
	
	funcVal = -m_stat->value(paramValues);
	
	std::vector<double> derivsVec;
	m_stat->getFreeDerivs(derivsVec);
	
	for (int i = 0; i < n; i++) 
	  {
	    gradient[i] = -derivsVec[i];
	  }
	
	if (m_verbose) {
	  std::cout << funcVal << "  ";
	  for (int i = 0; i < n; i++)
	    std::cout << x[i] << "  "
		      << gradient[i] << " : ";
	  std::cout << std::endl;
	}
	if (funcVal < fmin || m_numEvals == 0) {
	  fmin = funcVal;
	  for (int i=0; i < n; i++) {
	    xopt[i] = x[i];
	  }
	}
	m_numEvals++;
      }
      else if (strncmp(task, "NEW_X", strlen("NEW_X")) == 0) {
	for (int i=0; i < n; i++) {
	  x[i] = xopt[i];
	}
	if (isave[33] > m_maxIterations) {
	  m_retCode = LBFGS_TOOMANY;
	  break;
	}
      }
      else if (strncmp(task, "ABNO", strlen("ABNO")) == 0) {
	m_retCode = LBFGS_ABNO;
	break;
      }
      else if (strncmp(task, "ERROR", strlen("ERROR")) == 0) {
	m_retCode = LBFGS_ERROR;
	break;
      }
      else {
	m_retCode = LBFGS_NORMAL;
	break;
      }
    }
    delete [] nbd;
    delete [] lsave;
    delete [] dsave;
    delete [] intWorkArray;
    delete [] doubleWorkArray;
    delete [] gradient;
    delete [] xopt;
  }
}
