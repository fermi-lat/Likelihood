/** @file Lbfgs.cxx
 * @brief Lbfgs class implementation
 *
 */

#include "Likelihood/Lbfgs.h"
#include "Likelihood/Parameter.h"
#include "LikelihoodException.h"
#include <vector>
#include <algorithm>
#include <iostream>

namespace Likelihood {
  
  void Lbfgs::setMaxVarMetCorr(const int m)
  {m_maxVarMetCorr = m;}
  
  void Lbfgs::setPgtol(const double pgtol)
  {m_pgtol = pgtol;}

  void Lbfgs::setMaxIterations(const int mx)
  {m_maxIterations = mx;}

  int Lbfgs::getRetCode(void) const
  {return m_retCode;}

  std::string Lbfgs::getErrorString(void) const
  {return m_errorString;}

  void Lbfgs::find_min(int verbose, double tol) {

    m_numEvals = 0;
    m_errorString.erase();

    // Unpack model parameters into the arrays needed by LBFGS
    
    std::vector<Parameter> params;
    m_stat->getFreeParams(params);
    const int nparams = params.size();
    
    std::vector<double> paramVals(nparams);
    std::vector<double> paramMins(nparams);
    std::vector<double> paramMaxs(nparams);
    int i=0;
    for (std::vector<Parameter>::iterator p = params.begin();
	 p != params.end(); p++, i++) {
      paramVals[i] = p->getValue();
      paramMins[i] = p->getBounds().first;
      paramMaxs[i] = p->getBounds().second;
    }
    
    // Create the variables and arrays used by LBFGS
    // These serve as storage between calls to setulb_,
    // so they must be declared outside the loop.
    // Most of them don't need to be initialized.

    double funcVal;
    std::vector<double> gradient(nparams);
    std::vector<int> isave(44);
    std::vector<logical> lsave(4);
    std::vector<double> dsave(29);
    std::vector<char> csave(60, ' ');  // Blank-filled
    std::vector<int> intWorkArray(3*nparams);
    int workSize = (2*m_maxVarMetCorr + 4) * nparams
      + 12 * m_maxVarMetCorr * (m_maxVarMetCorr + 12);
    std::vector<double> doubleWorkArray(workSize);

    // Fortran-style string used for communication with LBFGS

    std::vector<char> task(60, ' '); // Blank-filled
    static const std::string strt("START");
    std::copy(strt.begin(), strt.end(), task.begin()); // Initialize
    
    // Call LBFGS in an infinite loop.  Break out when it's done.

    for (;;) {
      int iprint = verbose - 2;  
      double factr = tol * 1.0e+16; // One of the stopping criteria
      const std::vector<int> nbd(nparams, 2); // All params bounded for now
      setulb_(&nparams, &m_maxVarMetCorr, &paramVals[0], &paramMins[0], 
	      &paramMaxs[0], &nbd[0], &funcVal, &gradient[0], 
	      &factr, &m_pgtol, &doubleWorkArray[0], &intWorkArray[0], 
	      &task[0], &iprint, &csave[0], &lsave[0], &isave[0], &dsave[0], 
	      task.size(), csave.size());
      std::string taskString(task.begin(), task.end());
      int taskLength = taskString.find_last_not_of(' ') + 1;
      taskString.erase(taskLength); // Strip trailing blanks

      if (taskString.substr(0,2) == "FG") {
	// Request for values of function and gradient.
	// LBFGS is a minimizer, so we must flip the signs to maximize.
	funcVal = -m_stat->value(paramVals);
	m_stat->getFreeDerivs(gradient);
	for (int i = 0; i < nparams; i++) {
	  gradient[i] = -gradient[i];
	}
	m_numEvals++;
	
	if (verbose != 0) {
	  std::cout << "LBFGS " << funcVal << "  ";
	  for (int i = 0; i < nparams; i++) {
	    std::cout << paramVals[i] << "  " << gradient[i] << " : ";
	  }
	  std::cout << std::endl;
	}
      }  // Don't break.  Call setulb_ again
      else if (taskString.substr(0, 5) == "NEW_X") {
	// Ready to move to a new set of parameter values
	if (isave[33] > m_maxIterations) {
	  m_retCode = LBFGS_TOOMANY;
	  m_errorString = "Exceeded Specified Number of Iterations";
	  break;
	}
      }  // Otherwise don't break.  Call setulb_ again.
      else if (taskString.substr(0, 4) == "CONV") {
	// Normal convergence
	m_retCode = LBFGS_NORMAL;
	m_errorString = taskString;
	break;
      }
      else if (taskString.substr(0, 4) == "ABNO") {
	// Abnormal termination in line search
	m_retCode = LBFGS_ABNO;
	m_errorString = taskString;
	throw LikelihoodException(taskString, LBFGS_ABNO);
      }
      else if (taskString.substr(0, 5) == "ERROR") {
	// Error in input parameters
	m_retCode = LBFGS_ERROR;
	m_errorString = taskString;
	throw LikelihoodException(taskString, LBFGS_ERROR);
      }
      else {
	// Something else
	m_retCode = LBFGS_UNKNOWN;
	m_errorString = "LBFGS unknown condition";
	throw LikelihoodException(taskString, LBFGS_UNKNOWN);
      }
    }  // End of infinite loop

    // Get parameter values
    int j = 0;
    for (std::vector<Parameter>::iterator p = params.begin();
	 p != params.end(); p++, j++) {
      p->setValue(paramVals[j]);
    }
    // Put parameter values back into the Statistic
    (*m_stat)(paramVals);
  } // End of find_min
}


