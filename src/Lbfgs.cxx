/** @file Lbfgs.cxx
 * @brief Lbfgs class implementation
 *
 */

#include "Likelihood/Lbfgs.h"
#include "Likelihood/Parameter.h"
#include <vector>
#include <utility>
#include <cmath>

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

    int iprint = verbose - 2;
    m_numEvals = 0;
    double factr = tol * 1.0e+16;
    m_errorString.erase();

    // Unpack model parameters into the arrays needed by LBFGS
    
    std::vector<Parameter> Params;
    m_stat->getFreeParams(Params);
    const int nparams = Params.size();
    
    std::vector<double> paramVals;
    std::vector<double> paramMins;
    std::vector<double> paramMaxs;
    paramVals.reserve(nparams);
    paramMins.reserve(nparams);
    paramMaxs.reserve(nparams);

    for (std::vector<Parameter>::iterator p = Params.begin();
	   p != Params.end(); p++) {
      paramVals.push_back(p->getValue());
      paramMins.push_back(p->getBounds().first);
      paramMaxs.push_back(p->getBounds().second);
    }
    
    // Create the variables and arrays used by LBFGS

    const std::vector<int> nbd(nparams, 2); // All params bounded for now
    std::vector<logical> lsave(4);
    std::vector<int> intWorkArray(3*nparams);
    std::vector<int> isave(44);
    double funcVal;
    std::vector<double> gradient(nparams);
    std::vector<double> dsave(29);
    int workSize = (2*m_maxVarMetCorr + 4) * nparams
      + 12 * m_maxVarMetCorr * (m_maxVarMetCorr + 12);
    std::vector<double> doubleWorkArray(workSize);

    // Fortran-style strings used for communication with LBFGS

    std::vector<char> task(60, ' '), csave(60, ' '); // Blank-filled
    std::string strt("START"), fg("FG"), newx("NEW_X"), abno("ABNO");
    std::string errr("ERROR"), conv("CONV");
    std::copy(strt.begin(), strt.end(), task.begin());
    std::string taskString;
    
    // Call LBFGS in an infinite loop.  Break out when it's done.

    for (;;) {
      setulb_(&nparams, &m_maxVarMetCorr, &paramVals[0], &paramMins[0], 
	      &paramMaxs[0], &nbd[0], &funcVal, &gradient[0], 
	      &factr, &m_pgtol, &doubleWorkArray[0], &intWorkArray[0], 
	      &task[0], &iprint, &csave[0], &lsave[0], &isave[0], &dsave[0], 
	      task.size(), csave.size());
      taskString = std::string(task.begin(), task.end());
      int taskLength = taskString.find_last_not_of(' ') + 1;
      taskString = taskString.substr(0, taskLength);

      if (std::equal(fg.begin(), fg.end(), task.begin())) {
	// Request for values of function and gradient.
	// Flip signs so we are maximizing the function.
	funcVal = -m_stat->value(paramVals);
	m_stat->getFreeDerivs(gradient);
	std::transform(gradient.begin(), gradient.end(), gradient.begin(),
		       negate<double>());
	m_numEvals++;

	if (verbose != 0) {
	  std::cout << funcVal << "  ";
	  for (int i = 0; i < nparams; i++) {
	    std::cout << paramVals[i] << "  " << gradient[i] << " : ";
	  }
	  std::cout << std::endl;
	}
      }
      else if (std::equal(newx.begin(), newx.end(), task.begin())) {
	// Ready to move to a new set of parameter values
	if (isave[33] > m_maxIterations) {
	  m_retCode = LBFGS_TOOMANY;
	  m_errorString = "Exceeded Specified Number of Iterations";
	  break;
	}
      }
      else if (std::equal(abno.begin(), abno.end(), task.begin())) {
	// Abnormal termination in line search
	m_retCode = LBFGS_ABNO;
	m_errorString = taskString;
	break;
      }
      else if (std::equal(errr.begin(), errr.end(), task.begin())) {
	// Some other type of error
	m_retCode = LBFGS_ERROR;
	m_errorString = taskString;
	break;
      }
      else if (std::equal(conv.begin(), conv.end(), task.begin())) {
	// Normal convergence
	m_retCode = LBFGS_NORMAL;
	m_errorString = "Normal convergence";
	break;
      }
      else {
	// Something else
	m_retCode = LBFGS_UNKNOWN;
	m_errorString = "LBFGS unknown condition";
	break;
      }
    }  // End of infinite loop
  } // End of find_min
}


