/**
 * @file Lbfgs.h
 * @brief Declaration for the Lbfgs Optimizer subclass.
 * @author P. Nolan, J. Chiang
 *
 * $Header$
 */

#ifndef lbfgs_h
#define lbfgs_h

#include "Likelihood/Optimizer.h"
#include "Likelihood/Statistic.h"
#include <string>

namespace Likelihood {
  typedef long int logical;  // Copied from f2c.h
  typedef short ftnlen;

/** 
 * @class Lbfgs
 *
 * @brief Wrapper class for the Broyden-Fletcher-Goldfarb-Shanno
 * variable metric implementation of Byrd, Lu, Nocedal, & Zhu 1995,
 * SIAM, J. Sci. Comp., 16, 5 (http://www.netlib.org/opt/lbfgs_bcm.shar).
 *
 * @author P. Nolan, J. Chiang
 * 
 * $Header$
 */

class Lbfgs : public Optimizer {
    
public:
    
  Lbfgs(Statistic &stat) : m_maxVarMetCorr(5),
    m_maxIterations(100),
    m_pgtol(.00001),
    m_retCode(0)
    {m_stat = &stat;}

  virtual ~Lbfgs() {}

  void setMaxVarMetCorr(const int m);
  void setPgtol(const double pgtol);
  void setMaxIterations(const int iterations);
  
  int getRetCode(void) const;
  std::string getErrorString(void) const;
  
  void find_min(int verbose = 0, double tol = 1e-5);
  
  enum LbfgsReturnCodes {LBFGS_NORMAL, LBFGS_ABNO, LBFGS_ERROR,
			 LBFGS_TOOMANY, LBFGS_UNKNOWN};

 private:
  
  Statistic *m_stat;
  int m_maxVarMetCorr;  // # of variable metric corrections to save
  int m_maxIterations;  // Stop after this many function evaluations
  double m_pgtol;       // One of the stopping criteria
  int m_retCode;
  int m_numEvals;
  std::string m_errorString;
};

// The Fortran subroutine that controls the LBFGS minimizer
 extern "C" {
   void setulb_(const int *n, const int *m, double *x, const double *l, 
		const double *u, const int *nbd, double *f, double *g, 
		const double *factr, const double *pgtol, double *wa, 
		int *iwa, char *task, const int *iprint, char *csave, 
		logical *lsave, int *isave, double *dsave, 
		ftnlen task_len, ftnlen csave_len);
 }

} // namespace Likelihood

#endif // lbfgs_h


