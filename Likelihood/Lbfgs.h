/**
 * @file Lbfgs.h
 * @brief Declaration for the Lbfgs Optimizer subclass.
 * @author P. Nolan, J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Lbfgs.h,v 1.8 2003/06/05 21:39:48 pln Exp $
 */

#ifndef Likelihood_lbfgs_h
#define Likelihood_lbfgs_h

#include "Likelihood/Optimizer.h"
#include "Likelihood/Statistic.h"
#include <string>

namespace Likelihood {
  
  // These typedefs are copied from f2c.h.  It is problematic to 
  //include f2c.h in C++ code.  It has some macros which can have
  // strange effects.

  //! f2c implementation of the fortran LOGICAL type
  typedef long int logical; 

  //! type used by f2c to pass length of fortran string
  typedef long ftnlen;
  
  /** 
   * @class Lbfgs
   *
   * @brief Wrapper class for the Broyden-Fletcher-Goldfarb-Shanno
   * variable metric implementation of Byrd, Lu, Nocedal, & Zhu 1995,
   * SIAM, J. Sci. Comp., 16, 5 (http://www.netlib.org/opt/lbfgs_bcm.shar).
   *
   * @author P. Nolan, J. Chiang
   * 
   * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/Likelihood/Lbfgs.h,v 1.8 2003/06/05 21:39:48 pln Exp $
   */
  
  /**
   * @file "lbfgs_routines.c"
   *
   * @brief Fortran code for LBFGS-B translated by f2c
   *
   The original code for LBFGS-B, obtained from Netlib, is in Fortran.
   It was translated using f2c -C++.  The only hand modification required
   was to change \#include "f2c.h" to \#include "f2c/f2c.h" in conformity
   with the way CMT wants files to be laid out.
   
   The L in the name means "limited memory".  That means that the full
   approximate Hessian is not available, only some updates which can
   be used to generate relevant portions of it on the fly.  This can be
   a shortcoming for our purposes because when we maximize a likelihood
   function we would like to estimate uncertainties and correlations
   from the inverse Hessian.
   
   Lbfgs requires its calling routine to provide all the temporary
   storage arrays.  Thus there are no artificial limitations on the
   number of variables.  It uses the "reverse communication" style
   of API, so it isn't necessary to write the objective function in
   any specified style.
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

    //! Number of variable metric corrections to save
    int m_maxVarMetCorr; 

    //! Stop after this many function evaluations
    int m_maxIterations; 

    //! One of the stopping criteria
    double m_pgtol; 
    int m_retCode;

    int m_numEvals;
    std::string m_errorString;
  };
  
  extern "C" {
    //! The Fortran subroutine that controls the LBFGS minimizer
    void setulb_(const int *n, const int *m, double *x, const double *l, 
		 const double *u, const int *nbd, double *f, double *g, 
		 const double *factr, const double *pgtol, double *wa, 
		 int *iwa, char *task, const int *iprint, char *csave, 
		 logical *lsave, int *isave, double *dsave, 
		 ftnlen task_len, ftnlen csave_len);
  }
  
} // namespace Likelihood

#endif // Likelihood_lbfgs_h


