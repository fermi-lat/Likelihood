/**
 * @file Drmngb.h
 * @brief Declaration for the Drmngb Optimizer subclass.
 * @author P. Nolan, J. Chiang
 *
 */

#ifndef Likelihood_drmngb_h
#define Likelihood_drmngb_h

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
   * @class Drmngb
   *
   * @brief Wrapper class for the MNGB optimization method from 
   * Netlib's PORT library.
   *
   * Finds a local minimum of a continuously differentiable function subject 
   * to simple upper and lower bound constraints. User supplies gradient of 
   * objective function. Secant Hessian approximations are used. Uses a 
   * variant of Newton's method with a quasi-Newton (BFGS) Hessian updating 
   * method, and a model/trust-region technique to aid convergence from 
   * poor starting values.
   *
   * @author P. Nolan, J. Chiang
   * 
   */
  
  /**
   * @file "drmngb_routines.c"
   *
   * @brief Fortran code for DRMNGB translated by f2c
   *
   The original code for DRMNGBB, obtained from Netlib, is in Fortran.
   It was translated using f2c -C++.  The only hand modification required
   was to change \#include "f2c.h" to \#include "f2c/f2c.h" in conformity
   with the way CMT wants files to be laid out.
   
  */
  
  class Drmngb : public Optimizer {
    
  public:
    
    Drmngb(Statistic &stat) 
      {m_stat = &stat;}
    
    virtual ~Drmngb() {}
    
    void find_min(int verbose = 0, double tol = 1e-8);

    //! One-sigma confidence regions based on Hessian, assuming 
    // that this function is a likelihood
    std::vector<double> & getUncertainty(void);

    int getRetCode(void) const;

    enum DrmngbReturnCodes {
      XCONV = 3, RELCONV, BOTHCONV, ABSCONV, SINGCONV, FALSECONV,
      EVALLIM, ITLIM, STOPX, ALLOCATED = 14, LIVSMALL, LVSMALL,
      RESTARTATT, DNEG, VORANGE, CANTCOMPUTEF = 63, BADPAR, 
      CANTCOMPUTEG, BADARR, BADPAR1, BUGS
    };

  private:
    
    Statistic *m_stat;
    std::vector<double> m_uncertainty;
    int m_retCode;

  };

  extern "C" {
    void drmngb_(const double * bounds, double * scale, double * funcval,
		 double * grad, int * iv, const int * liv, const int *lv, 
		 const int * n, double * v, double *x);
    void divset_(const int * kind, int * iv, const int * liv, const int * lv, 
		double * v);
  }
  
} // namespace Likelihood

#endif // Likelihood_lbfgs_h


