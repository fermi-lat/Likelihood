/** @file Minuit.h
 * @brief Declaration for the Minuit Optimizer subclass.
 * @author P. Nolan
 *
 */

#ifndef MINUIT_H
#define MINUIT_H

#include "Likelihood/Optimizer.h"
#include "Likelihood/Statistic.h"

// These are copied from f2c.h.  Including f2c.h in C++
// code is problematic.  Some of its macros have weird
// side effects.  
//! f2c implementation of the Fortran LOGICAL type
typedef long int logical; 
//! The type used by f2c and g77 for passing the length 
//! of Fortran CHARACTER variables as an extra parameter
typedef long ftnlen;

namespace Likelihood {

  /** 
   * @class Minuit
   *
   * @brief Wrapper class for the Minuit optimizer from CERN
   *
   * @author P. Nolan
   *    
   */
  
  // Doxygen the C file here so it can be left as nearly as
  // possible in its pristine, machine-produced state.
  /**
   * @file ../src/minuit_routines.c
   *
   * @brief The Minuit package translated from Fortran by f2c
   *
   Minuit is a well-known optimizing/fitting package in the HEP
   community.  It has been developed for over 30 years at CERN.
   
   This file was produced from the CERN Fortran source code.
   First, g77 -E was used to insert all the \#include files and
   make a single, large Fortran file free of preprocessor commands.
   Then f2c -C++ produced this file.  The only hand modification
   required was to change \#include "f2c.h" to \#include "f2c/f2c.h"
   to conform to the way CMT wants files to be laid out.
   
   In non-interactive mode, the API for using Minuit involves the
   functions mninit_, mnseti_, mnparm_, mnpars_, mnexcm_, mncomd_,
   mnpout_, mnstat_, mnemat_, mnerrs_, mncont_, mnintr_, and mninpu_.
  */
  
  class Minuit : public Optimizer {
    
  public:
    
    Minuit(Statistic &stat);
    
    virtual ~Minuit() {}
    
    void find_min(int verbose = 0, double tol = 1e-3);
    //! Override the default maximum number of function evaluations
    void setMaxEval(int);
    //! Minuit return status.   3=OK, 2=forced positive def., 1= not accurate
    int getQuality(void) const; 
    //! Estimated vertical distance from minimum
    double getDistance(void) const; 
    //! One-sigma confidence regions based on Hessian, assuming this function is a likelihood
    std::vector<double> & getUncertainty(void);
    //! Symbolic form of the return codes for readability 
    enum MinuitQuality {MINUIT_NOTCALC, MINUIT_DIAG, MINUIT_FORCEDPOS, 
			MINUIT_NORMAL};

  private:
    
    Statistic *m_stat;
    //! Pass a command string to Minuit
    int doCmd(std::string command);
    int m_maxEval;
    int m_quality;
    double m_distance;
    std::vector<double> m_uncertainty;
  };
  
  //! The function which Minuit will minimize
  void fcn(int* npar, double* grad, double* fcnval,
	   double* xval, int* iflag, void* futil);
} // namespace Likelihood

// The Fortran subroutines which make up the Minuit API
extern "C" {
  //! Initialize Minuit with I/O unit numbers for in, out, save
  void mninit_(const int*, const int*, const int*);
  //! Define a parameter, assigning values and bounds
  void mnparm_( int *  num , const char * chnam, double * stval, 
	        double * step,  double * bnd1 , 
	        double * bnd2, int * ierror, ftnlen stringlen);
  //! Prototype of the function to be minimized.
  typedef void (mfcn)(int * npar, double * grad, double * fval, 
		      double * xval, int * iflag, void * futil);
  //! Execute a Minuit command specified as a character string
  void mncomd_(mfcn * fcn, const char * chstr, int * ierr, void * futil, 
	       ftnlen stringlen);
  //! Execute a Minuit command
  void mnexcm_(mfcn * fcn, char * chcom, double * arglis, int * narg, 
	       int * ierflg, void * futil, ftnlen strln);
  //! Get current value of a parameter
  void mnpout_(int * num, char * chnam, double * val, double * error, 
	       double * bnd1, double * bnd2, int * ivarbl, ftnlen strln);
  //! Get current status of minimization
  void mnstat_(double * fmin, double * fedm, double * errdef, int * npari, 
	       int * nparx, int * istat);
  //! Specify a title for a problem
  void mnseti_(char * ctitle, ftnlen strln);
  //! Define a parameter, assigning values and bounds from variables
  void mnpars_(char * chstr, int * icondn, ftnlen strln);
  //! Get current value of covariance matrix
  void mnemat_(double * emat, int * ndim);
  //! Access current parameter errors
  void mnerrs_(int * num, double * eplus, double * eminus, double * eparab, 
	       double * globcc);
  //! Find a function contour with the MNContour method
  void mncont_(mfcn * fcn, int * num1, int * num2, int * npt, double * xpt, 
	       double * ypt, int * nfound, void * futil);
  //! Utility function used by Minuit: interactive or batch mode
  logical intrac_(double *);
}

#endif // MINUIT_H


