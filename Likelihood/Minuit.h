/** @file Minuit.h
 * @brief Declaration for the Minuit Optimizer subclass.
 * @author P. Nolan
 *
 */

#ifndef MINUIT_H
#define MINUIT_H

#include "Likelihood/Optimizer.h"
#include "Likelihood/Statistic.h"
//#include "Likelihood/cfortran.h"
//#include "Likelihood/minuitfcn.h"
//#include "Likelihood/minuit.h"

typedef long int logical;  // Copied from f2c.h
typedef short ftnlen;

namespace Likelihood {

/** 
 * @class Minuit
 *
 * @brief Wrapper class for the Minuit optimizer from CERN
 *
 * @author P. Nolan
 *    
 */

class Minuit : public Optimizer {
    
public:
    
  Minuit(Statistic &stat);

  virtual ~Minuit() {}

  void find_min(int verbose = 0, double tol = 1e-3);

  void setMaxEval(int);
  int getQuality(void) const; 
  double getDistance(void) const; 
  std::vector<double> & getUncertainty(void) const;
  
  enum MinuitQuality {MINUIT_NOTCALC, MINUIT_DIAG, MINUIT_FORCEDPOS, 
		      MINUIT_NORMAL};

 private:
  
  Statistic *m_stat;
  int doCmd(std::string command);
  int m_maxEval;
  int m_quality;
  double m_distance;
  std::vector<double> m_uncertainty;
};

 void fcn(int* npar, double* grad, double* fcnval,
	  double* xval, int* iflag, void* futil);
} // namespace Likelihood

extern "C" {
  void mninit_(const int*, const int*, const int*);
  void mnparm_( int *  num , const char * chnam, double * stval, 
	        double * step,  double * bnd1 , 
	        double * bnd2, int * ierror, ftnlen stringlen);
  typedef void (mfcn)(int*, double *, double *, double *, int *, void *);
  void mncomd_(mfcn *, const char * chstr, int * ierr, void * futil, 
	       ftnlen stringlen);
  void mnexcm_(mfcn *, char *, double *, int *, int *, void *, ftnlen);
  void mnpout_(int *, char *, double *, double *, double *, double *, int *,
	       ftnlen);
  void mnstat_(double *, double *, double *, int *, int *, int *);
  void mnseti_(char *, ftnlen);
  void mnpars_(char *, int *, ftnlen);
  void mnemat_(double * emat, int *ndim);
  void mnerrs_(int *, double *, double *, double *, double *);
  void mncont_(mfcn *, int *, int *, int *, double *, double *, int *,
	       void *);
  logical intrac_(double *);
}

#endif // MINUIT_H


