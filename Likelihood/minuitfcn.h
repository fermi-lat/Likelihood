/* C interface to minuit function */
/* Andr\'as Major 1997 */

#include "cfortran.h"

void fcn(int npar, double* grad, double* fcnval, double* xval, int iflag, void* futil);

#ifdef NO_CFORTRAN_HACK
void minuitfcn(int npar, double* grad, double* fcnval, double* xval, int iflag, int* futil);
FCALLSCSUB6(minuitfcn,MINUITFCN,minuitfcn,INT,DOUBLEV,PDOUBLE,DOUBLEV,INT,PINT);
void minuitfcn(int npar, double* grad, double* fcnval, double* xval, int iflag, int* futil)
{
  fcn(npar, grad, fcnval, xval, iflag, (void *) futil);
}
#else
void minuitfcn(int* npar, double* grad, double* fcnval, double* xval, int* iflag, int* futil);
FCALLSCSUB6(minuitfcn,MINUITFCN,minuitfcn,PINT,DOUBLEV,PDOUBLE,DOUBLEV,PINT,PINT)
void minuitfcn(int* npar, double* grad, double* fcnval, double* xval, int* iflag, int* futil)
{
  fcn(*npar, grad, fcnval, xval, *iflag, (void *) futil);
}
#endif

int intrac_()
{
  return C2FLOGICAL(FALSE);
}
