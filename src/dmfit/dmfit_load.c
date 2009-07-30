/* dmfit_load.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c/f2c.h"

/* Common Block Declarations */

extern struct {
    real phidif[45360]	/* was [252][18][8] */;
    integer hasmooth;
} hasim_;

#define hasim_1 hasim_

/* Table of constant values */

static integer c__1 = 1;

integer dmfit_load__(filename, filename_len)
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_2000[] = "(1000(1x,e12.6))";

    /* System generated locals */
    integer ret_val, i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), f_clos();

    /* Local variables */
    static integer j, k, l;
    static doublereal milow[8];
    static integer zn;

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 13, 0, fmt_2000, 0 };


/* -----data tables */
/* -----lowest mass index */
    zn = 250;
/* ...lowest mass index for channel j */
    milow[0] = 1.;
/* c c-bar */
    milow[1] = 1.;
/* b b-bar */
    milow[2] = 8.;
/* t t-bar */
    milow[3] = 1.;
/* tau+ tau- */
    milow[4] = 4.;
/* w+ w- */
    milow[5] = 5.;
/* z z */
    milow[6] = 1.;
/* mu+ mu- */
    milow[7] = 1.;
/* -----clear the tables */
/* gluons */
    for (j = 1; j <= 8; ++j) {
	for (k = 1; k <= 18; ++k) {
	    for (l = 0; l <= 250; ++l) {
		hasim_1.phidif[l + (k + j * 18) * 252 - 4787] = (float)0.;
	    }
	}
    }
/* -----load the table for differential flux */
/*        open(unit=13,file='gammamc_dif.dat',status='old', */
    o__1.oerr = 0;
    o__1.ounit = 13;
    o__1.ofnmlen = filename_len;
    o__1.ofnm = filename;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = "formatted";
    o__1.oblnk = 0;
    f_open(&o__1);
    for (j = 1; j <= 8; ++j) {
	for (k = 1; k <= 18; ++k) {
	    if ((doublereal) k >= milow[j - 1]) {
		s_rsfe(&io___6);
		i__1 = zn - 1;
		for (l = 0; l <= i__1; ++l) {
		    do_fio(&c__1, (char *)&hasim_1.phidif[l + (k + j * 18) * 
			    252 - 4787], (ftnlen)sizeof(real));
		}
		e_rsfe();
	    }
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 13;
    cl__1.csta = 0;
    f_clos(&cl__1);
    for (j = 1; j <= 8; ++j) {
	for (k = 1; k <= 18; ++k) {
	    i__1 = zn - 1;
	    for (l = 0; l <= i__1; ++l) {
/* correct units of dyield / dz */
		hasim_1.phidif[l + (k + j * 18) * 252 - 4787] /= 1. / (
			doublereal) zn;
	    }
	    hasim_1.phidif[(k + j * 18) * 252 - 4788] = hasim_1.phidif[(k + j 
		    * 18) * 252 - 4787];
	    hasim_1.phidif[zn + (k + j * 18) * 252 - 4787] = hasim_1.phidif[
		    zn - 1 + (k + j * 18) * 252 - 4787];
	}
    }
    ret_val = 1;
    return ret_val;
} /* dmfit_load__ */

