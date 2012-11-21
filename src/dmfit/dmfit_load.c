/* dmfit_load.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c/f2c.h"

/* Common Block Declarations */

struct {
    real phidif[72576]	/* was [252][24][12] */;
    integer hasmooth;
} hasim_;

#define hasim_1 hasim_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;

/* Subroutine */ int dmfit_load__(char *filename, ftnlen filename_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer j, k, l, zn;
    static doublereal milow[12];

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 13, 0, 0, 0 };


/* -----data tables */
/* -----lowest mass index */
    zn = 250;
/* ...lowest mass index for channel j */
    milow[0] = 1.;
/* c c-bar */
    milow[1] = 1.;
/* b b-bar */
    milow[2] = 1.;
/* t t-bar */
    milow[3] = 1.;
/* tau+ tau- */
    milow[4] = 1.;
/* w+ w- */
    milow[5] = 1.;
/* z z */
    milow[6] = 1.;
/* mu+ mu- */
    milow[7] = 1.;
/* gluons */
    milow[8] = 1.;
/* gluons */
    milow[9] = 1.;
/* gluons */
    milow[10] = 1.;
/* gluons */
    milow[11] = 1.;
/* -----clear the tables */
/* gluons */
    for (j = 1; j <= 12; ++j) {
	for (k = 1; k <= 24; ++k) {
	    for (l = 0; l <= 250; ++l) {
		hasim_1.phidif[l + (k + j * 24) * 252 - 6299] = 0.f;
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
    for (j = 1; j <= 12; ++j) {
	for (k = 1; k <= 24; ++k) {
	    if ((doublereal) k >= milow[j - 1]) {
		s_rsle(&io___6);
		i__1 = zn - 1;
		for (l = 0; l <= i__1; ++l) {
		    do_lio(&c__4, &c__1, (char *)&hasim_1.phidif[l + (k + j * 
			    24) * 252 - 6299], (ftnlen)sizeof(real));
		}
		e_rsle();
	    }
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 13;
    cl__1.csta = 0;
    f_clos(&cl__1);
    for (j = 1; j <= 12; ++j) {
	for (k = 1; k <= 24; ++k) {
	    i__1 = zn - 1;
	    for (l = 0; l <= i__1; ++l) {
/* correct units of dyield / dz */
		hasim_1.phidif[l + (k + j * 24) * 252 - 6299] = 
			hasim_1.phidif[l + (k + j * 24) * 252 - 6299];
/*     &            (1.0d0/dble(zn)) */
/* / */
	    }
	    hasim_1.phidif[(k + j * 24) * 252 - 6300] = hasim_1.phidif[(k + j 
		    * 24) * 252 - 6299];
	    hasim_1.phidif[zn + (k + j * 24) * 252 - 6299] = hasim_1.phidif[
		    zn - 1 + (k + j * 24) * 252 - 6299];
	}
    }
/*         dmfit_load=1 */
/* 2000 format(1000(1x,e12.6)) */
    return 0;
} /* dmfit_load__ */

