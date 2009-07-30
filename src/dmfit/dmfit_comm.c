/* dmfit_comm.f -- translated by f2c (version 20031025).
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

extern struct {
    real phidif[45360]	/* was [252][18][10] */;
    integer hasmooth;
} hasim_;

#define hasim_1 hasim_

/* ******************************************************************** */
/* ******************************************************************** */
/* ** the smoothing is controlled by the parameter hasmooth in the */
/* ** following manner. */
/* ** hasmooth = 0 - no smoothing */
/* **            1 - smoothing of zi-1,zi and zi+1 bins if z>.3 */
/* **            2 - smoothing of zi-2,zi-1,zi,zi+1 and zi+2 if z>0.3 */
/* **************************************************************************** */
doublereal yieldget_(integer *zi, integer *mxi, integer *ch)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer zn;

    zn = 250;
/* -----only differential flux */
    if (hasim_1.hasmooth == 1) {
	if (*zi >= 1 && *zi <= zn - 1) {
	    ret_val = (doublereal) hasim_1.phidif[*zi - 1 + (*mxi + *ch * 18) 
		    * 252 - 4787] * .25 + (doublereal) hasim_1.phidif[*zi + (*
		    mxi + *ch * 18) * 252 - 4787] * .5 + (doublereal) 
		    hasim_1.phidif[*zi + 1 + (*mxi + *ch * 18) * 252 - 4787] *
		     .25;
	} else if (*zi == 0) {
	    ret_val = (doublereal) hasim_1.phidif[*zi + (*mxi + *ch * 18) * 
		    252 - 4787] * .75 + (doublereal) hasim_1.phidif[*zi + 1 + 
		    (*mxi + *ch * 18) * 252 - 4787] * .25;
	} else if (*zi == zn) {
	    ret_val = (doublereal) hasim_1.phidif[*zi + (*mxi + *ch * 18) * 
		    252 - 4787] * .75 + (doublereal) hasim_1.phidif[*zi - 1 + 
		    (*mxi + *ch * 18) * 252 - 4787] * .25;
	}
    } else if (hasim_1.hasmooth == 2) {
	if (*zi <= zn - 2) {
	    ret_val = (doublereal) hasim_1.phidif[*zi - 2 + (*mxi + *ch * 18) 
		    * 252 - 4787] * .1 + (doublereal) hasim_1.phidif[*zi - 1 
		    + (*mxi + *ch * 18) * 252 - 4787] * .225 + (doublereal) 
		    hasim_1.phidif[*zi + (*mxi + *ch * 18) * 252 - 4787] * 
		    .35 + (doublereal) hasim_1.phidif[*zi + 1 + (*mxi + *ch * 
		    18) * 252 - 4787] * .225 + (doublereal) hasim_1.phidif[*
		    zi + 2 + (*mxi + *ch * 18) * 252 - 4787] * .1;
	} else if (*zi == zn - 1) {
	    ret_val = (doublereal) hasim_1.phidif[*zi - 2 + (*mxi + *ch * 18) 
		    * 252 - 4787] * .1 + (doublereal) hasim_1.phidif[*zi - 1 
		    + (*mxi + *ch * 18) * 252 - 4787] * .225 + (doublereal) 
		    hasim_1.phidif[*zi + (*mxi + *ch * 18) * 252 - 4787] * 
		    .45 + (doublereal) hasim_1.phidif[*zi + 1 + (*mxi + *ch * 
		    18) * 252 - 4787] * .225;
	} else {
	    ret_val = (doublereal) hasim_1.phidif[*zi - 2 + (*mxi + *ch * 18) 
		    * 252 - 4787] * .1 + (doublereal) hasim_1.phidif[*zi - 1 
		    + (*mxi + *ch * 18) * 252 - 4787] * .225 + (doublereal) 
		    hasim_1.phidif[*zi + (*mxi + *ch * 18) * 252 - 4787] * 
		    .675;
	}
    } else {
	ret_val = (doublereal) hasim_1.phidif[*zi + (*mxi + *ch * 18) * 252 - 
		4787];
    }
    return ret_val;
} /* yieldget_ */

/* ******************************************************************** */
/* ******************************************************************** */
/* Subroutine */ int ifind_(doublereal *value, doublereal *array, doublereal *
	ipl, integer *ii, integer *imin, integer *imax)
{
    /* System generated locals */
    integer array_offset;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, iold, inew, imint, imaxt;

/* +10 to avoid sign pro */
    /* Parameter adjustments */
    array_offset = *imin + 10;
    array -= array_offset;

    /* Function Body */
    imint = *imin + 10;
    imaxt = *imax + 10;
    if (*value < array[imint] || *value >= array[imaxt + 1]) {
	*ii = -5;
	return 0;
    }
    iold = 0;
    inew = 0;
    i__ = (imaxt + imint) / 2;
L10:
    if (*value >= array[i__] && *value < array[i__ + 1]) {
	*ii = i__ - 10;
	*ipl = (*value - array[i__]) / (array[i__ + 1] - array[i__]);
	return 0;
    }
    if (*value > array[i__]) {
	inew = (i__ + 1 + imaxt) / 2;
	imint = i__;
    } else {
	inew = (imint + i__) / 2;
	imaxt = i__;
    }
    i__ = inew;
    if (iold == inew) {
	s_stop("", (ftnlen)0);
    }
    iold = inew;
    goto L10;
} /* ifind_ */

/* ******************************************************************** */
doublereal llg_(doublereal *x, doublereal *mx, doublereal *ml)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal pi, alpha;

    pi = 3.141592653589793238;
    alpha = .0078125;
    if (*x < .99999999999999989) {
/* Computing 2nd power */
	d__1 = *mx / *ml;
	ret_val = alpha / pi * (*x * *x - *x * 2. + 2.) / *x * log((1. - *x) *
		 (d__1 * d__1));
    } else {
	ret_val = 0.;
    }
    ret_val /= *mx;
    return ret_val;
} /* llg_ */

