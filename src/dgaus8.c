/* DGAUS8.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c/f2c.h"

/* Table of constant values */

static integer c__14 = 14;
static integer c__5 = 5;
static doublereal c_b6 = 1.;
static doublereal c_b8 = 2.;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__0 = 0;
static logical c_false = FALSE_;
static integer c__72 = 72;
static logical c_true = TRUE_;

/* DECK DGAUS8 */
/* Subroutine */ int dgaus8_(D_fp fun, doublereal *a, doublereal *b, 
	doublereal *err, doublereal *ans, integer *ierr)
{
    /* Initialized data */

    static doublereal x1 = .183434642495649805;
    static integer nlmn = 1;
    static integer kmx = 5000;
    static integer kml = 6;
    static doublereal x2 = .525532409916328986;
    static doublereal x3 = .79666647741362674;
    static doublereal x4 = .960289856497536232;
    static doublereal w1 = .362683783378361983;
    static doublereal w2 = .313706645877887287;
    static doublereal w3 = .222381034453374471;
    static doublereal w4 = .101228536290376259;
    static doublereal sq2 = 1.41421356;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), log(doublereal), pow_di(
	    doublereal *, integer *), sqrt(doublereal);

    /* Local variables */
    static doublereal area, anib;
    static integer nlmx;
    static doublereal c__;
    static integer k, l, nbits;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);
    static doublereal aa[60], ae, ce, ee, ef, hh[60], gl, gr[60];
    static integer lr[60];
    static doublereal vl[60], vr;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nib;
    static doublereal glr;
    static integer lmn;
    static doublereal eps, est, tol;
    static integer lmx, mxl;

/* ***BEGIN PROLOGUE  DGAUS8 */
/* ***PURPOSE  Integrate a real function of one variable over a finite */
/*            interval using an adaptive 8-point Legendre-Gauss */
/*            algorithm.  Intended primarily for high accuracy */
/*            integration or integration of smooth functions. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A1A1 */
/* ***TYPE      DOUBLE PRECISION (GAUS8-S, DGAUS8-D) */
/* ***KEYWORDS  ADAPTIVE QUADRATURE, AUTOMATIC INTEGRATOR, */
/*             GAUSS QUADRATURE, NUMERICAL INTEGRATION */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract  *** a DOUBLE PRECISION routine *** */
/*        DGAUS8 integrates real functions of one variable over finite */
/*        intervals using an adaptive 8-point Legendre-Gauss algorithm. */
/*        DGAUS8 is intended primarily for high accuracy integration */
/*        or integration of smooth functions. */

/*        The maximum number of significant digits obtainable in ANS */
/*        is the smaller of 18 and the number of digits carried in */
/*        double precision arithmetic. */

/*     Description of Arguments */

/*        Input--* FUN, A, B, ERR are DOUBLE PRECISION * */
/*        FUN - name of external function to be integrated.  This name */
/*              must be in an EXTERNAL statement in the calling program. */
/*              FUN must be a DOUBLE PRECISION function of one DOUBLE */
/*              PRECISION argument.  The value of the argument to FUN */
/*              is the variable of integration which ranges from A to B. */
/*        A   - lower limit of integration */
/*        B   - upper limit of integration (may be less than A) */
/*        ERR - is a requested pseudorelative error tolerance.  Normally */
/*              pick a value of ABS(ERR) so that DTOL .LT. ABS(ERR) .LE. */
/*              1.0D-3 where DTOL is the larger of 1.0D-18 and the */
/*              double precision unit roundoff D1MACH(4).  ANS will */
/*              normally have no more error than ABS(ERR) times the */
/*              integral of the absolute value of FUN(X).  Usually, */
/*              smaller values of ERR yield more accuracy and require */
/*              more function evaluations. */

/*              A negative value for ERR causes an estimate of the */
/*              absolute error in ANS to be returned in ERR.  Note that */
/*              ERR must be a variable (not a constant) in this case. */
/*              Note also that the user must reset the value of ERR */
/*              before making any more calls that use the variable ERR. */

/*        Output--* ERR,ANS are double precision * */
/*        ERR - will be an estimate of the absolute error in ANS if the */
/*              input value of ERR was negative.  (ERR is unchanged if */
/*              the input value of ERR was non-negative.)  The estimated */
/*              error is solely for information to the user and should */
/*              not be used as a correction to the computed integral. */
/*        ANS - computed value of integral */
/*        IERR- a status code */
/*            --Normal codes */
/*               1 ANS most likely meets requested error tolerance, */
/*                 or A=B. */
/*              -1 A and B are too nearly equal to allow normal */
/*                 integration.  ANS is set to zero. */
/*            --Abnormal code */
/*               2 ANS probably does not meet requested error tolerance. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, I1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810223  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  DGAUS8 */
/* ***FIRST EXECUTABLE STATEMENT  DGAUS8 */

/*     Initialize */

    k = i1mach_(&c__14);
    anib = d1mach_(&c__5) * k / .30102;
    nbits = (integer) anib;
/* Computing MIN */
    i__1 = 60, i__2 = nbits * 5 / 8;
    nlmx = min(i__1,i__2);
    *ans = 0.;
    *ierr = 1;
    ce = 0.;
    if (*a == *b) {
	goto L140;
    }
    lmx = nlmx;
    lmn = nlmn;
    if (*b == 0.) {
	goto L10;
    }
    if (d_sign(&c_b6, b) * *a <= 0.) {
	goto L10;
    }
    c__ = (d__1 = 1. - *a / *b, abs(d__1));
    if (c__ > .1) {
	goto L10;
    }
    if (c__ <= 0.) {
	goto L140;
    }
    anib = .5 - log(c__) / .69314718;
    nib = (integer) anib;
/* Computing MIN */
    i__1 = nlmx, i__2 = nbits - nib - 7;
    lmx = min(i__1,i__2);
    if (lmx < 1) {
	goto L130;
    }
    lmn = min(lmn,lmx);
L10:
/* Computing MAX */
    i__1 = 5 - nbits;
    d__1 = abs(*err), d__2 = pow_di(&c_b8, &i__1);
    tol = max(d__1,d__2) / 2.;
    if (*err == 0.) {
	tol = sqrt(d1mach_(&c__4));
    }
    eps = tol;
    hh[0] = (*b - *a) / 4.;
    aa[0] = *a;
    lr[0] = 1;
    l = 1;
    d__1 = aa[l - 1] + hh[l - 1] * 2.;
    d__2 = hh[l - 1] * 2.;
    d__3 = d__1 - x1 * d__2;
    d__4 = d__1 + x1 * d__2;
    d__5 = d__1 - x2 * d__2;
    d__6 = d__1 + x2 * d__2;
    d__7 = d__1 - x3 * d__2;
    d__8 = d__1 + x3 * d__2;
    d__9 = d__1 - x4 * d__2;
    d__10 = d__1 + x4 * d__2;
    est = d__2 * (w1 * ((*fun)(&d__3) + (*fun)(&d__4)) + w2 * ((*fun)(&d__5) 
	    + (*fun)(&d__6)) + (w3 * ((*fun)(&d__7) + (*fun)(&d__8)) + w4 * ((
	    *fun)(&d__9) + (*fun)(&d__10))));
    k = 8;
    area = abs(est);
    ef = .5;
    mxl = 0;

/*     Compute refined estimates, estimate the error, etc. */

L20:
    d__1 = aa[l - 1] + hh[l - 1];
    d__2 = d__1 - x1 * hh[l - 1];
    d__3 = d__1 + x1 * hh[l - 1];
    d__4 = d__1 - x2 * hh[l - 1];
    d__5 = d__1 + x2 * hh[l - 1];
    d__6 = d__1 - x3 * hh[l - 1];
    d__7 = d__1 + x3 * hh[l - 1];
    d__8 = d__1 - x4 * hh[l - 1];
    d__9 = d__1 + x4 * hh[l - 1];
    gl = hh[l - 1] * (w1 * ((*fun)(&d__2) + (*fun)(&d__3)) + w2 * ((*fun)(&
	    d__4) + (*fun)(&d__5)) + (w3 * ((*fun)(&d__6) + (*fun)(&d__7)) + 
	    w4 * ((*fun)(&d__8) + (*fun)(&d__9))));
    d__1 = aa[l - 1] + hh[l - 1] * 3.;
    d__2 = d__1 - x1 * hh[l - 1];
    d__3 = d__1 + x1 * hh[l - 1];
    d__4 = d__1 - x2 * hh[l - 1];
    d__5 = d__1 + x2 * hh[l - 1];
    d__6 = d__1 - x3 * hh[l - 1];
    d__7 = d__1 + x3 * hh[l - 1];
    d__8 = d__1 - x4 * hh[l - 1];
    d__9 = d__1 + x4 * hh[l - 1];
    gr[l - 1] = hh[l - 1] * (w1 * ((*fun)(&d__2) + (*fun)(&d__3)) + w2 * ((*
	    fun)(&d__4) + (*fun)(&d__5)) + (w3 * ((*fun)(&d__6) + (*fun)(&
	    d__7)) + w4 * ((*fun)(&d__8) + (*fun)(&d__9))));
    k += 16;
    area += abs(gl) + (d__1 = gr[l - 1], abs(d__1)) - abs(est);
/*     IF (L .LT .LMN) GO TO 11 */
    glr = gl + gr[l - 1];
    ee = (d__1 = est - glr, abs(d__1)) * ef;
/* Computing MAX */
    d__1 = eps * area, d__2 = tol * abs(glr);
    ae = max(d__1,d__2);
    if (ee - ae <= 0.) {
	goto L40;
    } else {
	goto L50;
    }
L30:
    mxl = 1;
L40:
    ce += est - glr;
    if (lr[l - 1] <= 0) {
	goto L60;
    } else {
	goto L80;
    }

/*     Consider the left half of this level */

L50:
    if (k > kmx) {
	lmx = kml;
    }
    if (l >= lmx) {
	goto L30;
    }
    ++l;
    eps *= .5;
    ef /= sq2;
    hh[l - 1] = hh[l - 2] * .5;
    lr[l - 1] = -1;
    aa[l - 1] = aa[l - 2];
    est = gl;
    goto L20;

/*     Proceed to right half at this level */

L60:
    vl[l - 1] = glr;
L70:
    est = gr[l - 2];
    lr[l - 1] = 1;
    aa[l - 1] += hh[l - 1] * 4.;
    goto L20;

/*     Return one level */

L80:
    vr = glr;
L90:
    if (l <= 1) {
	goto L120;
    }
    --l;
    eps *= 2.;
    ef *= sq2;
    if (lr[l - 1] <= 0) {
	goto L100;
    } else {
	goto L110;
    }
L100:
    vl[l - 1] = vl[l] + vr;
    goto L70;
L110:
    vr = vl[l] + vr;
    goto L90;

/*     Exit */

L120:
    *ans = vr;
    if (mxl == 0 || abs(ce) <= tol * 2. * area) {
	goto L140;
    }
    *ierr = 2;
    xermsg_("SLATEC", "DGAUS8", "ANS is probably insufficiently accurate.", &
	    c__3, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)40);
    goto L140;
L130:
    *ierr = -1;
    xermsg_("SLATEC", "DGAUS8", "A and B are too nearly equal to allow norma\
l integration. $$ANS is set to zero and IERR to -1.", &c__1, &c_n1, (ftnlen)6,
	     (ftnlen)6, (ftnlen)94);
L140:
    if (*err < 0.) {
	*err = ce;
    }
    return 0;
} /* dgaus8_ */

/* DECK I1MACH */
integer i1mach_(integer *i__)
{
    /* Initialized data */

    static integer imach[16] = { 5,6,6,6,32,4,2,31,2147483647,2,24,-126,127,
	    53,-1022,1023 };

    /* Format strings */
    static char fmt_9000[] = "(\0021ERROR    1 IN I1MACH - I OUT OF BOUND\
S\002)";

    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer output;

    /* Fortran I/O blocks */
    static cilist io___41 = { 0, 0, 0, fmt_9000, 0 };


/* ***BEGIN PROLOGUE  I1MACH */
/* ***PURPOSE  Return integer machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      INTEGER (I1MACH-I) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   I1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument and can be referenced as follows: */

/*        K = I1MACH(I) */

/*   where I=1,...,16.  The (output) value of K above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   I/O unit numbers: */
/*     I1MACH( 1) = the standard input unit. */
/*     I1MACH( 2) = the standard output unit. */
/*     I1MACH( 3) = the standard punch unit. */
/*     I1MACH( 4) = the standard error message unit. */

/*   Words: */
/*     I1MACH( 5) = the number of bits per integer storage unit. */
/*     I1MACH( 6) = the number of characters per integer storage unit. */

/*   Integers: */
/*     assume integers are represented in the S-digit, base-A form */

/*                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) */

/*                where 0 .LE. X(I) .LT. A for I=0,...,S-1. */
/*     I1MACH( 7) = A, the base. */
/*     I1MACH( 8) = S, the number of base-A digits. */
/*     I1MACH( 9) = A**S - 1, the largest magnitude. */

/*   Floating-Point Numbers: */
/*     Assume floating-point numbers are represented in the T-digit, */
/*     base-B form */
/*                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*                where 0 .LE. X(I) .LT. B for I=1,...,T, */
/*                0 .LT. X(1), and EMIN .LE. E .LE. EMAX. */
/*     I1MACH(10) = B, the base. */

/*   Single-Precision: */
/*     I1MACH(11) = T, the number of base-B digits. */
/*     I1MACH(12) = EMIN, the smallest exponent E. */
/*     I1MACH(13) = EMAX, the largest exponent E. */

/*   Double-Precision: */
/*     I1MACH(14) = T, the number of base-B digits. */
/*     I1MACH(15) = EMIN, the smallest exponent E. */
/*     I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   891012  Added VAX G-floating constants.  (WRB) */
/*   891012  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16. */
/*           (RWC) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added Convex -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/*   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler */
/*           options.  (DWL, RWC and WRB). */
/*   010817  Elevated IEEE to highest importance; see next set of */
/*           comments below.  (DWL) */
/* ***END PROLOGUE  I1MACH */

/* Initial data here correspond to the IEEE standard.  If one of the */
/* sets of initial data below is preferred, do the necessary commenting */
/* and uncommenting. (DWL) */
/* $$$      EQUIVALENCE (IMACH(4),OUTPUT) */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        129 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1025 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA IMACH( 1) /          7 / */
/*     DATA IMACH( 2) /          2 / */
/*     DATA IMACH( 3) /          2 / */
/*     DATA IMACH( 4) /          2 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         33 / */
/*     DATA IMACH( 9) / Z1FFFFFFFF / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -256 / */
/*     DATA IMACH(13) /        255 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /       -256 / */
/*     DATA IMACH(16) /        255 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /        -50 / */
/*     DATA IMACH(16) /         76 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /     -32754 / */
/*     DATA IMACH(16) /      32780 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -4095 / */
/*     DATA IMACH(13) /       4094 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -4095 / */
/*     DATA IMACH(16) /       4094 / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /    6LOUTPUT/ */
/*     DATA IMACH( 5) /         60 / */
/*     DATA IMACH( 6) /         10 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         48 / */
/*     DATA IMACH( 9) / 00007777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /       -929 / */
/*     DATA IMACH(13) /       1070 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /       -929 / */
/*     DATA IMACH(16) /       1069 / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / Z'7FFFFFFF' / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1023 / */
/*     DATA IMACH(13) /       1023 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16383 / */
/*     DATA IMACH(16) /      16383 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -pd8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1023 / */
/*     DATA IMACH(13) /       1023 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 46 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         46 / */
/*     DATA IMACH( 9) / 1777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 64 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 777777777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */

/*     DATA IMACH( 1) /         11 / */
/*     DATA IMACH( 2) /         12 / */
/*     DATA IMACH( 3) /          8 / */
/*     DATA IMACH( 4) /         10 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         24 / */
/*     DATA IMACH( 6) /          3 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         23 / */
/*     DATA IMACH( 9) /    8388607 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         38 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /         43 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         63 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         39 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         55 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          7 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1015 / */
/*     DATA IMACH(16) /       1017 / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) /  Z7FFFFFFF / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE IBM PC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         54 / */
/*     DATA IMACH(15) /       -101 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         62 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1021 / */
/*     DATA IMACH(13) /       1024 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16381 / */
/*     DATA IMACH(16) /      16384 / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          1 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /      -1024 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR */

/*     DATA IMACH( 1) /          1 / */
/*     DATA IMACH( 2) /          1 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/* ***FIRST EXECUTABLE STATEMENT  I1MACH */
    if (*i__ < 1 || *i__ > 16) {
	goto L10;
    }

    ret_val = imach[*i__ - 1];
    return ret_val;

L10:
    io___41.ciunit = output;
    s_wsfe(&io___41);
    e_wsfe();

/*     CALL FDUMP */

    s_stop("", (ftnlen)0);
    return ret_val;
} /* i1mach_ */

/* DECK XERMSG */
/* Subroutine */ int xermsg_(char *librar, char *subrou, char *messg, integer 
	*nerr, integer *level, ftnlen librar_len, ftnlen subrou_len, ftnlen 
	messg_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    char ch__1[87];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_len(char *, ftnlen), s_wsfi(icilist *), do_fio(integer *, char *
	    , ftnlen), e_wsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer lerr;
    static char temp[72];
    static integer i__;
    extern /* Subroutine */ int fdump_();
    static char xlibr[8];
    static integer ltemp, kount;
    static char xsubr[8];
    extern integer j4save_(integer *, integer *, logical *);
    static integer llevel, maxmes;
    static char lfirst[20];
    extern /* Subroutine */ int xercnt_(char *, char *, char *, integer *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer lkntrl, kdummy;
    extern /* Subroutine */ int xerhlt_(char *, ftnlen);
    static integer mkntrl;
    extern /* Subroutine */ int xersve_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen), xerprn_(
	    char *, integer *, char *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___55 = { 0, temp, 0, "('ERROR NUMBER = ', I8)", 72, 1 };


/* ***BEGIN PROLOGUE  XERMSG */
/* ***PURPOSE  Process error messages for SLATEC and other libraries. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERMSG-A) */
/* ***KEYWORDS  ERROR MESSAGE, XERROR */
/* ***AUTHOR  Fong, Kirby, (NMFECC at LLNL) */
/* ***DESCRIPTION */

/*   XERMSG processes a diagnostic message in a manner determined by the */
/*   value of LEVEL and the current value of the library error control */
/*   flag, KONTRL.  See subroutine XSETF for details. */

/*    LIBRAR   A character constant (or character variable) with the name */
/*             of the library.  This will be 'SLATEC' for the SLATEC */
/*             Common Math Library.  The error handling package is */
/*             general enough to be used by many libraries */
/*             simultaneously, so it is desirable for the routine that */
/*             detects and reports an error to identify the library name */
/*             as well as the routine name. */

/*    SUBROU   A character constant (or character variable) with the name */
/*             of the routine that detected the error.  Usually it is the */
/*             name of the routine that is calling XERMSG.  There are */
/*             some instances where a user callable library routine calls */
/*             lower level subsidiary routines where the error is */
/*             detected.  In such cases it may be more informative to */
/*             supply the name of the routine the user called rather than */
/*             the name of the subsidiary routine that detected the */
/*             error. */

/*    MESSG    A character constant (or character variable) with the text */
/*             of the error or warning message.  In the example below, */
/*             the message is a character constant that contains a */
/*             generic message. */

/*                   CALL XERMSG ('SLATEC', 'MMPY', */
/*                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION', */
/*                  *3, 1) */

/*             It is possible (and is sometimes desirable) to generate a */
/*             specific message--e.g., one that contains actual numeric */
/*             values.  Specific numeric values can be converted into */
/*             character strings using formatted WRITE statements into */
/*             character variables.  This is called standard Fortran */
/*             internal file I/O and is exemplified in the first three */
/*             lines of the following example.  You can also catenate */
/*             substrings of characters to construct the error message. */
/*             Here is an example showing the use of both writing to */
/*             an internal file and catenating character strings. */

/*                   CHARACTER*5 CHARN, CHARL */
/*                   WRITE (CHARN,10) N */
/*                   WRITE (CHARL,10) LDA */
/*                10 FORMAT(I5) */
/*                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN// */
/*                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'// */
/*                  *   CHARL, 3, 1) */

/*             There are two subtleties worth mentioning.  One is that */
/*             the // for character catenation is used to construct the */
/*             error message so that no single character constant is */
/*             continued to the next line.  This avoids confusion as to */
/*             whether there are trailing blanks at the end of the line. */
/*             The second is that by catenating the parts of the message */
/*             as an actual argument rather than encoding the entire */
/*             message into one large character variable, we avoid */
/*             having to know how long the message will be in order to */
/*             declare an adequate length for that large character */
/*             variable.  XERMSG calls XERPRN to print the message using */
/*             multiple lines if necessary.  If the message is very long, */
/*             XERPRN will break it into pieces of 72 characters (as */
/*             requested by XERMSG) for printing on multiple lines. */
/*             Also, XERMSG asks XERPRN to prefix each line with ' *  ' */
/*             so that the total line length could be 76 characters. */
/*             Note also that XERPRN scans the error message backwards */
/*             to ignore trailing blanks.  Another feature is that */
/*             the substring '$$' is treated as a new line sentinel */
/*             by XERPRN.  If you want to construct a multiline */
/*             message without having to count out multiples of 72 */
/*             characters, just use '$$' as a separator.  '$$' */
/*             obviously must occur within 72 characters of the */
/*             start of each line to have its intended effect since */
/*             XERPRN is asked to wrap around at 72 characters in */
/*             addition to looking for '$$'. */

/*    NERR     An integer value that is chosen by the library routine's */
/*             author.  It must be in the range -99 to 999 (three */
/*             printable digits).  Each distinct error should have its */
/*             own error number.  These error numbers should be described */
/*             in the machine readable documentation for the routine. */
/*             The error numbers need be unique only within each routine, */
/*             so it is reasonable for each routine to start enumerating */
/*             errors from 1 and proceeding to the next integer. */

/*    LEVEL    An integer value in the range 0 to 2 that indicates the */
/*             level (severity) of the error.  Their meanings are */

/*            -1  A warning message.  This is used if it is not clear */
/*                that there really is an error, but the user's attention */
/*                may be needed.  An attempt is made to only print this */
/*                message once. */

/*             0  A warning message.  This is used if it is not clear */
/*                that there really is an error, but the user's attention */
/*                may be needed. */

/*             1  A recoverable error.  This is used even if the error is */
/*                so serious that the routine cannot return any useful */
/*                answer.  If the user has told the error package to */
/*                return after recoverable errors, then XERMSG will */
/*                return to the Library routine which can then return to */
/*                the user's routine.  The user may also permit the error */
/*                package to terminate the program upon encountering a */
/*                recoverable error. */

/*             2  A fatal error.  XERMSG will not return to its caller */
/*                after it receives a fatal error.  This level should */
/*                hardly ever be used; it is much better to allow the */
/*                user a chance to recover.  An example of one of the few */
/*                cases in which it is permissible to declare a level 2 */
/*                error is a reverse communication Library routine that */
/*                is likely to be called repeatedly until it integrates */
/*                across some interval.  If there is a serious error in */
/*                the input such that another step cannot be taken and */
/*                the Library routine is called again without the input */
/*                error having been corrected by the caller, the Library */
/*                routine will probably be called forever with improper */
/*                input.  In this case, it is reasonable to declare the */
/*                error to be fatal. */

/*    Each of the arguments to XERMSG is input; none will be modified by */
/*    XERMSG.  A routine may make multiple calls to XERMSG with warning */
/*    level messages; however, after a call to XERMSG with a recoverable */
/*    error, the routine should return to the user.  Do not try to call */
/*    XERMSG with a second recoverable error after the first recoverable */
/*    error because the error package saves the error number.  The user */
/*    can retrieve this error number by calling another entry point in */
/*    the error handling package and then clear the error number when */
/*    recovering from the error.  Calling XERMSG in succession causes the */
/*    old error number to be overwritten by the latest error number. */
/*    This is considered harmless for error numbers associated with */
/*    warning messages but must not be done for error numbers of serious */
/*    errors.  After a call to XERMSG with a recoverable error, the user */
/*    must be given a chance to call NUMXER or XERCLR to retrieve or */
/*    clear the error number. */
/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   880101  DATE WRITTEN */
/*   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988. */
/*           THERE ARE TWO BASIC CHANGES. */
/*           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO */
/*               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES */
/*               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS */
/*               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE */
/*               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER */
/*               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY */
/*               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE */
/*               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76. */
/*           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE */
/*               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE */
/*               OF LOWER CASE. */
/*   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30. */
/*           THE PRINCIPAL CHANGES ARE */
/*           1.  CLARIFY COMMENTS IN THE PROLOGUES */
/*           2.  RENAME XRPRNT TO XERPRN */
/*           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES */
/*               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE / */
/*               CHARACTER FOR NEW RECORDS. */
/*   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO */
/*           CLEAN UP THE CODING. */
/*   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN */
/*           PREFIX. */
/*   891013  REVISED TO CORRECT COMMENTS. */
/*   891214  Prologue converted to Version 4.0 format.  (WRB) */
/*   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but */
/*           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added */
/*           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and */
/*           XERCTL to XERCNT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERMSG */
/* ***FIRST EXECUTABLE STATEMENT  XERMSG */
    lkntrl = j4save_(&c__2, &c__0, &c_false);
    maxmes = j4save_(&c__4, &c__0, &c_false);

/*       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL. */
/*       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE */
/*          SHOULD BE PRINTED. */

/*       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN */
/*          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE, */
/*          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2. */

    if (*nerr < -9999999 || *nerr > 99999999 || *nerr == 0 || *level < -1 || *
	    level > 2) {
	xerprn_(" ***", &c_n1, "FATAL ERROR IN...$$ XERMSG -- INVALID ERROR \
NUMBER OR LEVEL$$ JOB ABORT DUE TO FATAL ERROR.", &c__72, (ftnlen)4, (ftnlen)
		91);
	xersve_(" ", " ", " ", &c__0, &c__0, &c__0, &kdummy, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
	xerhlt_(" ***XERMSG -- INVALID INPUT", (ftnlen)27);
	return 0;
    }

/*       RECORD THE MESSAGE. */

    i__ = j4save_(&c__1, nerr, &c_true);
    xersve_(librar, subrou, messg, &c__1, nerr, level, &kount, librar_len, 
	    subrou_len, messg_len);

/*       HANDLE PRINT-ONCE WARNING MESSAGES. */

    if (*level == -1 && kount > 1) {
	return 0;
    }

/*       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG. */

    s_copy(xlibr, librar, (ftnlen)8, librar_len);
    s_copy(xsubr, subrou, (ftnlen)8, subrou_len);
    s_copy(lfirst, messg, (ftnlen)20, messg_len);
    lerr = *nerr;
    llevel = *level;
    xercnt_(xlibr, xsubr, lfirst, &lerr, &llevel, &lkntrl, (ftnlen)8, (ftnlen)
	    8, (ftnlen)20);

/* Computing MAX */
    i__1 = -2, i__2 = min(2,lkntrl);
    lkntrl = max(i__1,i__2);
    mkntrl = abs(lkntrl);

/*       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS */
/*       ZERO AND THE ERROR IS NOT FATAL. */

    if (*level < 2 && lkntrl == 0) {
	goto L30;
    }
    if (*level == 0 && kount > maxmes) {
	goto L30;
    }
    if (*level == 1 && kount > maxmes && mkntrl == 1) {
	goto L30;
    }
    if (*level == 2 && kount > max(1,maxmes)) {
	goto L30;
    }

/*       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A */
/*       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS) */
/*       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG */
/*       IS NOT ZERO. */

    if (lkntrl != 0) {
	s_copy(temp, "MESSAGE FROM ROUTINE ", (ftnlen)21, (ftnlen)21);
/* Computing MIN */
	i__1 = i_len(subrou, subrou_len);
	i__ = min(i__1,16);
	s_copy(temp + 21, subrou, i__, i__);
	i__1 = i__ + 21;
	s_copy(temp + i__1, " IN LIBRARY ", i__ + 33 - i__1, (ftnlen)12);
	ltemp = i__ + 33;
/* Computing MIN */
	i__1 = i_len(librar, librar_len);
	i__ = min(i__1,16);
	i__1 = ltemp;
	s_copy(temp + i__1, librar, ltemp + i__ - i__1, i__);
	i__1 = ltemp + i__;
	s_copy(temp + i__1, ".", ltemp + i__ + 1 - i__1, (ftnlen)1);
	ltemp = ltemp + i__ + 1;
	xerprn_(" ***", &c_n1, temp, &c__72, (ftnlen)4, ltemp);
    }

/*       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE */
/*       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE */
/*       FROM EACH OF THE FOLLOWING THREE OPTIONS. */
/*       1.  LEVEL OF THE MESSAGE */
/*              'INFORMATIVE MESSAGE' */
/*              'POTENTIALLY RECOVERABLE ERROR' */
/*              'FATAL ERROR' */
/*       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE */
/*              'PROG CONTINUES' */
/*              'PROG ABORTED' */
/*       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK */
/*           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS */
/*           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.) */
/*              'TRACEBACK REQUESTED' */
/*              'TRACEBACK NOT REQUESTED' */
/*       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT */
/*       EXCEED 74 CHARACTERS. */
/*       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED. */

    if (lkntrl > 0) {

/*       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL. */

	if (*level <= 0) {
	    s_copy(temp, "INFORMATIVE MESSAGE,", (ftnlen)20, (ftnlen)20);
	    ltemp = 20;
	} else if (*level == 1) {
	    s_copy(temp, "POTENTIALLY RECOVERABLE ERROR,", (ftnlen)30, (
		    ftnlen)30);
	    ltemp = 30;
	} else {
	    s_copy(temp, "FATAL ERROR,", (ftnlen)12, (ftnlen)12);
	    ltemp = 12;
	}

/*       THEN WHETHER THE PROGRAM WILL CONTINUE. */

	if (mkntrl == 2 && *level >= 1 || mkntrl == 1 && *level == 2) {
	    i__1 = ltemp;
	    s_copy(temp + i__1, " PROG ABORTED,", ltemp + 14 - i__1, (ftnlen)
		    14);
	    ltemp += 14;
	} else {
	    i__1 = ltemp;
	    s_copy(temp + i__1, " PROG CONTINUES,", ltemp + 16 - i__1, (
		    ftnlen)16);
	    ltemp += 16;
	}

/*       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK. */

	if (lkntrl > 0) {
	    i__1 = ltemp;
	    s_copy(temp + i__1, " TRACEBACK REQUESTED", ltemp + 20 - i__1, (
		    ftnlen)20);
	    ltemp += 20;
	} else {
	    i__1 = ltemp;
	    s_copy(temp + i__1, " TRACEBACK NOT REQUESTED", ltemp + 24 - i__1,
		     (ftnlen)24);
	    ltemp += 24;
	}
	xerprn_(" ***", &c_n1, temp, &c__72, (ftnlen)4, ltemp);
    }

/*       NOW SEND OUT THE MESSAGE. */

    xerprn_(" *  ", &c_n1, messg, &c__72, (ftnlen)4, messg_len);

/*       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A */
/*          TRACEBACK. */

    if (lkntrl > 0) {
	s_wsfi(&io___55);
	do_fio(&c__1, (char *)&(*nerr), (ftnlen)sizeof(integer));
	e_wsfi();
	for (i__ = 16; i__ <= 22; ++i__) {
	    if (*(unsigned char *)&temp[i__ - 1] != ' ') {
		goto L20;
	    }
/* L10: */
	}

L20:
/* Writing concatenation */
	i__3[0] = 15, a__1[0] = temp;
	i__3[1] = 23 - (i__ - 1), a__1[1] = temp + (i__ - 1);
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)87);
	xerprn_(" *  ", &c_n1, ch__1, &c__72, (ftnlen)4, 23 - (i__ - 1) + 15);
	fdump_();
    }

/*       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE. */

    if (lkntrl != 0) {
	xerprn_(" *  ", &c_n1, " ", &c__72, (ftnlen)4, (ftnlen)1);
	xerprn_(" ***", &c_n1, "END OF MESSAGE", &c__72, (ftnlen)4, (ftnlen)
		14);
	xerprn_("    ", &c__0, " ", &c__72, (ftnlen)4, (ftnlen)1);
    }

/*       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE */
/*       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN. */

L30:
    if (*level <= 0 || *level == 1 && mkntrl <= 1) {
	return 0;
    }

/*       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A */
/*       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR */
/*       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT. */

    if (lkntrl > 0 && kount < max(1,maxmes)) {
	if (*level == 1) {
	    xerprn_(" ***", &c_n1, "JOB ABORT DUE TO UNRECOVERED ERROR.", &
		    c__72, (ftnlen)4, (ftnlen)35);
	} else {
	    xerprn_(" ***", &c_n1, "JOB ABORT DUE TO FATAL ERROR.", &c__72, (
		    ftnlen)4, (ftnlen)29);
	}
	xersve_(" ", " ", " ", &c_n1, &c__0, &c__0, &kdummy, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
	xerhlt_(" ", (ftnlen)1);
    } else {
	xerhlt_(messg, messg_len);
    }
    return 0;
} /* xermsg_ */

/* DECK XERPRN */
/* Subroutine */ int xerprn_(char *prefix, integer *npref, char *messg, 
	integer *nwrap, ftnlen prefix_len, ftnlen messg_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    i_indx(char *, char *, ftnlen, ftnlen), s_cmp(char *, char *, 
	    ftnlen, ftnlen);

    /* Local variables */
    static integer i__, n;
    static char cbuff[148];
    static integer lpref, nextc, lwrap, nunit;
    extern integer i1mach_(integer *);
    static integer iu[5], lpiece, idelta, lenmsg;
    extern /* Subroutine */ int xgetua_(integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___64 = { 0, 0, 0, "(A)", 0 };
    static cilist io___68 = { 0, 0, 0, "(A)", 0 };


/* ***BEGIN PROLOGUE  XERPRN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Print error messages processed by XERMSG. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERPRN-A) */
/* ***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR */
/* ***AUTHOR  Fong, Kirby, (NMFECC at LLNL) */
/* ***DESCRIPTION */

/* This routine sends one or more lines to each of the (up to five) */
/* logical units to which error messages are to be sent.  This routine */
/* is called several times by XERMSG, sometimes with a single line to */
/* print and sometimes with a (potentially very long) message that may */
/* wrap around into multiple lines. */

/* PREFIX  Input argument of type CHARACTER.  This argument contains */
/*         characters to be put at the beginning of each line before */
/*         the body of the message.  No more than 16 characters of */
/*         PREFIX will be used. */

/* NPREF   Input argument of type INTEGER.  This argument is the number */
/*         of characters to use from PREFIX.  If it is negative, the */
/*         intrinsic function LEN is used to determine its length.  If */
/*         it is zero, PREFIX is not used.  If it exceeds 16 or if */
/*         LEN(PREFIX) exceeds 16, only the first 16 characters will be */
/*         used.  If NPREF is positive and the length of PREFIX is less */
/*         than NPREF, a copy of PREFIX extended with blanks to length */
/*         NPREF will be used. */

/* MESSG   Input argument of type CHARACTER.  This is the text of a */
/*         message to be printed.  If it is a long message, it will be */
/*         broken into pieces for printing on multiple lines.  Each line */
/*         will start with the appropriate prefix and be followed by a */
/*         piece of the message.  NWRAP is the number of characters per */
/*         piece; that is, after each NWRAP characters, we break and */
/*         start a new line.  In addition the characters '$$' embedded */
/*         in MESSG are a sentinel for a new line.  The counting of */
/*         characters up to NWRAP starts over for each new line.  The */
/*         value of NWRAP typically used by XERMSG is 72 since many */
/*         older error messages in the SLATEC Library are laid out to */
/*         rely on wrap-around every 72 characters. */

/* NWRAP   Input argument of type INTEGER.  This gives the maximum size */
/*         piece into which to break MESSG for printing on multiple */
/*         lines.  An embedded '$$' ends a line, and the count restarts */
/*         at the following character.  If a line break does not occur */
/*         on a blank (it would split a word) that word is moved to the */
/*         next line.  Values of NWRAP less than 16 will be treated as */
/*         16.  Values of NWRAP greater than 132 will be treated as 132. */
/*         The actual line length will be NPREF + NWRAP after NPREF has */
/*         been adjusted to fall between 0 and 16 and NWRAP has been */
/*         adjusted to fall between 16 and 132. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  I1MACH, XGETUA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   880621  DATE WRITTEN */
/*   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF */
/*           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK */
/*           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE */
/*           SLASH CHARACTER IN FORMAT STATEMENTS. */
/*   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO */
/*           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK */
/*           LINES TO BE PRINTED. */
/*   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF */
/*           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH. */
/*   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH. */
/*   891214  Prologue converted to Version 4.0 format.  (WRB) */
/*   900510  Added code to break messages between words.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERPRN */
/* ***FIRST EXECUTABLE STATEMENT  XERPRN */
    xgetua_(iu, &nunit);

/*       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD */
/*       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD */
/*       ERROR MESSAGE UNIT. */

    n = i1mach_(&c__4);
    i__1 = nunit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iu[i__ - 1] == 0) {
	    iu[i__ - 1] = n;
	}
/* L10: */
    }

/*       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE */
/*       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING */
/*       THE REST OF THIS ROUTINE. */

    if (*npref < 0) {
	lpref = i_len(prefix, prefix_len);
    } else {
	lpref = *npref;
    }
    lpref = min(16,lpref);
    if (lpref != 0) {
	s_copy(cbuff, prefix, lpref, prefix_len);
    }

/*       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE */
/*       TIME FROM MESSG TO PRINT ON ONE LINE. */

/* Computing MAX */
    i__1 = 16, i__2 = min(132,*nwrap);
    lwrap = max(i__1,i__2);

/*       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS. */

    lenmsg = i_len(messg, messg_len);
    n = lenmsg;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&messg[lenmsg - 1] != ' ') {
	    goto L30;
	}
	--lenmsg;
/* L20: */
    }
L30:

/*       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE. */

    if (lenmsg == 0) {
	i__1 = lpref;
	s_copy(cbuff + i__1, " ", lpref + 1 - i__1, (ftnlen)1);
	i__1 = nunit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___64.ciunit = iu[i__ - 1];
	    s_wsfe(&io___64);
	    do_fio(&c__1, cbuff, lpref + 1);
	    e_wsfe();
/* L40: */
	}
	return 0;
    }

/*       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING */
/*       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL. */
/*       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT. */
/*       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED. */

/*       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE */
/*       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE */
/*       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH */
/*       OF THE SECOND ARGUMENT. */

/*       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE */
/*       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER */
/*       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT */
/*       POSITION NEXTC. */

/*       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE */
/*                       REMAINDER OF THE CHARACTER STRING.  LPIECE */
/*                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC, */
/*                       WHICHEVER IS LESS. */

/*       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC: */
/*                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE */
/*                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY */
/*                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION */
/*                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF */
/*                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE */
/*                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC */
/*                       SHOULD BE INCREMENTED BY 2. */

/*       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP. */

/*       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1 */
/*                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS */
/*                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ. */
/*                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY */
/*                       AT THE END OF A LINE. */

    nextc = 1;
L50:
    lpiece = i_indx(messg + (nextc - 1), "$$", lenmsg - (nextc - 1), (ftnlen)
	    2);
    if (lpiece == 0) {

/*       THERE WAS NO NEW LINE SENTINEL FOUND. */

	idelta = 0;
/* Computing MIN */
	i__1 = lwrap, i__2 = lenmsg + 1 - nextc;
	lpiece = min(i__1,i__2);
	if (lpiece < lenmsg + 1 - nextc) {
	    for (i__ = lpiece + 1; i__ >= 2; --i__) {
		i__1 = nextc + i__ - 2;
		if (s_cmp(messg + i__1, " ", nextc + i__ - 1 - i__1, (ftnlen)
			1) == 0) {
		    lpiece = i__ - 1;
		    idelta = 1;
		    goto L54;
		}
/* L52: */
	    }
	}
L54:
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
	nextc = nextc + lpiece + idelta;
    } else if (lpiece == 1) {

/*       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1). */
/*       DON'T PRINT A BLANK LINE. */

	nextc += 2;
	goto L50;
    } else if (lpiece > lwrap + 1) {

/*       LPIECE SHOULD BE SET DOWN TO LWRAP. */

	idelta = 0;
	lpiece = lwrap;
	for (i__ = lpiece + 1; i__ >= 2; --i__) {
	    i__1 = nextc + i__ - 2;
	    if (s_cmp(messg + i__1, " ", nextc + i__ - 1 - i__1, (ftnlen)1) ==
		     0) {
		lpiece = i__ - 1;
		idelta = 1;
		goto L58;
	    }
/* L56: */
	}
L58:
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
	nextc = nextc + lpiece + idelta;
    } else {

/*       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1. */
/*       WE SHOULD DECREMENT LPIECE BY ONE. */

	--lpiece;
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
	nextc = nextc + lpiece + 2;
    }

/*       PRINT */

    i__1 = nunit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___68.ciunit = iu[i__ - 1];
	s_wsfe(&io___68);
	do_fio(&c__1, cbuff, lpref + lpiece);
	e_wsfe();
/* L60: */
    }

    if (nextc <= lenmsg) {
	goto L50;
    }
    return 0;
} /* xerprn_ */

/* DECK XERSVE */
/* Subroutine */ int xersve_(char *librar, char *subrou, char *messg, integer 
	*kflag, integer *nerr, integer *level, integer *icount, ftnlen 
	librar_len, ftnlen subrou_len, ftnlen messg_len)
{
    /* Initialized data */

    static integer kountx = 0;
    static integer nmsg = 0;

    /* Format strings */
    static char fmt_9000[] = "(\0020          ERROR MESSAGE SUMMARY\002/\002\
 LIBRARY    SUBROUTINE MESSAGE START             NERR\002,\002     LEVEL    \
 COUNT\002)";
    static char fmt_9010[] = "(1x,a,3x,a,3x,a,3i10)";
    static char fmt_9020[] = "(\0020OTHER ERRORS NOT INDIVIDUALLY TABULATED \
= \002,i10)";
    static char fmt_9030[] = "(1x)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, iunit, kunit, nunit, kount[10];
    extern integer i1mach_(integer *);
    static char libtab[8*10], mestab[20*10];
    static integer nertab[10], levtab[10];
    static char subtab[8*10];
    extern /* Subroutine */ int xgetua_(integer *, integer *);
    static char lib[8], mes[20], sub[8];
    static integer lun[5];

    /* Fortran I/O blocks */
    static cilist io___75 = { 0, 0, 0, fmt_9000, 0 };
    static cilist io___77 = { 0, 0, 0, fmt_9010, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_9030, 0 };


/* ***BEGIN PROLOGUE  XERSVE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Record that an error has occurred. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3 */
/* ***TYPE      ALL (XERSVE-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/* *Usage: */

/*        INTEGER  KFLAG, NERR, LEVEL, ICOUNT */
/*        CHARACTER * (len) LIBRAR, SUBROU, MESSG */

/*        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT) */

/* *Arguments: */

/*        LIBRAR :IN    is the library that the message is from. */
/*        SUBROU :IN    is the subroutine that the message is from. */
/*        MESSG  :IN    is the message to be saved. */
/*        KFLAG  :IN    indicates the action to be performed. */
/*                      when KFLAG > 0, the message in MESSG is saved. */
/*                      when KFLAG=0 the tables will be dumped and */
/*                      cleared. */
/*                      when KFLAG < 0, the tables will be dumped and */
/*                      not cleared. */
/*        NERR   :IN    is the error number. */
/*        LEVEL  :IN    is the error severity. */
/*        ICOUNT :OUT   the number of times this message has been seen, */
/*                      or zero if the table has overflowed and does not */
/*                      contain this message specifically.  When KFLAG=0, */
/*                      ICOUNT will not be altered. */

/* *Description: */

/*   Record that this error occurred and possibly dump and clear the */
/*   tables. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  I1MACH, XGETUA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800319  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900413  Routine modified to remove reference to KFLAG.  (WRB) */
/*   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling */
/*           sequence, use IF-THEN-ELSE, make number of saved entries */
/*           easily changeable, changed routine name from XERSAV to */
/*           XERSVE.  (RWC) */
/*   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERSVE */
/* ***FIRST EXECUTABLE STATEMENT  XERSVE */

    if (*kflag <= 0) {

/*        Dump the table. */

	if (nmsg == 0) {
	    return 0;
	}

/*        Print to each unit. */

	xgetua_(lun, &nunit);
	i__1 = nunit;
	for (kunit = 1; kunit <= i__1; ++kunit) {
	    iunit = lun[kunit - 1];
	    if (iunit == 0) {
		iunit = i1mach_(&c__4);
	    }

/*           Print the table header. */

	    io___75.ciunit = iunit;
	    s_wsfe(&io___75);
	    e_wsfe();

/*           Print body of table. */

	    i__2 = nmsg;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		io___77.ciunit = iunit;
		s_wsfe(&io___77);
		do_fio(&c__1, libtab + (i__ - 1 << 3), (ftnlen)8);
		do_fio(&c__1, subtab + (i__ - 1 << 3), (ftnlen)8);
		do_fio(&c__1, mestab + (i__ - 1) * 20, (ftnlen)20);
		do_fio(&c__1, (char *)&nertab[i__ - 1], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&levtab[i__ - 1], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&kount[i__ - 1], (ftnlen)sizeof(integer)
			);
		e_wsfe();
/* L10: */
	    }

/*           Print number of other errors. */

	    if (kountx != 0) {
		io___84.ciunit = iunit;
		s_wsfe(&io___84);
		do_fio(&c__1, (char *)&kountx, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    io___85.ciunit = iunit;
	    s_wsfe(&io___85);
	    e_wsfe();
/* L20: */
	}

/*        Clear the error tables. */

	if (*kflag == 0) {
	    nmsg = 0;
	    kountx = 0;
	}
    } else {

/*        PROCESS A MESSAGE... */
/*        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG, */
/*        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL. */

	s_copy(lib, librar, (ftnlen)8, librar_len);
	s_copy(sub, subrou, (ftnlen)8, subrou_len);
	s_copy(mes, messg, (ftnlen)20, messg_len);
	i__1 = nmsg;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (s_cmp(lib, libtab + (i__ - 1 << 3), (ftnlen)8, (ftnlen)8) == 
		    0 && s_cmp(sub, subtab + (i__ - 1 << 3), (ftnlen)8, (
		    ftnlen)8) == 0 && s_cmp(mes, mestab + (i__ - 1) * 20, (
		    ftnlen)20, (ftnlen)20) == 0 && *nerr == nertab[i__ - 1] &&
		     *level == levtab[i__ - 1]) {
		++kount[i__ - 1];
		*icount = kount[i__ - 1];
		return 0;
	    }
/* L30: */
	}

	if (nmsg < 10) {

/*           Empty slot found for new message. */

	    ++nmsg;
	    s_copy(libtab + (i__ - 1 << 3), lib, (ftnlen)8, (ftnlen)8);
	    s_copy(subtab + (i__ - 1 << 3), sub, (ftnlen)8, (ftnlen)8);
	    s_copy(mestab + (i__ - 1) * 20, mes, (ftnlen)20, (ftnlen)20);
	    nertab[i__ - 1] = *nerr;
	    levtab[i__ - 1] = *level;
	    kount[i__ - 1] = 1;
	    *icount = 1;
	} else {

/*           Table is full. */

	    ++kountx;
	    *icount = 0;
	}
    }
    return 0;

/*     Formats. */

} /* xersve_ */

/* DECK D1MACH */
doublereal d1mach_(integer *i__)
{
    /* Initialized data */

    static doublereal dmach[5] = { 2.23e-308,1.79e308,1.111e-16,2.222e-16,
	    .30102999566398119521 };

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D1MACH */
/* ***PURPOSE  Return floating point machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   D1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument, and can be referenced as follows: */

/*        D = D1MACH(I) */

/*   where I=1,...,5.  The (output) value of D above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude. */
/*   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   D1MACH( 3) = B**(-T), the smallest relative spacing. */
/*   D1MACH( 4) = B**(1-T), the largest relative spacing. */
/*   D1MACH( 5) = LOG10(B) */

/*   Assume double precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(14) = T, the number of base-B digits. */
/*   I1MACH(15) = EMIN, the smallest exponent E. */
/*   I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890213  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   900911  Added SUN 386i constants.  (WRB) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added CONVEX -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/*   010817  Elevated IEEE to highest importance; see next set of */
/*           comments below.  (DWL) */
/* ***END PROLOGUE  D1MACH */


/* Initial data here correspond to the IEEE standard.  The values for */
/* DMACH(1), DMACH(3) and DMACH(4) are slight upper bounds.  The value */
/* for DMACH(2) is a slight lower bound.  The value for DMACH(5) is */
/* a 20-digit approximation.  If one of the sets of initial data below */
/* is preferred, do the necessary commenting and uncommenting. (DWL) */

/* $$$      EQUIVALENCE (DMACH(1),SMALL(1)) */
/* $$$      EQUIVALENCE (DMACH(2),LARGE(1)) */
/* $$$      EQUIVALENCE (DMACH(3),RIGHT(1)) */
/* $$$      EQUIVALENCE (DMACH(4),DIVER(1)) */
/* $$$      EQUIVALENCE (DMACH(5),LOG10(1)) */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 / */
/*     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 / */
/*     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 / */
/*     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA SMALL(1) / ZC00800000 / */
/*     DATA SMALL(2) / Z000000000 / */
/*     DATA LARGE(1) / ZDFFFFFFFF / */
/*     DATA LARGE(2) / ZFFFFFFFFF / */
/*     DATA RIGHT(1) / ZCC5800000 / */
/*     DATA RIGHT(2) / Z000000000 / */
/*     DATA DIVER(1) / ZCC6800000 / */
/*     DATA DIVER(2) / Z000000000 / */
/*     DATA LOG10(1) / ZD00E730E7 / */
/*     DATA LOG10(2) / ZC77800DC0 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O0000000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O0007777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O7770000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O7777777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA SMALL(1) / Z"3001800000000000" / */
/*     DATA SMALL(2) / Z"3001000000000000" / */
/*     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" / */
/*     DATA LARGE(2) / Z"4FFE000000000000" / */
/*     DATA RIGHT(1) / Z"3FD2800000000000" / */
/*     DATA RIGHT(2) / Z"3FD2000000000000" / */
/*     DATA DIVER(1) / Z"3FD3800000000000" / */
/*     DATA DIVER(2) / Z"3FD3000000000000" / */
/*     DATA LOG10(1) / Z"3FFF9A209A84FBCF" / */
/*     DATA LOG10(2) / Z"3FFFF7988F8959AC" / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA SMALL(1) / 00564000000000000000B / */
/*     DATA SMALL(2) / 00000000000000000000B / */
/*     DATA LARGE(1) / 37757777777777777777B / */
/*     DATA LARGE(2) / 37157777777777777777B / */
/*     DATA RIGHT(1) / 15624000000000000000B / */
/*     DATA RIGHT(2) / 00000000000000000000B / */
/*     DATA DIVER(1) / 15634000000000000000B / */
/*     DATA DIVER(2) / 00000000000000000000B / */
/*     DATA LOG10(1) / 17164642023241175717B / */
/*     DATA LOG10(2) / 16367571421742254654B / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn OR -pd8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CC0000000000000' / */
/*     DATA DMACH(4) / Z'3CD0000000000000' / */
/*     DATA DMACH(5) / Z'3FF34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F900000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F910000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE CRAY */

/*     DATA SMALL(1) / 201354000000000000000B / */
/*     DATA SMALL(2) / 000000000000000000000B / */
/*     DATA LARGE(1) / 577767777777777777777B / */
/*     DATA LARGE(2) / 000007777777777777774B / */
/*     DATA RIGHT(1) / 376434000000000000000B / */
/*     DATA RIGHT(2) / 000000000000000000000B / */
/*     DATA DIVER(1) / 376444000000000000000B / */
/*     DATA DIVER(2) / 000000000000000000000B / */
/*     DATA LOG10(1) / 377774642023241175717B / */
/*     DATA LOG10(2) / 000007571421742254654B / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */
/*     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD - */
/*     STATIC DMACH(5) */

/*     DATA SMALL /    20K, 3*0 / */
/*     DATA LARGE / 77777K, 3*177777K / */
/*     DATA RIGHT / 31420K, 3*0 / */
/*     DATA DIVER / 32020K, 3*0 / */
/*     DATA LOG10 / 40423K, 42023K, 50237K, 74776K / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA DMACH(1) / '0000000000000010'X / */
/*     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X / */
/*     DATA DMACH(3) / '0000000000003CC0'X / */
/*     DATA DMACH(4) / '0000000000003CD0'X / */
/*     DATA DMACH(5) / '79FF509F44133FF3'X / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FORMAT */

/*     DATA DMACH(1) / '0010000000000000'X / */
/*     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X / */
/*     DATA DMACH(3) / '3CA0000000000000'X / */
/*     DATA DMACH(4) / '3CB0000000000000'X / */
/*     DATA DMACH(5) / '3FD34413509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/ */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/ */
/*     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/ */
/*     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/ */
/*     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/ */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /        128,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /       9344,           0 / */
/*     DATA DIVER(1), DIVER(2) /       9472,           0 / */
/*     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 / */

/*     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /         16,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /      15552,           0 / */
/*     DATA DIVER(1), DIVER(2) /      15568,           0 / */
/*     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 / */

/*     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */
/*     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION) */

/*     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X / */
/*     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X / */
/*     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X / */
/*     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X / */
/*     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA SMALL(1), SMALL(2) / '20000000, '00000201 / */
/*     DATA LARGE(1), LARGE(2) / '37777777, '37777577 / */
/*     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 / */
/*     DATA DIVER(1), DIVER(2) / '20000000, '00000334 / */
/*     DATA LOG10(1), LOG10(2) / '23210115, '10237777 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     THREE WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 / */
/*     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B / */
/*     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B / */
/*     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) /  40000B,       0 / */
/*     DATA SMALL(3), SMALL(4) /       0,       1 / */
/*     DATA LARGE(1), LARGE(2) /  77777B, 177777B / */
/*     DATA LARGE(3), LARGE(4) / 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) /  40000B,       0 / */
/*     DATA RIGHT(3), RIGHT(4) /       0,    225B / */
/*     DATA DIVER(1), DIVER(2) /  40000B,       0 / */
/*     DATA DIVER(3), DIVER(4) /       0,    227B / */
/*     DATA LOG10(1), LOG10(2) /  46420B,  46502B / */
/*     DATA LOG10(3), LOG10(4) /  76747B, 176377B / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B / */
/*     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B / */
/*     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B / */
/*     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B / */
/*     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF / */

/*     MACHINE CONSTANTS FOR THE IBM PC */
/*     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION */
/*     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087. */

/*     DATA SMALL(1) / 2.23D-308  / */
/*     DATA LARGE(1) / 1.79D+308  / */
/*     DATA RIGHT(1) / 1.11D-16   / */
/*     DATA DIVER(1) / 2.22D-16   / */
/*     DATA LOG10(1) / 0.301029995663981195D0 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    8388608,           0 / */
/*     DATA LARGE(1), LARGE(2) / 2147483647,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /  612368384,           0 / */
/*     DATA DIVER(1), DIVER(2) /  620756992,           0 / */
/*     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 / */

/*     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 / */
/*     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 / */
/*     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 / */
/*     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    128,      0 / */
/*     DATA SMALL(3), SMALL(4) /      0,      0 / */
/*     DATA LARGE(1), LARGE(2) /  32767,     -1 / */
/*     DATA LARGE(3), LARGE(4) /     -1,     -1 / */
/*     DATA RIGHT(1), RIGHT(2) /   9344,      0 / */
/*     DATA RIGHT(3), RIGHT(4) /      0,      0 / */
/*     DATA DIVER(1), DIVER(2) /   9472,      0 / */
/*     DATA DIVER(3), DIVER(4) /      0,      0 / */
/*     DATA LOG10(1), LOG10(2) /  16282,   8346 / */
/*     DATA LOG10(3), LOG10(4) / -31493, -12296 / */

/*     DATA SMALL(1), SMALL(2) / O000200, O000000 / */
/*     DATA SMALL(3), SMALL(4) / O000000, O000000 / */
/*     DATA LARGE(1), LARGE(2) / O077777, O177777 / */
/*     DATA LARGE(3), LARGE(4) / O177777, O177777 / */
/*     DATA RIGHT(1), RIGHT(2) / O022200, O000000 / */
/*     DATA RIGHT(3), RIGHT(4) / O000000, O000000 / */
/*     DATA DIVER(1), DIVER(2) / O022400, O000000 / */
/*     DATA DIVER(3), DIVER(4) / O000000, O000000 / */
/*     DATA LOG10(1), LOG10(2) / O037632, O020232 / */
/*     DATA LOG10(3), LOG10(4) / O102373, O147770 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE SUN 386i */

/*     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' / */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' / */
/*     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF' */
/*     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 / */

/* ***FIRST EXECUTABLE STATEMENT  D1MACH */
    if (*i__ < 1 || *i__ > 5) {
	xermsg_("SLATEC", "D1MACH", "I OUT OF BOUNDS", &c__1, &c__2, (ftnlen)
		6, (ftnlen)6, (ftnlen)15);
    }

    ret_val = dmach[*i__ - 1];
    return ret_val;

} /* d1mach_ */

/* DECK XGETUA */
/* Subroutine */ int xgetua_(integer *iunita, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, index;
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XGETUA */
/* ***PURPOSE  Return unit number(s) to which error messages are being */
/*            sent. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XGETUA-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XGETUA may be called to determine the unit number or numbers */
/*        to which error messages are being sent. */
/*        These unit numbers may have been set by a call to XSETUN, */
/*        or a call to XSETUA, or may be a default value. */

/*     Description of Parameters */
/*      --Output-- */
/*        IUNIT - an array of one to five unit numbers, depending */
/*                on the value of N.  A value of zero refers to the */
/*                default unit, as defined by the I1MACH machine */
/*                constant routine.  Only IUNIT(1),...,IUNIT(N) are */
/*                defined by XGETUA.  The values of IUNIT(N+1),..., */
/*                IUNIT(5) are not defined (for N .LT. 5) or altered */
/*                in any way by XGETUA. */
/*        N     - the number of units to which copies of the */
/*                error messages are being sent.  N will be in the */
/*                range from 1 to 5. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XGETUA */
/* ***FIRST EXECUTABLE STATEMENT  XGETUA */
    /* Parameter adjustments */
    --iunita;

    /* Function Body */
    *n = j4save_(&c__5, &c__0, &c_false);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	index = i__ + 4;
	if (i__ == 1) {
	    index = 3;
	}
	iunita[i__] = j4save_(&index, &c__0, &c_false);
/* L30: */
    }
    return 0;
} /* xgetua_ */

/* DECK FDUMP */
/* Subroutine */ int fdump_()
{
/* ***BEGIN PROLOGUE  FDUMP */
/* ***PURPOSE  Symbolic dump (should be locally written). */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3 */
/* ***TYPE      ALL (FDUMP-A) */
/* ***KEYWORDS  ERROR, XERMSG */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*        ***Note*** Machine Dependent Routine */
/*        FDUMP is intended to be replaced by a locally written */
/*        version which produces a symbolic dump.  Failing this, */
/*        it should be replaced by a version which prints the */
/*        subprogram nesting list.  Note that this dump must be */
/*        printed on each of up to five files, as indicated by the */
/*        XGETUA routine.  See XSETUA and XGETUA for details. */

/*     Written by Ron Jones, with SLATEC Common Math Library Subcommittee */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  FDUMP */
/* ***FIRST EXECUTABLE STATEMENT  FDUMP */
    return 0;
} /* fdump_ */

/* DECK J4SAVE */
integer j4save_(integer *iwhich, integer *ivalue, logical *iset)
{
    /* Initialized data */

    static integer iparam[9] = { 0,2,0,10,1,0,0,0,0 };

    /* System generated locals */
    integer ret_val;

/* ***BEGIN PROLOGUE  J4SAVE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Save or recall global variables needed by error */
/*            handling routines. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***TYPE      INTEGER (J4SAVE-I) */
/* ***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        J4SAVE saves and recalls several global variables needed */
/*        by the library error handling routines. */

/*     Description of Parameters */
/*      --Input-- */
/*        IWHICH - Index of item desired. */
/*                = 1 Refers to current error number. */
/*                = 2 Refers to current error control flag. */
/*                = 3 Refers to current unit number to which error */
/*                    messages are to be sent.  (0 means use standard.) */
/*                = 4 Refers to the maximum number of times any */
/*                     message is to be printed (as set by XERMAX). */
/*                = 5 Refers to the total number of units to which */
/*                     each error message is to be written. */
/*                = 6 Refers to the 2nd unit for error messages */
/*                = 7 Refers to the 3rd unit for error messages */
/*                = 8 Refers to the 4th unit for error messages */
/*                = 9 Refers to the 5th unit for error messages */
/*        IVALUE - The value to be set for the IWHICH-th parameter, */
/*                 if ISET is .TRUE. . */
/*        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE */
/*                 given the value, IVALUE.  If ISET=.FALSE., the */
/*                 IWHICH-th parameter will be unchanged, and IVALUE */
/*                 is a dummy parameter. */
/*      --Output-- */
/*        The (old) value of the IWHICH-th parameter will be returned */
/*        in the function value, J4SAVE. */

/* ***SEE ALSO  XERMSG */
/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900205  Minor modifications to prologue.  (WRB) */
/*   900402  Added TYPE section.  (WRB) */
/*   910411  Added KEYWORDS section.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  J4SAVE */
/* ***FIRST EXECUTABLE STATEMENT  J4SAVE */
    ret_val = iparam[(0 + (0 + (*iwhich - 1 << 2))) / 4];
    if (*iset) {
	iparam[*iwhich - 1] = *ivalue;
    }
    return ret_val;
} /* j4save_ */

/* DECK XERCNT */
/* Subroutine */ int xercnt_(char *librar, char *subrou, char *messg, integer 
	*nerr, integer *level, integer *kontrl, ftnlen librar_len, ftnlen 
	subrou_len, ftnlen messg_len)
{
/* ***BEGIN PROLOGUE  XERCNT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Allow user control over handling of errors. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERCNT-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        Allows user control over handling of individual errors. */
/*        Just after each message is recorded, but before it is */
/*        processed any further (i.e., before it is printed or */
/*        a decision to abort is made), a call is made to XERCNT. */
/*        If the user has provided his own version of XERCNT, he */
/*        can then override the value of KONTROL used in processing */
/*        this message by redefining its value. */
/*        KONTRL may be set to any value from -2 to 2. */
/*        The meanings for KONTRL are the same as in XSETF, except */
/*        that the value of KONTRL changes only for this message. */
/*        If KONTRL is set to a value outside the range from -2 to 2, */
/*        it will be moved back into that range. */

/*     Description of Parameters */

/*      --Input-- */
/*        LIBRAR - the library that the routine is in. */
/*        SUBROU - the subroutine that XERMSG is being called from */
/*        MESSG  - the first 20 characters of the error message. */
/*        NERR   - same as in the call to XERMSG. */
/*        LEVEL  - same as in the call to XERMSG. */
/*        KONTRL - the current value of the control flag as set */
/*                 by a call to XSETF. */

/*      --Output-- */
/*        KONTRL - the new value of KONTRL.  If KONTRL is not */
/*                 defined, it will remain at its original value. */
/*                 This changed value of control affects only */
/*                 the current occurrence of the current message. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900206  Routine changed from user-callable to subsidiary.  (WRB) */
/*   900510  Changed calling sequence to include LIBRARY and SUBROUTINE */
/*           names, changed routine name from XERCTL to XERCNT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERCNT */
/* ***FIRST EXECUTABLE STATEMENT  XERCNT */
    return 0;
} /* xercnt_ */

/* DECK XERHLT */
/* Subroutine */ int xerhlt_(char *messg, ftnlen messg_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

/* ***BEGIN PROLOGUE  XERHLT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Abort program execution and print error message. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERHLT-A) */
/* ***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        ***Note*** machine dependent routine */
/*        XERHLT aborts the execution of the program. */
/*        The error message causing the abort is given in the calling */
/*        sequence, in case one needs it for printing on a dayfile, */
/*        for example. */

/*     Description of Parameters */
/*        MESSG is as in XERMSG. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900206  Routine changed from user-callable to subsidiary.  (WRB) */
/*   900510  Changed calling sequence to delete length of character */
/*           and changed routine name from XERABT to XERHLT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERHLT */
/* ***FIRST EXECUTABLE STATEMENT  XERHLT */
    s_stop("", (ftnlen)0);
    return 0;
} /* xerhlt_ */

#ifdef __cplusplus
	}
#endif
