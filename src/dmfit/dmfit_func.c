/* dmfit_func.f -- translated by f2c (version 20090411).
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
    real phidif[45360]	/* was [252][18][10] */;
    integer hasmooth;
} hasim_;

#define hasim_1 hasim_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static doublereal c_b20 = 5.11e-4;
static integer c_n1 = -1;
static integer c__17 = 17;

/* -----DMFIT dN/dE routine */
doublereal dmfit_de__(doublereal *mx, integer *ch, doublereal *ee)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double d_lg10(doublereal *);

    /* Local variables */
    extern doublereal yieldget_(integer *, integer *, integer *);
    static integer i__;
    static doublereal z__, mi[18], mt;
    static integer zi;
    static doublereal mw, mz;
    static integer zn, m1i, m2i;
    static doublereal mp1, mp2;
    extern doublereal llg_(doublereal *, doublereal *, doublereal *);
    static integer hsm;
    static doublereal tmp, zpl, phi1, phi2, flux;
    static integer chref;
    static doublereal eeold;
    extern /* Subroutine */ int ifind_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *);
    static doublereal mxold, zindex[504]	/* was [252][2] */;

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };


/* -----ARGUMENTS: */
/* -----WIMP mass */
/*     since the MC simulations run between 10 and 5000 GeV, we impose */
/*     the corresponding cuts */
/* -----DM annihilation mode */
/* -----the following channels, corresponding to WIMP annihilation modes */
/*     are implemented: */

/*     1: e^+e^- */
/*     2: mu^+ mu^- */
/*     3: tau^+ tau^- */
/*     4: b \bar b */
/*     5: t \bar t */
/*     6: gluon gluon */
/*     7: W+ W- */
/*     8: Z Z */
/*     9: c \bar c */
/*     10: cosmo b \bar b */
/*     11: cosmo gam gam Line */

/*     Notice that if a channel is not kinematically available */
/*     (e.g., ch 7 with MX < M_W), then automatically we fit */
/*     with default channel 4. This applies to CH=5,7,8 */
/* -----Smoothing procedure (hard-wired in this routine!) */
/*     HSM = 0: simply takes the closest point in the table */
/*     HSM = 1: interpolates between 3 points */
/*     HSM = 2: interpolates between 5 points */
/* -----User-supplied photon energy, in GeV */
/* -----Photon flux (differential) */
/* -----whether or not to smooth the interpolation */
/* -----masses of top, W, Z */
/* -----lowest mass index */
/*      real*8 milow(10) */
/* -----this is the reference channel */
/* -----variables to initialize tables */
/* ,ntype */
/*      integer j,k,l */
/* -----number of decades tabulated */
/* -----backup masses and energies for low energy extrapolation */
/* -----function that computes the differential \gamma flux from e+e- */
/* -----data tables */
/* ******************************************************************* */
/*     backup the input energy and mass values */
    eeold = *ee;
    mxold = *mx;
/* -----mass cut: lower limit, all channels but not e+e- */
    if (*mx < 10. && *ch != 1) {
/*         write(*,*) 'WARNING: MASS IS TOO LOW!' */
/*         write(*,*) 'USING EXTRAPOLATION OF MC DATA!' */
/* -----here we do the barbaric extrapolation to lower masses */
	*ee = *ee / *mx * 10.;
	*mx = 10.;
/*         stop */
    }
    hsm = 0;
/* -----imposes the smoothing to use in the interpolation */
    if (hsm == 0) {
	hasim_1.hasmooth = 0;
    } else if (hsm == 1) {
	hasim_1.hasmooth = 1;
    } else if (hsm == 2) {
	hasim_1.hasmooth = 2;
    } else {
	hasim_1.hasmooth = 0;
    }
/* -----set the masses for top, W, Z */
    mt = 174.3;
    mw = 80.33;
    mz = 91.187;
/* -----switches to CH 4 if mass limit inconsistent */
    if (*ch == 5) {
	if (*mx < mt) {
	    *ch = 4;
	    s_wsle(&io___7);
	    do_lio(&c__9, &c__1, "ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!", 
		    (ftnlen)41);
	    e_wsle();
	    s_wsle(&io___8);
	    do_lio(&c__9, &c__1, "       SWITHCING TO DEFAULT, B\bar B!", (
		    ftnlen)36);
	    e_wsle();
	}
    }
    if (*ch == 7) {
	if (*mx < mw) {
	    *ch = 4;
	    s_wsle(&io___9);
	    do_lio(&c__9, &c__1, "ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!", 
		    (ftnlen)41);
	    e_wsle();
	    s_wsle(&io___10);
	    do_lio(&c__9, &c__1, "       SWITHCING TO DEFAULT, B\bar B!", (
		    ftnlen)36);
	    e_wsle();
	}
    }
    if (*ch == 8) {
	if (*mx < mz) {
	    *ch = 4;
	    s_wsle(&io___11);
	    do_lio(&c__9, &c__1, "ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!", 
		    (ftnlen)41);
	    e_wsle();
	    s_wsle(&io___12);
	    do_lio(&c__9, &c__1, "       SWITHCING TO DEFAULT, B\bar B!", (
		    ftnlen)36);
	    e_wsle();
	}
    }
/* -----translate the channels */
/* -----this is the default channel, b \bar b */
    chref = 4;
    if (*ch == 1) {
/* -----for the e+e- channel, go ahead and compute it! */
	d__1 = *ee / *mx;
	ret_val = llg_(&d__1, mx, &c_b20);
	return ret_val;
    } else if (*ch == 2) {
	chref = 7;
    } else if (*ch == 3) {
	chref = 4;
    } else if (*ch == 4) {
	chref = 2;
    } else if (*ch == 5) {
	chref = 3;
    } else if (*ch == 6) {
	chref = 8;
    } else if (*ch == 7) {
	chref = 5;
    } else if (*ch == 8) {
	chref = 6;
    } else if (*ch == 9) {
	chref = 1;
    } else if (*ch == 10) {
	chref = 9;
    } else if (*ch == 11) {
	chref = 10;
    } else {
	chref = 4;
    }
/* ******************************************************************** */
/* -----this initializes and loads the MC simulation data for gammas */
/* -----masses for simulation corresponding to mass index i */
    mi[0] = 10.f;
    mi[1] = 25.f;
    mi[2] = 50.f;
    mi[3] = 80.3f;
    mi[4] = 91.2f;
    mi[5] = 100.f;
    mi[6] = 150.f;
    mi[7] = 176.f;
    mi[8] = 200.f;
    mi[9] = 250.f;
    mi[10] = 350.f;
    mi[11] = 500.f;
    mi[12] = 750.f;
    mi[13] = 1e3f;
    mi[14] = 1500.f;
    mi[15] = 2e3f;
    mi[16] = 3e3f;
    mi[17] = 5e3f;
/* -----initialize eindex array where the energies for the bins are stored */
/* -----integrated yields (lower end of each bin) */
    zn = 250;
    i__1 = zn;
    for (i__ = -1; i__ <= i__1; ++i__) {
	zindex[i__ + 1] = (doublereal) i__ / (doublereal) zn;
    }
    i__1 = zn;
    for (i__ = -1; i__ <= i__1; ++i__) {
	zindex[i__ + 253] = (doublereal) i__ / (doublereal) zn + .5 / (
		doublereal) zn;
    }
/* ******************************************************************** */
/* -----interpolation: lower energy */
    if (*ee >= *mx) {
	ret_val = 0.;
	return ret_val;
    }
    d__1 = *ee / *mx;
    z__ = (d_lg10(&d__1) + 10.) / 10.;
    i__1 = zn - 1;
    ifind_(&z__, &zindex[252], &zpl, &zi, &c_n1, &i__1);
    if (zi == -5 || zi >= zn) {
	ret_val = 0.;
	return ret_val;
    }
    ifind_(mx, mi, &tmp, &m1i, &c__1, &c__17);
    mp1 = mi[m1i - 1];
    m2i = m1i + 1;
    mp2 = mi[m2i - 1];
    if (*mx >= mi[17]) {
	m1i = 18;
	m2i = 18;
	mp1 = mi[17];
	mp2 = mp1;
	i__1 = zi + 1;
	flux = (1.f - zpl) * yieldget_(&zi, &m1i, &chref) + zpl * yieldget_(&
		i__1, &m1i, &chref);
	flux = flux * .434294481903 / (*ee * 10.);
    } else {
/*           write(*,*) CH,zi,m1i,EE,yieldget(zi,m1i,chref) */
	i__1 = zi + 1;
	phi1 = (1.f - zpl) * yieldget_(&zi, &m1i, &chref) + zpl * yieldget_(&
		i__1, &m1i, &chref);
	i__1 = zi + 1;
	phi2 = (1.f - zpl) * yieldget_(&zi, &m2i, &chref) + zpl * yieldget_(&
		i__1, &m2i, &chref);
	flux = phi1 + (phi2 - phi1) * (*mx - mp1) * (*mx + mp1) / ((mp2 - mp1)
		 * (mp2 + mp1));
	flux = flux * .434294481903 / (*ee * 10.);
    }
/* L105: */
    if (mxold < 10. && *ch != 1) {
	flux /= mxold / 10.;
    }
    ret_val = flux;
/* -----here we restore the input values if they were altered */
    *ee = eeold;
    *mx = mxold;
    return ret_val;
} /* dmfit_de__ */

/* -----DMFIT d(dN/dE)/dM routine */
doublereal dmfit_dm__(doublereal *mx, integer *ch, doublereal *ee)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern doublereal dmfit_de__(doublereal *, integer *, doublereal *);
    static doublereal d1, d2, m1, m2, inc1, inc2, delta;

/* -----ARGUMENTS: */
/* -----the other function... */
/* -----auxiliary variables */
/* ******************************************************************* */
/*     this sets the "delta" for the computation of the derivative */
    delta = 1e-4;
    m1 = *mx / (delta + 1.);
    d1 = *mx - m1;
    m2 = *mx * (delta + 1.);
    d2 = m2 - *mx;
    if (m1 > *ee) {
	inc1 = dmfit_de__(mx, ch, ee) - dmfit_de__(&m1, ch, ee);
	inc1 /= d1;
	inc2 = dmfit_de__(&m2, ch, ee) - dmfit_de__(mx, ch, ee);
	inc2 /= d2;
	ret_val = (inc1 + inc2) / 2.;
    } else {
	inc2 = dmfit_de__(&m2, ch, ee) - dmfit_de__(mx, ch, ee);
	inc2 /= d2;
	ret_val = inc2;
    }
    return ret_val;
} /* dmfit_dm__ */

/* ******************************************************************** */
/* -----DMFIT \int_E^MX (dN/dE) dE */
doublereal dmfit_deint__(doublereal *mx, integer *ch, doublereal *ee)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    extern doublereal dmfit_de__(doublereal *, integer *, doublereal *);
    static doublereal integral;
    static integer ii;
    static doublereal ene1, ene2, dnde1, dnde2, delta;
    static integer nstep;

/* -----ARGUMENTS: */
/* -----the other function, which we want to integrate... */
/* -----auxiliary variables */
/*     sets to 0 the integral */
    integral = 0.;
/*     defines the number of steps - say 100 */
    nstep = 100;
/*     case where energy is above the mass - flux is 0 */
    if (*ee >= *mx) {
	ret_val = 0.;
	return ret_val;
    }
/*     integration */
    i__1 = nstep;
    for (ii = 0; ii <= i__1; ++ii) {
	ene1 = *ee * exp((log(*mx) - log(*ee)) * (doublereal) ii / (
		doublereal) nstep);
	ene2 = *ee * exp((log(*mx) - log(*ee)) * (doublereal) (ii + 1) / (
		doublereal) nstep);
	delta = ene2 - ene1;
	dnde1 = dmfit_de__(mx, ch, &ene1);
	dnde2 = dmfit_de__(mx, ch, &ene2);
	integral += (dnde1 + dnde2) * delta / 2.;
    }
    ret_val = integral;
    return ret_val;
} /* dmfit_deint__ */

/* ******************************************************************** */
/* -----DMFIT \int_E^MX (dN/dM) dE */
doublereal dmfit_dmint__(doublereal *mx, integer *ch, doublereal *ee)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    extern doublereal dmfit_dm__(doublereal *, integer *, doublereal *);
    static doublereal integral;
    static integer ii;
    static doublereal ene1, ene2, dnde1, dnde2, delta;
    static integer nstep;

/* -----ARGUMENTS: */
/* -----the other function, which we want to integrate... */
/* -----auxiliary variables */
/*     sets to 0 the integral */
    integral = 0.;
/*     defines the number of steps - say 100 */
    nstep = 100;
/*     case where energy is above the mass - flux is 0 */
    if (*ee >= *mx) {
	ret_val = 0.;
	return ret_val;
    }
/*     integration */
    i__1 = nstep;
    for (ii = 0; ii <= i__1; ++ii) {
	ene1 = *ee * exp((log(*mx) - log(*ee)) * (doublereal) ii / (
		doublereal) nstep);
	ene2 = *ee * exp((log(*mx) - log(*ee)) * (doublereal) (ii + 1) / (
		doublereal) nstep);
	delta = ene2 - ene1;
	dnde1 = dmfit_dm__(mx, ch, &ene1);
	dnde2 = dmfit_dm__(mx, ch, &ene2);
	integral += (dnde1 + dnde2) * delta / 2.;
    }
    ret_val = integral;
    return ret_val;
} /* dmfit_dmint__ */

