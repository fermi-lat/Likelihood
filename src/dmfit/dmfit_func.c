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
    real phidif[72576]	/* was [252][24][12] */;
    integer hasmooth;
} hasim_;

#define hasim_1 hasim_

/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__23 = 23;

/* -----DMFIT dN/dE routine */
doublereal dmfit_de__(doublereal *mx, integer *ch, doublereal *ee)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double d_lg10(doublereal *);

    /* Local variables */
    extern doublereal yieldget_(integer *, integer *, integer *);
    static integer i__;
    static doublereal z__, mi[24];
    static integer zi, zn, m1i, m2i;
    static doublereal mp1, mp2;
    static integer hsm;
    static doublereal tmp, zpl, phi1, phi2, flux;
    static integer chref;
    static doublereal eeold;
    extern /* Subroutine */ int ifind_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *);
    static doublereal mxold, zindex[504]	/* was [252][2] */;

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
/*      real*8 mt,mw,mz */
/* -----lowest mass index */
/*      real*8 milow(10) */
/* -----this is the reference channel */
/* -----variables to initialize tables */
/* ,ntype */
/*      integer j,k,l */
/* -----number of decades tabulated */
/* -----backup masses and energies for low energy extrapolation */
/* -----function that computes the differential \gamma flux from e+e- */
/*      real*8 llg */
/* -----data tables */
/* ******************************************************************* */
/*     backup the input energy and mass values */
    eeold = *ee;
    mxold = *mx;
/* -----mass cut: lower limit, all channels but not e+e- */
/* $$$      if(MX.lt.10.d0.and.CH.ne.1) then */
/* $$$c         write(*,*) 'WARNING: MASS IS TOO LOW!' */
/* $$$c         write(*,*) 'USING EXTRAPOLATION OF MC DATA!' */
/* $$$c-----here we do the barbaric extrapolation to lower masses */
/* $$$         EE=EE/MX*10.d0 */
/* $$$         MX=10.d0 */
/* $$$c         stop */
/* $$$      endif */
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
/* $$$      mt=174.3d0 */
/* $$$      mw=80.33d0 */
/* $$$      mz=91.187d0 */
/* -----switches to CH 4 if mass limit inconsistent */
/* $$$      if(CH.eq.5) then */
/* $$$         if(MX.lt.mt) then */
/* $$$            CH=4 */
/* $$$            write(*,*) 'ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!' */
/* $$$            write(*,*) '       SWITHCING TO DEFAULT, B\bar B!' */
/* $$$         endif */
/* $$$      endif */
/* $$$ */
/* $$$      if(CH.eq.7) then */
/* $$$         if(MX.lt.mw) then */
/* $$$            CH=4 */
/* $$$            write(*,*) 'ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!' */
/* $$$            write(*,*) '       SWITHCING TO DEFAULT, B\bar B!' */
/* $$$         endif */
/* $$$      endif */
/* $$$ */
/* $$$      if(CH.eq.8) then */
/* $$$         if(MX.lt.mz) then */
/* $$$            CH=4 */
/* $$$            write(*,*) 'ERROR: CHANNEL NOT KINEMATICALLY ALLOWED!' */
/* $$$            write(*,*) '       SWITHCING TO DEFAULT, B\bar B!' */
/* $$$         endif */
/* $$$      endif */
/* -----translate the channels */
/* -----this is the default channel, b \bar b */
/* $$$      chref=4 */
/* $$$      if(CH.eq.1) then */
/* $$$c-----for the e+e- channel, go ahead and compute it! */
/* $$$         dmfit_de=llg(EE/MX,MX,0.511d-3) */
/* $$$         return */
/* $$$ */
/* $$$c-----for muon channel at low mass, the MonteCarlo statistics */
/* $$$c-----resulting in the DMFIT table is too low, so we switch to */
/* $$$c-----the same equation as the electron channel instead. */
/* $$$      elseif(CH.eq.2.and.MX.lt.10.d0) then */
/* $$$         dmfit_de=llg(EE/MX,MX,0.1057d0) */
/* $$$         return */
/* $$$      elseif(CH.eq.2) then */
/* $$$         chref=7 */
/* $$$      elseif(CH.eq.3) then */
/* $$$         chref=4 */
/* $$$      elseif(CH.eq.4) then */
/* $$$         chref=2 */
/* $$$      elseif(CH.eq.5) then */
/* $$$         chref=3 */
/* $$$      elseif(CH.eq.6) then */
/* $$$         chref=8 */
/* $$$      elseif(CH.eq.7) then */
/* $$$         chref=5 */
/* $$$      elseif(CH.eq.8) then */
/* $$$         chref=6 */
/* $$$      elseif(CH.eq.9) then */
/* $$$         chref=1 */
/* $$$      elseif(CH.eq.10) then */
/* $$$         chref=9 */
/* $$$      elseif(CH.eq.11) then */
/* $$$         chref=10 */
/* $$$      else */
/* $$$         chref=4 */
/* $$$      endif */
    if (*ch == 1) {
/* -----e+e- */
	chref = 9;
    } else if (*ch == 2) {
/* -----mu+mu- */
	chref = 7;
    } else if (*ch == 3) {
/* -----tau+tau- */
	chref = 4;
    } else if (*ch == 4) {
/* -----b bbar */
	chref = 2;
    } else if (*ch == 5) {
/* ------t tbar */
	chref = 3;
    } else if (*ch == 6) {
/* -----g g */
	chref = 8;
    } else if (*ch == 7) {
/* -----W+W- */
	chref = 5;
    } else if (*ch == 8) {
/* -----ZZ */
	chref = 6;
    } else if (*ch == 9) {
/* -----c cbar */
	chref = 1;
    } else if (*ch == 10) {
/* -----s sbar */
	chref = 10;
    } else if (*ch == 11) {
/* -----u ubar */
	chref = 11;
    } else if (*ch == 12) {
/* -----d dbar */
	chref = 12;
    } else {
	chref = 4;
    }
/* ******************************************************************** */
/* -----this initializes and loads the MC simulation data for gammas */
/* -----masses for simulation corresponding to mass index i */
/* $$$        mi(1)=10.0 */
/* $$$        mi(2)=25.0 */
/* $$$        mi(3)=50.0 */
/* $$$        mi(4)=80.3 */
/* $$$        mi(5)=91.2 */
/* $$$        mi(6)=100.0 */
/* $$$        mi(7)=150.0 */
/* $$$        mi(8)=176.0 */
/* $$$        mi(9)=200.0 */
/* $$$        mi(10)=250.0 */
/* $$$        mi(11)=350.0 */
/* $$$        mi(12)=500.0 */
/* $$$        mi(13)=750.0 */
/* $$$        mi(14)=1000.0 */
/* $$$        mi(15)=1500.0 */
/* $$$        mi(16)=2000.0 */
/* $$$        mi(17)=3000.0 */
/* $$$        mi(18)=5000.0 */
    mi[0] = 2.;
    mi[1] = 4.;
    mi[2] = 6.;
    mi[3] = 8.;
    mi[4] = 10.;
    mi[5] = 25.;
    mi[6] = 50.;
    mi[7] = 80.3;
    mi[8] = 91.2;
    mi[9] = 100.;
    mi[10] = 150.;
    mi[11] = 176.;
    mi[12] = 200.;
    mi[13] = 250.;
    mi[14] = 350.;
    mi[15] = 500.;
    mi[16] = 750.;
    mi[17] = 1e3;
    mi[18] = 1500.;
    mi[19] = 2e3;
    mi[20] = 3e3;
    mi[21] = 5e3;
    mi[22] = 7e3;
    mi[23] = 1e4;
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
    if (*ee >= *mx || *mx > mi[23]) {
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
    ifind_(mx, mi, &tmp, &m1i, &c__1, &c__23);
    mp1 = mi[m1i - 1];
    m2i = m1i + 1;
    mp2 = mi[m2i - 1];
    if (*mx == mi[23]) {
	m1i = 24;
	m2i = 24;
	mp1 = mi[23];
	mp2 = mp1;
	i__1 = zi + 1;
	flux = (1.f - zpl) * yieldget_(&zi, &m1i, &chref) + zpl * yieldget_(&
		i__1, &m1i, &chref);
/* $$$          FLUX=FLUX*lge/(ndec*EE) */
	flux /= *ee;
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
/* $$$          FLUX=FLUX*lge/(ndec*EE) */
	flux /= *ee;
    }
/* L105: */
/* $$$        if(MXOLD.lt.10.d0.and.CH.ne.1) then */
/* $$$           FLUX=FLUX/(MXOLD/10.d0) */
/* $$$        endif */
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

