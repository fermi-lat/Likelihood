/* foo.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c/f2c.h"

/* Common Block Declarations */

struct {
    char cpnam[1000];
} mn7nam_;

#define mn7nam_1 mn7nam_

struct {
    doublereal u[100], alim[100], blim[100];
} mn7ext_;

#define mn7ext_1 mn7ext_

struct {
    doublereal erp[100], ern[100], werr[100], globcc[100];
} mn7err_;

#define mn7err_1 mn7err_

struct {
    integer nvarl[100], niofex[100], nexofi[100];
} mn7inx_;

#define mn7inx_1 mn7inx_

struct {
    doublereal x[100], xt[100], dirin[100];
} mn7int_;

#define mn7int_1 mn7int_

struct {
    doublereal xs[100], xts[100], dirins[100];
} mn7fx2_;

#define mn7fx2_1 mn7fx2_

struct {
    doublereal grd[100], g2[100], gstep[100], gin[100], dgrd[100];
} mn7der_;

#define mn7der_1 mn7der_

struct {
    doublereal grds[100], g2s[100], gsteps[100];
} mn7fx3_;

#define mn7fx3_1 mn7fx3_

struct {
    integer ipfix[100], npfix;
} mn7fx1_;

#define mn7fx1_1 mn7fx1_

struct {
    doublereal vhmat[5050];
} mn7var_;

#define mn7var_1 mn7var_

struct {
    doublereal vthmat[5050];
} mn7vat_;

#define mn7vat_1 mn7vat_

struct {
    doublereal p[10100]	/* was [100][101] */, pstar[100], pstst[100], pbar[
	    100], prho[100];
} mn7sim_;

#define mn7sim_1 mn7sim_

struct {
    integer maxint, npar, maxext, nu;
} mn7npr_;

#define mn7npr_1 mn7npr_

struct {
    integer isysrd, isyswr, isyssa, npagwd, npagln, newpag;
} mn7iou_;

#define mn7iou_1 mn7iou_

struct {
    integer istkrd[10], nstkrd, istkwr[10], nstkwr;
} mn7io2_;

#define mn7io2_1 mn7io2_

struct {
    char cfrom[8], cstatu[10], ctitl[50], cword[20], cundef[10], cvrsn[6], 
	    covmes[88];
} mn7tit_;

#define mn7tit_1 mn7tit_

struct {
    integer isw[7], idbg[11], nblock, icomnd;
} mn7flg_;

#define mn7flg_1 mn7flg_

struct {
    doublereal amin, up, edm, fval3, epsi, apsi, dcovar;
} mn7min_;

#define mn7min_1 mn7min_

struct {
    integer nfcn, nfcnmx, nfcnlc, nfcnfr, itaur, istrat, nwrmes[2];
} mn7cnv_;

#define mn7cnv_1 mn7cnv_

struct {
    doublereal word7[30];
} mn7arg_;

#define mn7arg_1 mn7arg_

struct {
    logical lwarn, lrepor, limset, lnolim, lnewmn, lphead;
} mn7log_;

#define mn7log_1 mn7log_

struct {
    doublereal epsmac, epsma2, vlimlo, vlimhi, undefi, bigedm, updflt;
} mn7cns_;

#define mn7cns_1 mn7cns_

struct {
    doublereal xpt[101], ypt[101];
} mn7rpt_;

#define mn7rpt_1 mn7rpt_

struct {
    char chpt[101];
} mn7cpt_;

#define mn7cpt_1 mn7cpt_

struct {
    doublereal xmidcr, ymidcr, xdircr, ydircr;
    integer ke1cr, ke2cr;
} mn7xcr_;

#define mn7xcr_1 mn7xcr_

struct {
    char origin[200]	/* was [10][2] */, warmes[1200]	/* was [10][2] */;
} mn7wrc_;

#define mn7wrc_1 mn7wrc_

struct {
    integer nfcwar[20]	/* was [10][2] */, icirc[2];
} mn7wri_;

#define mn7wri_1 mn7wri_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__3 = 3;
static real c_b66 = (float)10.;
static integer c__20 = 20;
static integer c__30 = 30;
static integer c__0 = 0;
static integer c__5 = 5;
static integer c__100 = 100;
static doublereal c_b1209 = .05;
static integer c__10 = 10;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:18  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni */
/* Minuit */


/* Subroutine */ int minuit_0_(int n__, S_fp fcn, U_fp futil, integer *i1, 
	integer *i2, integer *i3)
{
    /* Initialized data */

    static char cwhyxt[40+1] = "FOR UNKNOWN REASONS                     ";
    static integer jsysrd = 5;
    static integer jsyswr = 6;
    static integer jsyssa = 7;

    /* Format strings */
    static char fmt_280[] = "(/\002 MINUIT WARNING: PROBABLE ERROR IN USER F\
UNCTION.\002/\002 FOR FIXED VALUES OF PARAMETERS, FCN IS TIME-DEPENDENT\002\
/\002 F =\002,e22.14,\002 FOR FIRST CALL\002/\002 F =\002,e22.14,\002 FOR SE\
COND CALL.\002/)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2];

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer i_indx(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static doublereal fnew;
    static integer nparx;
    static doublereal fzero, first;
    extern /* Subroutine */ int mnread_(S_fp, integer *, integer *, U_fp), 
	    mncler_(), mninit_(integer *, integer *, integer *);
    static integer iflgut;
    extern /* Subroutine */ int mninex_(doublereal *), mnprin_(integer *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 0, 0, "(1X,75(1H*))", 0 };
    static cilist io___6 = { 0, 0, 0, "(1X,75(1H*))", 0 };
    static cilist io___7 = { 0, 0, 0, "(26X,A,I4)", 0 };
    static cilist io___8 = { 0, 0, 0, "(1X,75(1H*))", 0 };
    static cilist io___10 = { 0, 0, 0, "(/A,A)", 0 };
    static cilist io___14 = { 0, 0, 0, "(/A,A/)", 0 };
    static cilist io___16 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___17 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___18 = { 0, 0, 0, "(A,A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



/*  CPNAM   Parameter name (10 characters) */
/*  U       External (visible to user in FCN) value of parameter */
/*  ALIM, BLIM Lower and upper parameter limits. If both zero, no limits. */
/*  ERP,ERN Positive and negative MINOS errors, if calculated. */
/*  WERR    External parameter error (standard deviation, defined by UP) */
/*  GLOBCC  Global Correlation Coefficient */
/*  NVARL   =-1 if parameter undefined,      =0 if constant, */
/*          = 1 if variable without limits,  =4 if variable with limits */
/*   (Note that if parameter has been fixed, NVARL=1 or =4, and NIOFEX=0) */
/*  NIOFEX  Internal parameter number, or zero if not currently variable */
/*  NEXOFI  External parameter number for currently variable parameters */
/*  X, XT   Internal parameter values (X are sometimes saved in XT) */
/*  DIRIN   (Internal) step sizes for current step */
/*  variables with names ending in ..S are saved values for fixed params */
/*  VHMAT   (Internal) error matrix stored as Half MATrix, since */
/*                it is symmetric */
/*  VTHMAT  VHMAT is sometimes saved in VTHMAT, especially in MNMNOT */

/*  ISW definitions: */
/*      ISW(1) =0 normally, =1 means CALL LIMIT EXCEEDED */
/*      ISW(2) =0 means no error matrix */
/*             =1 means only approximate error matrix */
/*             =2 means full error matrix, but forced pos-def. */
/*             =3 means good normal full error matrix exists */
/*      ISW(3) =0 if Minuit is calculating the first derivatives */
/*             =1 if first derivatives calculated inside FCN */
/*      ISW(4) =-1 if most recent minimization did not converge. */
/*             = 0 if problem redefined since most recent minimization. */
/*             =+1 if most recent minimization did converge. */
/*      ISW(5) is the PRInt level.  See SHO PRIntlevel */
/*      ISW(6) = 0 for batch mode, =1 for interactive mode */
/*                      =-1 for originally interactive temporarily batch */

/*  LWARN is true if warning messges are to be put out (default=true) */
/*            SET WARN turns it on, set NOWarn turns it off */
/*  LREPOR is true if exceptional conditions are put out (default=false) */
/*            SET DEBUG turns it on, SET NODebug turns it off */
/*  LIMSET is true if a parameter is up against limits (for MINOS) */
/*  LNOLIM is true if there are no limits on any parameters (not yet used) */
/*  LNEWMN is true if the previous process has unexpectedly improved FCN */
/*  LPHEAD is true if a heading should be put out for the next parameter */
/*        definition, false if a parameter has just been defined */

    switch(n__) {
	case 1: goto L_mintio;
	}

/*                                 . . . . . . . . . . initialize minuit */
    io___5.ciunit = jsyswr;
    s_wsfe(&io___5);
    e_wsfe();
    mninit_(&jsysrd, &jsyswr, &jsyssa);
/*                                      . . . . initialize new data block */
L100:
    io___6.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___6);
    e_wsfe();
    ++mn7flg_1.nblock;
    io___7.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___7);
    do_fio(&c__1, "MINUIT DATA BLOCK NO.", (ftnlen)21);
    do_fio(&c__1, (char *)&mn7flg_1.nblock, (ftnlen)sizeof(integer));
    e_wsfe();
    io___8.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___8);
    e_wsfe();
/*               . . . . . . . . . . .   set parameter lists to undefined */
    mncler_();
/*                                             . . . . . . . . read title */
    mnread_((S_fp)fcn, &c__1, &iflgut, (U_fp)futil);
    if (iflgut == 2) {
	goto L500;
    }
    if (iflgut == 3) {
	goto L600;
    }
/*                                        . . . . . . . . read parameters */
    mnread_((S_fp)fcn, &c__2, &iflgut, (U_fp)futil);
    if (iflgut == 2) {
	goto L500;
    }
    if (iflgut == 3) {
	goto L600;
    }
    if (iflgut == 4) {
	goto L700;
    }
/*                              . . . . . . verify FCN not time-dependent */
    io___10.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___10);
    do_fio(&c__1, " MINUIT: FIRST CALL TO USER FUNCTION,", (ftnlen)37);
    do_fio(&c__1, " WITH IFLAG=1", (ftnlen)13);
    e_wsfe();
    nparx = mn7npr_1.npar;
    mninex_(mn7int_1.x);
    fzero = mn7cns_1.undefi;
    (*fcn)(&nparx, mn7der_1.gin, &fzero, mn7ext_1.u, &c__1, (U_fp)futil);
    first = mn7cns_1.undefi;
    (*fcn)(&nparx, mn7der_1.gin, &first, mn7ext_1.u, &c__4, (U_fp)futil);
    mn7cnv_1.nfcn = 2;
    if (fzero == mn7cns_1.undefi && first == mn7cns_1.undefi) {
	s_copy(cwhyxt, "BY ERROR IN USER FUNCTION.  ", (ftnlen)40, (ftnlen)28)
		;
	io___14.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___14);
	do_fio(&c__1, " USER HAS NOT CALCULATED FUNCTION", (ftnlen)33);
	do_fio(&c__1, " VALUE WHEN IFLAG=1 OR 4", (ftnlen)24);
	e_wsfe();
	goto L800;
    }
    mn7min_1.amin = first;
    if (first == mn7cns_1.undefi) {
	mn7min_1.amin = fzero;
    }
    mnprin_(&c__1, &mn7min_1.amin);
    mn7cnv_1.nfcn = 2;
    if (first == fzero) {
	goto L300;
    }
    fnew = (float)0.;
    (*fcn)(&nparx, mn7der_1.gin, &fnew, mn7ext_1.u, &c__4, (U_fp)futil);
    if (fnew != mn7min_1.amin) {
	io___16.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___16);
	do_fio(&c__1, (char *)&mn7min_1.amin, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&fnew, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    mn7cnv_1.nfcn = 3;
L300:
    mn7min_1.fval3 = mn7min_1.amin * (float)2. + (float)1.;
/*                                   . . . . . . . . . . . read commands */
    mnread_((S_fp)fcn, &c__3, &iflgut, (U_fp)futil);
    if (iflgut == 2) {
	goto L500;
    }
    if (iflgut == 3) {
	goto L600;
    }
    if (iflgut == 4) {
	goto L700;
    }
/* Writing concatenation */
    i__1[0] = 19, a__1[0] = "BY MINUIT COMMAND: ";
    i__1[1] = 20, a__1[1] = mn7tit_1.cword;
    s_cat(cwhyxt, a__1, i__1, &c__2, (ftnlen)40);
    if (i_indx(mn7tit_1.cword, "STOP", (ftnlen)20, (ftnlen)4) > 0) {
	goto L800;
    }
    if (i_indx(mn7tit_1.cword, "EXI", (ftnlen)20, (ftnlen)3) > 0) {
	goto L800;
    }
    if (i_indx(mn7tit_1.cword, "RET", (ftnlen)20, (ftnlen)3) == 0) {
	goto L100;
    }
    s_copy(cwhyxt, "AND RETURNS TO USER PROGRAM.    ", (ftnlen)40, (ftnlen)32)
	    ;
    io___17.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___17);
    do_fio(&c__1, " ..........MINUIT TERMINATED ", (ftnlen)29);
    do_fio(&c__1, cwhyxt, (ftnlen)40);
    e_wsfe();
    return 0;
/*                                           . . . . . . stop conditions */
L500:
    s_copy(cwhyxt, "BY END-OF-DATA ON PRIMARY INPUT FILE.   ", (ftnlen)40, (
	    ftnlen)40);
    goto L800;
L600:
    s_copy(cwhyxt, "BY UNRECOVERABLE READ ERROR ON INPUT.   ", (ftnlen)40, (
	    ftnlen)40);
    goto L800;
L700:
    s_copy(cwhyxt, ": FATAL ERROR IN PARAMETER DEFINITIONS. ", (ftnlen)40, (
	    ftnlen)40);
L800:
    io___18.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___18);
    do_fio(&c__1, " ..........MINUIT TERMINATED ", (ftnlen)29);
    do_fio(&c__1, cwhyxt, (ftnlen)40);
    e_wsfe();
    s_stop("", (ftnlen)0);

/*  ......................entry to set unit numbers  - - - - - - - - - - */

L_mintio:
    jsysrd = *i1;
    jsyswr = *i2;
    jsyssa = *i3;
    return 0;
} /* minuit_ */

/* Subroutine */ int minuit_(S_fp fcn, U_fp futil)
{
    return minuit_0_(0, fcn, futil, (integer *)0, (integer *)0, (integer *)0);
    }

/* Subroutine */ int mintio_(integer *i1, integer *i2, integer *i3)
{
    return minuit_0_(1, (S_fp)0, (U_fp)0, i1, i2, i3);
    }


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:18  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni */
/* Minuit */


/* Subroutine */ int mnamin_(S_fp fcn, U_fp futil)
{
    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static doublereal fnew;
    static integer nparx;
    extern /* Subroutine */ int mnexin_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___20 = { 0, 0, 0, "(/A,A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called  from many places.  Initializes the value of AMIN by */
/* C        calling the user function. Prints out the function value and */
/* C        parameter values if Print Flag value is high enough. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    nparx = mn7npr_1.npar;
    if (mn7flg_1.isw[4] >= 1) {
	io___20.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___20);
	do_fio(&c__1, " FIRST CALL TO ", (ftnlen)15);
	do_fio(&c__1, "USER FUNCTION AT NEW START POINT, WITH IFLAG=4.", (
		ftnlen)47);
	e_wsfe();
    }
    mnexin_(mn7int_1.x);
    (*fcn)(&nparx, mn7der_1.gin, &fnew, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    mn7min_1.amin = fnew;
    mn7min_1.edm = mn7cns_1.bigedm;
    return 0;
} /* mnamin_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:18  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni */
/* Minuit */


/* Subroutine */ int mnbins_(doublereal *a1, doublereal *a2, integer *naa, 
	doublereal *bl, doublereal *bh, integer *nb, doublereal *bwid)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_lg10(doublereal *), pow_ri(real *, integer *);

    /* Local variables */
    static doublereal ah, al;
    static integer na;
    static doublereal alb;
    static integer log__;
    static doublereal awid;
    static integer kwid, lwid;
    static doublereal sigfig, sigrnd;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/*         SUBROUTINE TO DETERMINE REASONABLE HISTOGRAM INTERVALS */
/*         GIVEN ABSOLUTE UPPER AND LOWER BOUNDS  A1 AND A2 */
/*         AND DESIRED MAXIMUM NUMBER OF BINS NAA */
/*         PROGRAM MAKES REASONABLE BINNING FROM BL TO BH OF WIDTH BWID */
/*         F. JAMES,   AUGUST, 1974 , stolen for Minuit, 1988 */
    al = min(*a1,*a2);
    ah = max(*a1,*a2);
    if (al == ah) {
	ah = al + (float)1.;
    }
/*         IF NAA .EQ. -1 , PROGRAM USES BWID INPUT FROM CALLING ROUTINE */
    if (*naa == -1) {
	goto L150;
    }
L10:
    na = *naa - 1;
    if (na < 1) {
	na = 1;
    }
/*          GET NOMINAL BIN WIDTH IN EXPON FORM */
L20:
    awid = (ah - al) / (real) na;
    d__1 = awid;
    log__ = (integer) d_lg10(&d__1);
    if (awid <= 1.) {
	--log__;
    }
    i__1 = -log__;
    sigfig = awid * pow_ri(&c_b66, &i__1);
/*         ROUND MANTISSA UP TO 2, 2.5, 5, OR 10 */
    if (sigfig > (float)2.) {
	goto L40;
    }
    sigrnd = (float)2.;
    goto L100;
L40:
    if (sigfig > (float)2.5) {
	goto L50;
    }
    sigrnd = (float)2.5;
    goto L100;
L50:
    if (sigfig > (float)5.) {
	goto L60;
    }
    sigrnd = (float)5.;
    goto L100;
L60:
    sigrnd = (float)1.;
    ++log__;
L100:
    *bwid = sigrnd * pow_ri(&c_b66, &log__);
    goto L200;
/*         GET NEW BOUNDS FROM NEW WIDTH BWID */
L150:
    if (*bwid <= 0.) {
	goto L10;
    }
L200:
    alb = al / *bwid;
    lwid = (integer) alb;
    if (alb < 0.) {
	--lwid;
    }
    *bl = *bwid * (real) lwid;
    alb = ah / *bwid + (float)1.;
    kwid = (integer) alb;
    if (alb < 0.) {
	--kwid;
    }
    *bh = *bwid * (real) kwid;
    *nb = kwid - lwid;
    if (*naa > 5) {
	goto L240;
    }
    if (*naa == -1) {
	return 0;
    }
/*          REQUEST FOR ONE BIN IS DIFFICULT CASE */
    if (*naa > 1 || *nb == 1) {
	return 0;
    }
    *bwid *= (float)2.;
    *nb = 1;
    return 0;
L240:
    if (*nb << 1 != *naa) {
	return 0;
    }
    ++na;
    goto L20;
} /* mnbins_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni */
/* Minuit */


/* Subroutine */ int mncalf_(S_fp fcn, doublereal *pvec, doublereal *ycalf, 
	U_fp futil)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal f;
    static integer i__, j, m, n, ndex;
    static doublereal denom;
    static integer nparx;
    extern /* Subroutine */ int mninex_(doublereal *);


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called only from MNIMPR.  Transforms the function FCN */
/* C        by dividing out the quadratic part in order to find further */
/* C        minima.    Calculates  ycalf = (f-fmin)/(x-xmin)*v*(x-xmin) */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    /* Parameter adjustments */
    --pvec;

    /* Function Body */
    nparx = mn7npr_1.npar;
    mninex_(&pvec[1]);
    (*fcn)(&nparx, mn7der_1.gin, &f, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7der_1.grd[i__ - 1] = (float)0.;
	i__2 = mn7npr_1.npar;
	for (j = 1; j <= i__2; ++j) {
	    m = max(i__,j);
	    n = min(i__,j);
	    ndex = m * (m - 1) / 2 + n;
/* L200: */
	    mn7der_1.grd[i__ - 1] += mn7vat_1.vthmat[ndex - 1] * (mn7int_1.xt[
		    j - 1] - pvec[j]);
	}
    }
    denom = (float)0.;
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L210: */
	denom += mn7der_1.grd[i__ - 1] * (mn7int_1.xt[i__ - 1] - pvec[i__]);
    }
    if (denom <= 0.) {
	mn7min_1.dcovar = (float)1.;
	mn7flg_1.isw[1] = 0;
	denom = (float)1.;
    }
    *ycalf = (f - mn7min_1.apsi) / denom;
    return 0;
} /* mncalf_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni */
/* Minuit */


/* Subroutine */ int mncler_()
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int mnrset_(integer *);


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called from MINUIT and by option from MNEXCM */
/* C        Resets the parameter list to UNDEFINED */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    mn7fx1_1.npfix = 0;
    mn7npr_1.nu = 0;
    mn7npr_1.npar = 0;
    mn7cnv_1.nfcn = 0;
    mn7cnv_1.nwrmes[0] = 0;
    mn7cnv_1.nwrmes[1] = 0;
    i__1 = mn7npr_1.maxext;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7ext_1.u[i__ - 1] = (float)0.;
	s_copy(mn7nam_1.cpnam + (i__ - 1) * 10, mn7tit_1.cundef, (ftnlen)10, (
		ftnlen)10);
	mn7inx_1.nvarl[i__ - 1] = -1;
/* L10: */
	mn7inx_1.niofex[i__ - 1] = 0;
    }
    mnrset_(&c__1);
    s_copy(mn7tit_1.cfrom, "CLEAR   ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "UNDEFINED ", (ftnlen)10, (ftnlen)10);
    mn7log_1.lnolim = TRUE_;
    mn7log_1.lphead = TRUE_;
    return 0;
} /* mncler_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni */
/* Minuit */


/* Subroutine */ int mncntr_(S_fp fcn, integer *ke1, integer *ke2, integer *
	ierrf, U_fp futil)
{
    /* Initialized data */

    static char clabel[20+1] = "0123456789ABCDEFGHIJ";

    /* Format strings */
    static char fmt_1351[] = "(\002 INVALID PARAMETER NUMBER(S) REQUESTED.  \
IGNORED.\002/)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal ff;
    static integer nl, ix, iy, nx, ny, ki1, ki2, nl2;
    static doublereal xb4;
    static integer ics;
    static doublereal fmn, fmx, xlo, ylo, xup, yup, fcna[115], fcnb[115];
    static char chln[115];
    static doublereal devs, xsav, ysav;
    static char chmid[115];
    static integer ngrid, ixmid;
    static doublereal bwidx;
    static integer nparx;
    static doublereal bwidy, unext, ylabel;
    extern /* Subroutine */ int mnamin_(S_fp, U_fp);
    static doublereal contur[20];
    static char chzero[115];
    extern /* Subroutine */ int mnhess_(S_fp, U_fp), mnwerr_();
    static integer ixzero;

    /* Fortran I/O blocks */
    static cilist io___67 = { 0, 0, 0, "(A,I3,A,A)", 0 };
    static cilist io___69 = { 0, 0, 0, "(12X,A,A)", 0 };
    static cilist io___77 = { 0, 0, 0, "(1X,G12.4,1X,A)", 0 };
    static cilist io___78 = { 0, 0, 0, "(14X,A)", 0 };
    static cilist io___81 = { 0, 0, 0, "(8X,G12.4,A,G12.4)", 0 };
    static cilist io___82 = { 0, 0, 0, "(14X,A,G12.4)", 0 };
    static cilist io___83 = { 0, 0, 0, "(8X,G12.4,A,G12.4,A,G12.4)", 0 };
    static cilist io___84 = { 0, 0, 0, "(6X,A,I3,A,A,A,G12.4)", 0 };
    static cilist io___85 = { 0, 0, 0, "(A,G12.4,A,G12.4,A)", 0 };
    static cilist io___86 = { 0, 0, 0, fmt_1351, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       to print function contours in two variables, on line printer */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/*                 input arguments: parx, pary, devs, ngrid */
    if (*ke1 <= 0 || *ke2 <= 0) {
	goto L1350;
    }
    if (*ke1 > mn7npr_1.nu || *ke2 > mn7npr_1.nu) {
	goto L1350;
    }
    ki1 = mn7inx_1.niofex[*ke1 - 1];
    ki2 = mn7inx_1.niofex[*ke2 - 1];
    if (ki1 <= 0 || ki2 <= 0) {
	goto L1350;
    }
    if (ki1 == ki2) {
	goto L1350;
    }

    if (mn7flg_1.isw[1] < 1) {
	mnhess_((S_fp)fcn, (U_fp)futil);
	mnwerr_();
    }
    nparx = mn7npr_1.npar;
    xsav = mn7ext_1.u[*ke1 - 1];
    ysav = mn7ext_1.u[*ke2 - 1];
    devs = mn7arg_1.word7[2];
    if (devs <= 0.) {
	devs = (float)2.;
    }
    xlo = mn7ext_1.u[*ke1 - 1] - devs * mn7err_1.werr[ki1 - 1];
    xup = mn7ext_1.u[*ke1 - 1] + devs * mn7err_1.werr[ki1 - 1];
    ylo = mn7ext_1.u[*ke2 - 1] - devs * mn7err_1.werr[ki2 - 1];
    yup = mn7ext_1.u[*ke2 - 1] + devs * mn7err_1.werr[ki2 - 1];
    ngrid = (integer) mn7arg_1.word7[3];
    if (ngrid <= 0) {
	ngrid = 25;
/* Computing MIN */
	i__1 = mn7iou_1.npagwd - 15;
	nx = min(i__1,ngrid);
/* Computing MIN */
	i__1 = mn7iou_1.npagln - 7;
	ny = min(i__1,ngrid);
    } else {
	nx = ngrid;
	ny = ngrid;
    }
    if (nx < 11) {
	nx = 11;
    }
    if (ny < 11) {
	ny = 11;
    }
    if (nx >= 115) {
	nx = 114;
    }
/*         ask if parameter outside limits */
    if (mn7inx_1.nvarl[*ke1 - 1] > 1) {
	if (xlo < mn7ext_1.alim[*ke1 - 1]) {
	    xlo = mn7ext_1.alim[*ke1 - 1];
	}
	if (xup > mn7ext_1.blim[*ke1 - 1]) {
	    xup = mn7ext_1.blim[*ke1 - 1];
	}
    }
    if (mn7inx_1.nvarl[*ke2 - 1] > 1) {
	if (ylo < mn7ext_1.alim[*ke2 - 1]) {
	    ylo = mn7ext_1.alim[*ke2 - 1];
	}
	if (yup > mn7ext_1.blim[*ke2 - 1]) {
	    yup = mn7ext_1.blim[*ke2 - 1];
	}
    }
    bwidx = (xup - xlo) / (real) nx;
    bwidy = (yup - ylo) / (real) ny;
    ixmid = (integer) ((xsav - xlo) * (real) nx / (xup - xlo)) + 1;
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    for (i__ = 1; i__ <= 20; ++i__) {
/* Computing 2nd power */
	r__1 = (real) (i__ - 1);
	contur[i__ - 1] = mn7min_1.amin + mn7min_1.up * (r__1 * r__1);
/* L185: */
    }
    contur[0] += mn7min_1.up * (float).01;
/*                fill FCNB to prepare first row, and find column zero */
    mn7ext_1.u[*ke2 - 1] = yup;
    ixzero = 0;
    xb4 = 1.;
    i__1 = nx + 1;
    for (ix = 1; ix <= i__1; ++ix) {
	mn7ext_1.u[*ke1 - 1] = xlo + (real) (ix - 1) * bwidx;
	(*fcn)(&nparx, mn7der_1.gin, &ff, mn7ext_1.u, &c__4, (U_fp)futil);
	fcnb[ix - 1] = ff;
	if (xb4 < 0. && mn7ext_1.u[*ke1 - 1] > 0.) {
	    ixzero = ix - 1;
	}
	xb4 = mn7ext_1.u[*ke1 - 1];
	*(unsigned char *)&chmid[ix - 1] = '*';
	*(unsigned char *)&chzero[ix - 1] = '-';
/* L200: */
    }
    io___67.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___67);
    do_fio(&c__1, " Y-AXIS: PARAMETER ", (ftnlen)19);
    do_fio(&c__1, (char *)&(*ke2), (ftnlen)sizeof(integer));
    do_fio(&c__1, ": ", (ftnlen)2);
    do_fio(&c__1, mn7nam_1.cpnam + (*ke2 - 1) * 10, (ftnlen)10);
    e_wsfe();
    if (ixzero > 0) {
	*(unsigned char *)&chzero[ixzero - 1] = '+';
	s_copy(chln, " ", (ftnlen)115, (ftnlen)1);
	io___69.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___69);
	do_fio(&c__1, chln, ixzero);
	do_fio(&c__1, "X=0", (ftnlen)3);
	e_wsfe();
    }
/*                 loop over rows */
    i__1 = ny;
    for (iy = 1; iy <= i__1; ++iy) {
	unext = mn7ext_1.u[*ke2 - 1] - bwidy;
/*                 prepare this line's background pattern for contour */
	s_copy(chln, " ", (ftnlen)115, (ftnlen)1);
	*(unsigned char *)&chln[ixmid - 1] = '*';
	if (ixzero != 0) {
	    *(unsigned char *)&chln[ixzero - 1] = ':';
	}
	if (mn7ext_1.u[*ke2 - 1] > ysav && unext < ysav) {
	    s_copy(chln, chmid, (ftnlen)115, (ftnlen)115);
	}
	if (mn7ext_1.u[*ke2 - 1] > 0. && unext < 0.) {
	    s_copy(chln, chzero, (ftnlen)115, (ftnlen)115);
	}
	mn7ext_1.u[*ke2 - 1] = unext;
	ylabel = mn7ext_1.u[*ke2 - 1] + bwidy * (float).5;
/*                 move FCNB to FCNA and fill FCNB with next row */
	i__2 = nx + 1;
	for (ix = 1; ix <= i__2; ++ix) {
	    fcna[ix - 1] = fcnb[ix - 1];
	    mn7ext_1.u[*ke1 - 1] = xlo + (real) (ix - 1) * bwidx;
	    (*fcn)(&nparx, mn7der_1.gin, &ff, mn7ext_1.u, &c__4, (U_fp)futil);
	    fcnb[ix - 1] = ff;
/* L220: */
	}
/*                 look for contours crossing the FCNxy squares */
	i__2 = nx;
	for (ix = 1; ix <= i__2; ++ix) {
/* Computing MAX */
	    d__1 = fcna[ix - 1], d__2 = fcnb[ix - 1], d__1 = max(d__1,d__2), 
		    d__2 = fcna[ix], d__1 = max(d__1,d__2), d__2 = fcnb[ix];
	    fmx = max(d__1,d__2);
/* Computing MIN */
	    d__1 = fcna[ix - 1], d__2 = fcnb[ix - 1], d__1 = min(d__1,d__2), 
		    d__2 = fcna[ix], d__1 = min(d__1,d__2), d__2 = fcnb[ix];
	    fmn = min(d__1,d__2);
	    for (ics = 1; ics <= 20; ++ics) {
		if (contur[ics - 1] > fmn) {
		    goto L240;
		}
/* L230: */
	    }
	    goto L250;
L240:
	    if (contur[ics - 1] < fmx) {
		*(unsigned char *)&chln[ix - 1] = *(unsigned char *)&clabel[
			ics - 1];
	    }
L250:
	    ;
	}
/*                 print a row of the contour plot */
	io___77.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___77);
	do_fio(&c__1, (char *)&ylabel, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, chln, nx);
	e_wsfe();
/* L280: */
    }
/*                 contours printed, label x-axis */
    s_copy(chln, " ", (ftnlen)115, (ftnlen)1);
    *(unsigned char *)chln = 'I';
    *(unsigned char *)&chln[ixmid - 1] = 'I';
    *(unsigned char *)&chln[nx - 1] = 'I';
    io___78.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___78);
    do_fio(&c__1, chln, nx);
    e_wsfe();
/*                the hardest of all: print x-axis scale! */
    s_copy(chln, " ", (ftnlen)115, (ftnlen)1);
    if (nx <= 26) {
/* Computing MAX */
	i__1 = nx - 12;
	nl = max(i__1,2);
	nl2 = nl / 2;
	io___81.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___81);
	do_fio(&c__1, (char *)&xlo, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, chln, nl);
	do_fio(&c__1, (char *)&xup, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___82.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___82);
	do_fio(&c__1, chln, nl2);
	do_fio(&c__1, (char *)&xsav, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
/* Computing MAX */
	i__1 = nx - 24;
	nl = max(i__1,2) / 2;
	nl2 = nl;
	if (nl > 10) {
	    nl2 = nl - 6;
	}
	io___83.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___83);
	do_fio(&c__1, (char *)&xlo, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, chln, nl);
	do_fio(&c__1, (char *)&xsav, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, chln, nl2);
	do_fio(&c__1, (char *)&xup, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    io___84.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___84);
    do_fio(&c__1, " X-AXIS: PARAMETER", (ftnlen)18);
    do_fio(&c__1, (char *)&(*ke1), (ftnlen)sizeof(integer));
    do_fio(&c__1, ": ", (ftnlen)2);
    do_fio(&c__1, mn7nam_1.cpnam + (*ke1 - 1) * 10, (ftnlen)10);
    do_fio(&c__1, "  ONE COLUMN=", (ftnlen)13);
    do_fio(&c__1, (char *)&bwidx, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___85.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___85);
    do_fio(&c__1, " FUNCTION VALUES: F(I)=", (ftnlen)23);
    do_fio(&c__1, (char *)&mn7min_1.amin, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, " +", (ftnlen)2);
    do_fio(&c__1, (char *)&mn7min_1.up, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, " *I**2", (ftnlen)6);
    e_wsfe();
/*                 finished.  reset input values */
    mn7ext_1.u[*ke1 - 1] = xsav;
    mn7ext_1.u[*ke2 - 1] = ysav;
    *ierrf = 0;
    return 0;
L1350:
    io___86.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___86);
    e_wsfe();
    *ierrf = 1;
    return 0;
} /* mncntr_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mncomd_(U_fp fcn, char *crdbin, integer *icondn, U_fp 
	futil, ftnlen crdbin_len)
{
    /* Initialized data */

    static char clower[26+1] = "abcdefghijklmnopqrstuvwxyz";
    static char cupper[26+1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, ic, lnc, ierr, ipos, llist;
    static doublereal plist[30];
    static logical leader;
    static char comand[20], crdbuf[100];
    static integer lenbuf;
    extern /* Subroutine */ int mncrck_(char *, integer *, char *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, ftnlen, 
	    ftnlen), mnexcm_(U_fp, char *, doublereal *, integer *, integer *,
	     U_fp, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___95 = { 0, 0, 0, "(A)", 0 };
    static cilist io___101 = { 0, 0, 0, "(A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called by user.  'Reads' a command string and executes. */
/* C     Equivalent to MNEXCM except that the command is given as a */
/* C          character string. */
/* C */
/* C     ICONDN = 0: command executed normally */
/* C              1: command is blank, ignored */
/* C              2: command line unreadable, ignored */
/* C              3: unknown command, ignored */
/* C              4: abnormal termination (e.g., MIGRAD not converged) */
/* C              5: command is a request to read PARAMETER definitions */
/* C              6: 'SET INPUT' command */
/* C              7: 'SET TITLE' command */
/* C              8: 'SET COVAR' command */
/* C              9: reserved */
/* C             10: END command */
/* C             11: EXIT or STOP command */
/* C             12: RETURN command */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */




    lenbuf = i_len(crdbin, crdbin_len);
    s_copy(crdbuf, crdbin, (ftnlen)100, crdbin_len);
    *icondn = 0;
/*     record not case-sensitive, get upper case, strip leading blanks */
    leader = TRUE_;
    ipos = 1;
    i__1 = min(20,lenbuf);
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&crdbuf[i__ - 1] == '\'') {
	    goto L111;
	}
	if (*(unsigned char *)&crdbuf[i__ - 1] == ' ') {
	    if (leader) {
		++ipos;
	    }
	    goto L110;
	}
	leader = FALSE_;
	for (ic = 1; ic <= 26; ++ic) {
	    if (*(unsigned char *)&crdbuf[i__ - 1] == *(unsigned char *)&
		    clower[ic - 1]) {
		*(unsigned char *)&crdbuf[i__ - 1] = *(unsigned char *)&
			cupper[ic - 1];
	    }
/* L108: */
	}
L110:
	;
    }
L111:
/*                     blank or null command */
    if (ipos > lenbuf) {
	io___95.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___95);
	do_fio(&c__1, " BLANK COMMAND IGNORED.", (ftnlen)23);
	e_wsfe();
	*icondn = 1;
	goto L900;
    }
/*                                           . .   preemptive commands */
/*               if command is 'PARAMETER' */
    if (s_cmp(crdbuf + (ipos - 1), "PAR", (ftnlen)3, (ftnlen)3) == 0) {
	*icondn = 5;
	mn7log_1.lphead = TRUE_;
	goto L900;
    }
/*               if command is 'SET INPUT' */
    if (s_cmp(crdbuf + (ipos - 1), "SET INP", (ftnlen)7, (ftnlen)7) == 0) {
	*icondn = 6;
	mn7log_1.lphead = TRUE_;
	goto L900;
    }
/*              if command is 'SET TITLE' */
    if (s_cmp(crdbuf + (ipos - 1), "SET TIT", (ftnlen)7, (ftnlen)7) == 0) {
	*icondn = 7;
	mn7log_1.lphead = TRUE_;
	goto L900;
    }
/*               if command is 'SET COVARIANCE' */
    if (s_cmp(crdbuf + (ipos - 1), "SET COV", (ftnlen)7, (ftnlen)7) == 0) {
	*icondn = 8;
	mn7log_1.lphead = TRUE_;
	goto L900;
    }
/*               crack the command . . . . . . . . . . . . . . . . */
    mncrck_(crdbuf + (ipos - 1), &c__20, comand, &lnc, &c__30, plist, &llist, 
	    &ierr, &mn7iou_1.isyswr, lenbuf - (ipos - 1), (ftnlen)20);
    if (ierr > 0) {
	io___101.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___101);
	do_fio(&c__1, " COMMAND CANNOT BE INTERPRETED", (ftnlen)30);
	e_wsfe();
	*icondn = 2;
	goto L900;
    }

    mnexcm_((U_fp)fcn, comand, plist, &llist, &ierr, (U_fp)futil, lnc);
    *icondn = ierr;
L900:
    return 0;
} /* mncomd_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mncont_(U_fp fcn, integer *ke1, integer *ke2, integer *
	nptu, doublereal *xptu, doublereal *yptu, integer *ierrf, U_fp futil)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer i__, j;
    static doublereal w[100], a1, a2;
    static integer i1, i2;
    static doublereal dc;
    static integer lr, ke3, ki1, ki2, ki3;
    static doublereal gcc[100];
    static integer isw2, isw4, nall, iold, line, mpar, ierr, inew;
    static doublereal dist, xdir, ydir, aopt;
    static integer move, next;
    static doublereal u1min, u2min, abest;
    static integer nfcol, iercr;
    static doublereal scalx, scaly;
    static integer idist, npcol, kints;
    static doublereal val2mi, val2pl, sclfac, bigdis;
    static logical ldebug;
    static integer nfcnco;
    extern /* Subroutine */ int mnfree_(integer *), mncuve_(U_fp, U_fp), 
	    mnmnot_(U_fp, integer *, integer *, doublereal *, doublereal *, 
	    U_fp), mnwarn_(char *, char *, char *, ftnlen, ftnlen, ftnlen), 
	    mnplot_(doublereal *, doublereal *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen);
    static integer nowpts;
    static doublereal sigsav;
    static integer istrav, nfmxin;
    extern /* Subroutine */ int mnfixp_(integer *, integer *), mncros_(U_fp, 
	    doublereal *, integer *, U_fp), mninex_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___108 = { 0, 0, 0, "(1X,A,I4,A)", 0 };
    static cilist io___111 = { 0, 0, 0, "(1X,A,I3,2X,A)", 0 };
    static cilist io___112 = { 0, 0, 0, "(1X,A,I3,A)", 0 };
    static cilist io___119 = { 0, 0, 0, "(A)", 0 };
    static cilist io___149 = { 0, 0, 0, "(A,A,I3,A)", 0 };
    static cilist io___151 = { 0, 0, 0, "(A,I3,2X,A)", 0 };
    static cilist io___152 = { 0, 0, 0, "(25X,A,I3,2X,A)", 0 };
    static cilist io___155 = { 0, 0, 0, "(/I5,A,G13.5,A,G11.3)", 0 };
    static cilist io___156 = { 0, 0, 0, "(9X,A,3X,A,18X,A,3X,A)", 0 };
    static cilist io___159 = { 0, 0, 0, "(1X,I5,2G13.5,10X,I5,2G13.5)", 0 };
    static cilist io___160 = { 0, 0, 0, "(1X,I5,2G13.5)", 0 };
    static cilist io___161 = { 0, 0, 0, "(A)", 0 };
    static cilist io___162 = { 0, 0, 0, "(A)", 0 };
    static cilist io___163 = { 0, 0, 0, "(A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       Find NPTU points along a contour where the function */
/* C             FMIN (X(KE1),X(KE2)) =  AMIN+UP */
/* C       where FMIN is the minimum of FCN with respect to all */
/* C       the other NPAR-2 variable parameters (if any). */
/* C   IERRF on return will be equal to the number of points found: */
/* C     NPTU if normal termination with NPTU points found */
/* C     -1   if errors in the calling sequence (KE1, KE2 not variable) */
/* C      0   if less than four points can be found (using MNMNOT) */
/* C     n>3  if only n points can be found (n < NPTU) */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/*                 input arguments: parx, pary, devs, ngrid */
    /* Parameter adjustments */
    --yptu;
    --xptu;

    /* Function Body */
    ldebug = mn7flg_1.idbg[6] >= 1;
    if (*ke1 <= 0 || *ke2 <= 0) {
	goto L1350;
    }
    if (*ke1 > mn7npr_1.nu || *ke2 > mn7npr_1.nu) {
	goto L1350;
    }
    ki1 = mn7inx_1.niofex[*ke1 - 1];
    ki2 = mn7inx_1.niofex[*ke2 - 1];
    if (ki1 <= 0 || ki2 <= 0) {
	goto L1350;
    }
    if (ki1 == ki2) {
	goto L1350;
    }
    if (*nptu < 4) {
	goto L1400;
    }

    nfcnco = mn7cnv_1.nfcn;
    mn7cnv_1.nfcnmx = (*nptu + 5) * 100 * (mn7npr_1.npar + 1);
/*           The minimum */
    mncuve_((U_fp)fcn, (U_fp)futil);
    u1min = mn7ext_1.u[*ke1 - 1];
    u2min = mn7ext_1.u[*ke2 - 1];
    *ierrf = 0;
    s_copy(mn7tit_1.cfrom, "MNContour ", (ftnlen)8, (ftnlen)10);
    mn7cnv_1.nfcnfr = nfcnco;
    if (mn7flg_1.isw[4] >= 0) {
	io___108.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___108);
	do_fio(&c__1, "START MNCONTOUR CALCULATION OF", (ftnlen)30);
	do_fio(&c__1, (char *)&(*nptu), (ftnlen)sizeof(integer));
	do_fio(&c__1, " POINTS ON CONTOUR.", (ftnlen)19);
	e_wsfe();
	if (mn7npr_1.npar > 2) {
	    if (mn7npr_1.npar == 3) {
		ki3 = 6 - ki1 - ki2;
		ke3 = mn7inx_1.nexofi[ki3 - 1];
		io___111.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___111);
		do_fio(&c__1, "EACH POINT IS A MINIMUM WITH RESPECT TO PARAM\
ETER ", (ftnlen)50);
		do_fio(&c__1, (char *)&ke3, (ftnlen)sizeof(integer));
		do_fio(&c__1, mn7nam_1.cpnam + (ke3 - 1) * 10, (ftnlen)10);
		e_wsfe();
	    } else {
		io___112.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___112);
		do_fio(&c__1, "EACH POINT IS A MINIMUM WITH RESPECT TO THE O\
THER", (ftnlen)49);
		i__1 = mn7npr_1.npar - 2;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		do_fio(&c__1, " VARIABLE PARAMETERS.", (ftnlen)21);
		e_wsfe();
	    }
	}
    }

/*           Find the first four points using MNMNOT */
/*              ........................ first two points */
    mnmnot_((U_fp)fcn, ke1, ke2, &val2pl, &val2mi, (U_fp)futil);
    if (mn7err_1.ern[ki1 - 1] == mn7cns_1.undefi) {
	xptu[1] = mn7ext_1.alim[*ke1 - 1];
	mnwarn_("W", "MNContour ", "Contour squeezed by parameter limits.", (
		ftnlen)1, (ftnlen)10, (ftnlen)37);
    } else {
	if (mn7err_1.ern[ki1 - 1] >= 0.) {
	    goto L1500;
	}
	xptu[1] = u1min + mn7err_1.ern[ki1 - 1];
    }
    yptu[1] = val2mi;

    if (mn7err_1.erp[ki1 - 1] == mn7cns_1.undefi) {
	xptu[3] = mn7ext_1.blim[*ke1 - 1];
	mnwarn_("W", "MNContour ", "Contour squeezed by parameter limits.", (
		ftnlen)1, (ftnlen)10, (ftnlen)37);
    } else {
	if (mn7err_1.erp[ki1 - 1] <= 0.) {
	    goto L1500;
	}
	xptu[3] = u1min + mn7err_1.erp[ki1 - 1];
    }
    yptu[3] = val2pl;
    scalx = (float)1. / (xptu[3] - xptu[1]);
/*              ........................... next two points */
    mnmnot_((U_fp)fcn, ke2, ke1, &val2pl, &val2mi, (U_fp)futil);
    if (mn7err_1.ern[ki2 - 1] == mn7cns_1.undefi) {
	yptu[2] = mn7ext_1.alim[*ke2 - 1];
	mnwarn_("W", "MNContour ", "Contour squeezed by parameter limits.", (
		ftnlen)1, (ftnlen)10, (ftnlen)37);
    } else {
	if (mn7err_1.ern[ki2 - 1] >= 0.) {
	    goto L1500;
	}
	yptu[2] = u2min + mn7err_1.ern[ki2 - 1];
    }
    xptu[2] = val2mi;
    if (mn7err_1.erp[ki2 - 1] == mn7cns_1.undefi) {
	yptu[4] = mn7ext_1.blim[*ke2 - 1];
	mnwarn_("W", "MNContour ", "Contour squeezed by parameter limits.", (
		ftnlen)1, (ftnlen)10, (ftnlen)37);
    } else {
	if (mn7err_1.erp[ki2 - 1] <= 0.) {
	    goto L1500;
	}
	yptu[4] = u2min + mn7err_1.erp[ki2 - 1];
    }
    xptu[4] = val2pl;
    scaly = (float)1. / (yptu[4] - yptu[2]);
    nowpts = 4;
    next = 5;
    if (ldebug) {
	io___119.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___119);
	do_fio(&c__1, " Plot of four points found by MINOS", (ftnlen)35);
	e_wsfe();
	mn7rpt_1.xpt[0] = u1min;
	mn7rpt_1.ypt[0] = u2min;
	*(unsigned char *)&mn7cpt_1.chpt[0] = ' ';
/* Computing MIN */
	i__1 = nowpts + 1;
	nall = min(i__1,101);
	i__1 = nall;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    mn7rpt_1.xpt[i__ - 1] = xptu[i__ - 1];
	    mn7rpt_1.ypt[i__ - 1] = yptu[i__ - 1];
/* L85: */
	}
	*(unsigned char *)&mn7cpt_1.chpt[1] = 'A';
	*(unsigned char *)&mn7cpt_1.chpt[2] = 'B';
	*(unsigned char *)&mn7cpt_1.chpt[3] = 'C';
	*(unsigned char *)&mn7cpt_1.chpt[4] = 'D';
	mnplot_(mn7rpt_1.xpt, mn7rpt_1.ypt, mn7cpt_1.chpt, &nall, &
		mn7iou_1.isyswr, &mn7iou_1.npagwd, &mn7iou_1.npagln, (ftnlen)
		1);
    }

/*               ..................... save some values before fixing */
    isw2 = mn7flg_1.isw[1];
    isw4 = mn7flg_1.isw[3];
    sigsav = mn7min_1.edm;
    istrav = mn7cnv_1.istrat;
    dc = mn7min_1.dcovar;
    mn7min_1.apsi = mn7min_1.epsi * (float).5;
    abest = mn7min_1.amin;
    mpar = mn7npr_1.npar;
    nfmxin = mn7cnv_1.nfcnmx;
    i__1 = mpar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L125: */
	mn7int_1.xt[i__ - 1] = mn7int_1.x[i__ - 1];
    }
    i__1 = mpar * (mpar + 1) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L130: */
	mn7vat_1.vthmat[j - 1] = mn7var_1.vhmat[j - 1];
    }
    i__1 = mpar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gcc[i__ - 1] = mn7err_1.globcc[i__ - 1];
/* L135: */
	w[i__ - 1] = mn7err_1.werr[i__ - 1];
    }
/*                           fix the two parameters in question */
    kints = mn7inx_1.niofex[*ke1 - 1];
    mnfixp_(&kints, &ierr);
    kints = mn7inx_1.niofex[*ke2 - 1];
    mnfixp_(&kints, &ierr);
/*               ......................Fill in the rest of the points */
    i__1 = *nptu;
    for (inew = next; inew <= i__1; ++inew) {
/*            find the two neighbouring points with largest separation */
	bigdis = (float)0.;
	i__2 = inew - 1;
	for (iold = 1; iold <= i__2; ++iold) {
	    i2 = iold + 1;
	    if (i2 == inew) {
		i2 = 1;
	    }
/* Computing 2nd power */
	    d__1 = scalx * (xptu[iold] - xptu[i2]);
/* Computing 2nd power */
	    d__2 = scaly * (yptu[iold] - yptu[i2]);
	    dist = d__1 * d__1 + d__2 * d__2;
	    if (dist > bigdis) {
		bigdis = dist;
		idist = iold;
	    }
/* L200: */
	}
	i1 = idist;
	i2 = i1 + 1;
	if (i2 == inew) {
	    i2 = 1;
	}
/*                   next point goes between I1 and I2 */
	a1 = .5;
	a2 = .5;
L300:
	mn7xcr_1.xmidcr = a1 * xptu[i1] + a2 * xptu[i2];
	mn7xcr_1.ymidcr = a1 * yptu[i1] + a2 * yptu[i2];
	xdir = yptu[i2] - yptu[i1];
	ydir = xptu[i1] - xptu[i2];
/* Computing MAX */
	d__3 = (d__1 = xdir * scalx, abs(d__1)), d__4 = (d__2 = ydir * scaly, 
		abs(d__2));
	sclfac = max(d__3,d__4);
	mn7xcr_1.xdircr = xdir / sclfac;
	mn7xcr_1.ydircr = ydir / sclfac;
	mn7xcr_1.ke1cr = *ke1;
	mn7xcr_1.ke2cr = *ke2;
/*                Find the contour crossing point along DIR */
	mn7min_1.amin = abest;
	mncros_((U_fp)fcn, &aopt, &iercr, (U_fp)futil);
	if (iercr > 1) {
/*              If cannot find mid-point, try closer to point 1 */
	    if (a1 > .5) {
		if (mn7flg_1.isw[4] >= 0) {
		    io___149.ciunit = mn7iou_1.isyswr;
		    s_wsfe(&io___149);
		    do_fio(&c__1, " MNCONT CANNOT FIND NEXT", (ftnlen)24);
		    do_fio(&c__1, " POINT ON CONTOUR.  ONLY ", (ftnlen)25);
		    do_fio(&c__1, (char *)&nowpts, (ftnlen)sizeof(integer));
		    do_fio(&c__1, " POINTS FOUND.", (ftnlen)14);
		    e_wsfe();
		}
		goto L950;
	    }
	    mnwarn_("W", "MNContour ", "Cannot find midpoint, try closer.", (
		    ftnlen)1, (ftnlen)10, (ftnlen)33);
	    a1 = (float).75;
	    a2 = (float).25;
	    goto L300;
	}
/*                Contour has been located, insert new point in list */
	i__2 = i1 + 1;
	for (move = nowpts; move >= i__2; --move) {
	    xptu[move + 1] = xptu[move];
	    yptu[move + 1] = yptu[move];
/* L830: */
	}
	++nowpts;
	xptu[i1 + 1] = mn7xcr_1.xmidcr + mn7xcr_1.xdircr * aopt;
	yptu[i1 + 1] = mn7xcr_1.ymidcr + mn7xcr_1.ydircr * aopt;
/* L900: */
    }
L950:

    *ierrf = nowpts;
    s_copy(mn7tit_1.cstatu, "SUCCESSFUL", (ftnlen)10, (ftnlen)10);
    if (nowpts < *nptu) {
	s_copy(mn7tit_1.cstatu, "INCOMPLETE", (ftnlen)10, (ftnlen)10);
    }
/*                make a lineprinter plot of the contour */
    if (mn7flg_1.isw[4] >= 0) {
	mn7rpt_1.xpt[0] = u1min;
	mn7rpt_1.ypt[0] = u2min;
	*(unsigned char *)&mn7cpt_1.chpt[0] = ' ';
/* Computing MIN */
	i__1 = nowpts + 1;
	nall = min(i__1,101);
	i__1 = nall;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    mn7rpt_1.xpt[i__ - 1] = xptu[i__ - 1];
	    mn7rpt_1.ypt[i__ - 1] = yptu[i__ - 1];
	    *(unsigned char *)&mn7cpt_1.chpt[i__ - 1] = 'X';
/* L1000: */
	}
	io___151.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___151);
	do_fio(&c__1, " Y-AXIS: PARAMETER ", (ftnlen)19);
	do_fio(&c__1, (char *)&(*ke2), (ftnlen)sizeof(integer));
	do_fio(&c__1, mn7nam_1.cpnam + (*ke2 - 1) * 10, (ftnlen)10);
	e_wsfe();
	mnplot_(mn7rpt_1.xpt, mn7rpt_1.ypt, mn7cpt_1.chpt, &nall, &
		mn7iou_1.isyswr, &mn7iou_1.npagwd, &mn7iou_1.npagln, (ftnlen)
		1);
	io___152.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___152);
	do_fio(&c__1, "X-AXIS: PARAMETER ", (ftnlen)18);
	do_fio(&c__1, (char *)&(*ke1), (ftnlen)sizeof(integer));
	do_fio(&c__1, mn7nam_1.cpnam + (*ke1 - 1) * 10, (ftnlen)10);
	e_wsfe();
    }
/*                 print out the coordinates around the contour */
    if (mn7flg_1.isw[4] >= 1) {
	npcol = (nowpts + 1) / 2;
	nfcol = nowpts / 2;
	io___155.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___155);
	do_fio(&c__1, (char *)&nowpts, (ftnlen)sizeof(integer));
	do_fio(&c__1, " POINTS ON CONTOUR.   FMIN=", (ftnlen)27);
	do_fio(&c__1, (char *)&abest, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "   ERRDEF=", (ftnlen)10);
	do_fio(&c__1, (char *)&mn7min_1.up, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___156.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___156);
	do_fio(&c__1, mn7nam_1.cpnam + (*ke1 - 1) * 10, (ftnlen)10);
	do_fio(&c__1, mn7nam_1.cpnam + (*ke2 - 1) * 10, (ftnlen)10);
	do_fio(&c__1, mn7nam_1.cpnam + (*ke1 - 1) * 10, (ftnlen)10);
	do_fio(&c__1, mn7nam_1.cpnam + (*ke2 - 1) * 10, (ftnlen)10);
	e_wsfe();
	i__1 = nfcol;
	for (line = 1; line <= i__1; ++line) {
	    lr = line + npcol;
	    io___159.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___159);
	    do_fio(&c__1, (char *)&line, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xptu[line], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&yptu[line], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&lr, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xptu[lr], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&yptu[lr], (ftnlen)sizeof(doublereal));
	    e_wsfe();
/* L1050: */
	}
	if (nfcol < npcol) {
	    io___160.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___160);
	    do_fio(&c__1, (char *)&npcol, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xptu[npcol], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&yptu[npcol], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/*                                    . . contour finished. reset v */
    mn7cnv_1.itaur = 1;
    mnfree_(&c__1);
    mnfree_(&c__1);
    i__1 = mpar * (mpar + 1) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L1100: */
	mn7var_1.vhmat[j - 1] = mn7vat_1.vthmat[j - 1];
    }
    i__1 = mpar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7err_1.globcc[i__ - 1] = gcc[i__ - 1];
	mn7err_1.werr[i__ - 1] = w[i__ - 1];
/* L1120: */
	mn7int_1.x[i__ - 1] = mn7int_1.xt[i__ - 1];
    }
    mninex_(mn7int_1.x);
    mn7min_1.edm = sigsav;
    mn7min_1.amin = abest;
    mn7flg_1.isw[1] = isw2;
    mn7flg_1.isw[3] = isw4;
    mn7min_1.dcovar = dc;
    mn7cnv_1.itaur = 0;
    mn7cnv_1.nfcnmx = nfmxin;
    mn7cnv_1.istrat = istrav;
    mn7ext_1.u[*ke1 - 1] = u1min;
    mn7ext_1.u[*ke2 - 1] = u2min;
    goto L2000;
/*                                     Error returns */
L1350:
    io___161.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___161);
    do_fio(&c__1, " INVALID PARAMETER NUMBERS.", (ftnlen)27);
    e_wsfe();
    goto L1450;
L1400:
    io___162.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___162);
    do_fio(&c__1, " LESS THAN FOUR POINTS REQUESTED.", (ftnlen)33);
    e_wsfe();
L1450:
    *ierrf = -1;
    s_copy(mn7tit_1.cstatu, "USER ERROR", (ftnlen)10, (ftnlen)10);
    goto L2000;
L1500:
    io___163.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___163);
    do_fio(&c__1, " MNCONT UNABLE TO FIND FOUR POINTS.", (ftnlen)35);
    e_wsfe();
    mn7ext_1.u[*ke1 - 1] = u1min;
    mn7ext_1.u[*ke2 - 1] = u2min;
    *ierrf = 0;
    s_copy(mn7tit_1.cstatu, "FAILED", (ftnlen)10, (ftnlen)6);
L2000:
    s_copy(mn7tit_1.cfrom, "MNContour ", (ftnlen)8, (ftnlen)10);
    mn7cnv_1.nfcnfr = nfcnco;
    return 0;
} /* mncont_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mncrck_(char *crdbuf, integer *maxcwd, char *comand, 
	integer *lnc, integer *mxp, doublereal *plist, integer *llist, 
	integer *ierr, integer *isyswr, ftnlen crdbuf_len, ftnlen comand_len)
{
    /* Initialized data */

    static char cnull[15+1] = ")NULL STRING   ";
    static char cnumer[13+1] = "123456789-.0+";

    /* Format strings */
    static char fmt_253[] = "(\002 MINUIT WARNING: INPUT DATA WORD TOO LONG\
.\002/\002     ORIGINAL:\002,a/\002 TRUNCATED TO:\002,a)";
    static char fmt_511[] = "(/\002 MINUIT WARNING IN MNCRCK: \002/\002 COMM\
AND HAS INPUT\002,i5,\002 NUMERIC FIELDS, BUT MINUIT CAN ACCEPT ONLY\002,i3)";

    /* System generated locals */
    integer i__1, i__2;
    icilist ici__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_cmp(char *, char *, ftnlen, ftnlen), s_rsfi(icilist *), e_rsfi()
	    ;

    /* Local variables */
    static integer ic, ifld, iend, lend, left, nreq, ipos, kcmnd, nextb, 
	    ibegin, ltoadd;
    static char celmnt[19*25];
    static integer ielmnt, lelmnt[25], nelmnt;

    /* Fortran I/O blocks */
    static cilist io___174 = { 0, 0, 0, fmt_253, 0 };
    static cilist io___182 = { 0, 0, 0, fmt_511, 0 };
    static cilist io___183 = { 0, 0, 0, "(A,A,A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C */
/* C       Called from MNREAD. */
/* C       Cracks the free-format input, expecting zero or more */
/* C         alphanumeric fields (which it joins into COMAND(1:LNC)) */
/* C         followed by one or more numeric fields separated by */
/* C         blanks and/or one comma.  The numeric fields are put into */
/* C         the LLIST (but at most MXP) elements of PLIST. */
/* C      IERR = 0 if no errors, */
/* C           = 1 if error(s). */
/* C      Diagnostic messages are written to ISYSWR */
/* C */
    /* Parameter adjustments */
    --plist;

    /* Function Body */
    ielmnt = 0;
    lend = i_len(crdbuf, crdbuf_len);
    nextb = 1;
    *ierr = 0;
/*                                   . . . .  loop over words CELMNT */
L10:
    i__1 = lend;
    for (ipos = nextb; ipos <= i__1; ++ipos) {
	ibegin = ipos;
	if (*(unsigned char *)&crdbuf[ipos - 1] == ' ') {
	    goto L100;
	}
	if (*(unsigned char *)&crdbuf[ipos - 1] == ',') {
	    goto L250;
	}
	goto L150;
L100:
	;
    }
    goto L300;
L150:
/*               found beginning of word, look for end */
    i__1 = lend;
    for (ipos = ibegin + 1; ipos <= i__1; ++ipos) {
	if (*(unsigned char *)&crdbuf[ipos - 1] == ' ') {
	    goto L250;
	}
	if (*(unsigned char *)&crdbuf[ipos - 1] == ',') {
	    goto L250;
	}
/* L180: */
    }
    ipos = lend + 1;
L250:
    iend = ipos - 1;
    ++ielmnt;
    if (iend >= ibegin) {
	s_copy(celmnt + (ielmnt - 1) * 19, crdbuf + (ibegin - 1), (ftnlen)19, 
		iend - (ibegin - 1));
    } else {
	s_copy(celmnt + (ielmnt - 1) * 19, cnull, (ftnlen)19, (ftnlen)15);
    }
    lelmnt[ielmnt - 1] = iend - ibegin + 1;
    if (lelmnt[ielmnt - 1] > 19) {
	io___174.ciunit = *isyswr;
	s_wsfe(&io___174);
	do_fio(&c__1, crdbuf + (ibegin - 1), iend - (ibegin - 1));
	do_fio(&c__1, celmnt + (ielmnt - 1) * 19, (ftnlen)19);
	e_wsfe();
	lelmnt[ielmnt - 1] = 19;
    }
    if (ipos >= lend) {
	goto L300;
    }
    if (ielmnt >= 25) {
	goto L300;
    }
/*                     look for comma or beginning of next word */
    i__1 = lend;
    for (ipos = iend + 1; ipos <= i__1; ++ipos) {
	if (*(unsigned char *)&crdbuf[ipos - 1] == ' ') {
	    goto L280;
	}
	nextb = ipos;
	if (*(unsigned char *)&crdbuf[ipos - 1] == ',') {
	    nextb = ipos + 1;
	}
	goto L10;
L280:
	;
    }
/*                 All elements found, join the alphabetic ones to */
/*                                form a command */
L300:
    nelmnt = ielmnt;
    s_copy(comand, " ", comand_len, (ftnlen)1);
    *lnc = 1;
    plist[1] = (float)0.;
    *llist = 0;
    if (ielmnt == 0) {
	goto L900;
    }
    kcmnd = 0;
    i__1 = nelmnt;
    for (ielmnt = 1; ielmnt <= i__1; ++ielmnt) {
	if (s_cmp(celmnt + (ielmnt - 1) * 19, cnull, (ftnlen)19, (ftnlen)15) 
		== 0) {
	    goto L450;
	}
	for (ic = 1; ic <= 13; ++ic) {
	    if (*(unsigned char *)&celmnt[(ielmnt - 1) * 19] == *(unsigned 
		    char *)&cnumer[ic - 1]) {
		goto L450;
	    }
/* L350: */
	}
	if (kcmnd >= *maxcwd) {
	    goto L400;
	}
	left = *maxcwd - kcmnd;
	ltoadd = lelmnt[ielmnt - 1];
	if (ltoadd > left) {
	    ltoadd = left;
	}
	i__2 = kcmnd;
	s_copy(comand + i__2, celmnt + (ielmnt - 1) * 19, kcmnd + ltoadd - 
		i__2, ltoadd);
	kcmnd += ltoadd;
	if (kcmnd == *maxcwd) {
	    goto L400;
	}
	++kcmnd;
	*(unsigned char *)&comand[kcmnd - 1] = ' ';
L400:
	;
    }
    *lnc = kcmnd;
    goto L900;
L450:
    *lnc = kcmnd;
/*                      . . . .  we have come to a numeric field */
    *llist = 0;
    i__1 = nelmnt;
    for (ifld = ielmnt; ifld <= i__1; ++ifld) {
	++(*llist);
	if (*llist > *mxp) {
	    nreq = nelmnt - ielmnt + 1;
	    io___182.ciunit = *isyswr;
	    s_wsfe(&io___182);
	    do_fio(&c__1, (char *)&nreq, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*mxp), (ftnlen)sizeof(integer));
	    e_wsfe();
	    goto L900;
	}
	if (s_cmp(celmnt + (ifld - 1) * 19, cnull, (ftnlen)19, (ftnlen)15) == 
		0) {
	    plist[*llist] = (float)0.;
	} else {
	    ici__1.icierr = 1;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 19;
	    ici__1.iciunit = celmnt + (ifld - 1) * 19;
	    ici__1.icifmt = "(BN,F19.0)";
	    i__2 = s_rsfi(&ici__1);
	    if (i__2 != 0) {
		goto L575;
	    }
	    i__2 = do_fio(&c__1, (char *)&plist[*llist], (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L575;
	    }
	    i__2 = e_rsfi();
	    if (i__2 != 0) {
		goto L575;
	    }
	}
	goto L600;
L575:
	io___183.ciunit = *isyswr;
	s_wsfe(&io___183);
	do_fio(&c__1, " FORMAT ERROR IN NUMERIC FIELD: \"", (ftnlen)33);
	do_fio(&c__1, celmnt + (ifld - 1) * 19, lelmnt[ifld - 1]);
	do_fio(&c__1, "\"", (ftnlen)1);
	e_wsfe();
	*ierr = 1;
	plist[*llist] = (float)0.;
L600:
	;
    }
/*                                  end loop over numeric fields */
L900:
    if (*lnc <= 0) {
	*lnc = 1;
    }
    return 0;
} /* mncrck_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mncros_(U_fp fcn, doublereal *aopt, integer *iercr, U_fp 
	futil)
{
    /* Initialized data */

    static char charal[28+1] = " .ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double sqrt(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, x1, x2;
    static integer ik, it;
    static doublereal rt, aim, tla, tlf;
    static integer kex, ipt;
    static doublereal dfda, alsb[3], flsb[3], bmin, bmax;
    static integer inew;
    static doublereal zmid, sdev, zdir, zlim;
    static integer iout;
    static doublereal coeff[3], aleft, ecart;
    static integer ileft;
    static doublereal aulim, fdist;
    static integer ierev;
    static doublereal adist;
    static integer maxlk, ibest;
    static doublereal anext, fnext, slope;
    static logical ldebug;
    static doublereal ecarmn;
    static char chsign[4];
    extern /* Subroutine */ int mneval_(U_fp, doublereal *, doublereal *, 
	    integer *, U_fp);
    static doublereal ecarmx, determ, aminsv;
    extern /* Subroutine */ int mnwarn_(char *, char *, char *, ftnlen, 
	    ftnlen, ftnlen);
    static integer noless, iworst;
    extern /* Subroutine */ int mnpfit_(doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *);
    static integer iright;
    static doublereal smalla, aright;
    static integer itoohi;
    extern /* Subroutine */ int mnplot_(doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___200 = { 0, 0, 0, "(A,I8,A,F10.5,A,2F10.5)", 0 };
    static cilist io___203 = { 0, 0, 0, "(A,I8,A,F10.5,A,2F10.5)", 0 };
    static cilist io___207 = { 0, 0, 0, "(A,I8,A,F10.5,A,2F10.5)", 0 };
    static cilist io___212 = { 0, 0, 0, "(A,I8,A,F10.5,A,2F10.5)", 0 };
    static cilist io___229 = { 0, 0, 0, "(A)", 0 };
    static cilist io___237 = { 0, 0, 0, "(A,I8,A,F10.5,A,2F10.5)", 0 };
    static cilist io___240 = { 0, 0, 0, "(2X,A,A,I3)", 0 };
    static cilist io___241 = { 0, 0, 0, "(10X,A)", 0 };
    static cilist io___242 = { 0, 0, 0, "(10X,A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       Find point where MNEVAL=AMIN+UP, along the line through */
/* C       XMIDCR,YMIDCR with direction XDIRCR,YDIRCR,   where X and Y */
/* C       are parameters KE1CR and KE2CR.  If KE2CR=0 (from MINOS), */
/* C       only KE1CR is varied.  From MNCONT, both are varied. */
/* C       Crossing point is at */
/* C        (U(KE1),U(KE2)) = (XMID,YMID) + AOPT*(XDIR,YDIR) */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    ldebug = mn7flg_1.idbg[6] >= 1;
    aminsv = mn7min_1.amin;
/*        convergence when F is within TLF of AIM and next prediction */
/*        of AOPT is within TLA of previous value of AOPT */
    aim = mn7min_1.amin + mn7min_1.up;
    tlf = mn7min_1.up * .01;
    tla = .01;
    mn7rpt_1.xpt[0] = (float)0.;
    mn7rpt_1.ypt[0] = aim;
    *(unsigned char *)&mn7cpt_1.chpt[0] = ' ';
    ipt = 1;
    if (mn7xcr_1.ke2cr == 0) {
	mn7rpt_1.xpt[1] = (float)-1.;
	mn7rpt_1.ypt[1] = mn7min_1.amin;
	*(unsigned char *)&mn7cpt_1.chpt[1] = '.';
	ipt = 2;
    }
/*                    find the largest allowed A */
    aulim = (float)100.;
    for (ik = 1; ik <= 2; ++ik) {
	if (ik == 1) {
	    kex = mn7xcr_1.ke1cr;
	    zmid = mn7xcr_1.xmidcr;
	    zdir = mn7xcr_1.xdircr;
	} else {
	    if (mn7xcr_1.ke2cr == 0) {
		goto L100;
	    }
	    kex = mn7xcr_1.ke2cr;
	    zmid = mn7xcr_1.ymidcr;
	    zdir = mn7xcr_1.ydircr;
	}
	if (mn7inx_1.nvarl[kex - 1] <= 1) {
	    goto L100;
	}
	if (zdir == 0.) {
	    goto L100;
	}
	zlim = mn7ext_1.alim[kex - 1];
	if (zdir > 0.) {
	    zlim = mn7ext_1.blim[kex - 1];
	}
/* Computing MIN */
	d__1 = aulim, d__2 = (zlim - zmid) / zdir;
	aulim = min(d__1,d__2);
L100:
	;
    }
/*                  LSB = Line Search Buffer */
/*          first point */
    anext = (float)0.;
    *aopt = anext;
    mn7log_1.limset = FALSE_;
    if (aulim < *aopt + tla) {
	mn7log_1.limset = TRUE_;
    }
    mneval_((U_fp)fcn, &anext, &fnext, &ierev, (U_fp)futil);
/* debug printout: */
    if (ldebug) {
	io___200.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___200);
	do_fio(&c__1, " MNCROS: calls=", (ftnlen)15);
	do_fio(&c__1, (char *)&mn7cnv_1.nfcn, (ftnlen)sizeof(integer));
	do_fio(&c__1, "   AIM=", (ftnlen)7);
	do_fio(&c__1, (char *)&aim, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "  F,A=", (ftnlen)6);
	do_fio(&c__1, (char *)&fnext, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*aopt), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (ierev > 0) {
	goto L900;
    }
    if (mn7log_1.limset && fnext <= aim) {
	goto L930;
    }
    ++ipt;
    mn7rpt_1.xpt[ipt - 1] = anext;
    mn7rpt_1.ypt[ipt - 1] = fnext;
    *(unsigned char *)&mn7cpt_1.chpt[ipt - 1] = *(unsigned char *)&charal[ipt 
	    - 1];
    alsb[0] = anext;
    flsb[0] = fnext;
/* Computing MAX */
    d__1 = fnext, d__2 = aminsv + mn7min_1.up * (float).1;
    fnext = max(d__1,d__2);
    *aopt = sqrt(mn7min_1.up / (fnext - aminsv)) - (float)1.;
    if ((d__1 = fnext - aim, abs(d__1)) < tlf) {
	goto L800;
    }

    if (*aopt < -.5) {
	*aopt = -.5;
    }
    if (*aopt > 1.) {
	*aopt = 1.;
    }
    mn7log_1.limset = FALSE_;
    if (*aopt > aulim) {
	*aopt = aulim;
	mn7log_1.limset = TRUE_;
    }
    mneval_((U_fp)fcn, aopt, &fnext, &ierev, (U_fp)futil);
/* debug printout: */
    if (ldebug) {
	io___203.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___203);
	do_fio(&c__1, " MNCROS: calls=", (ftnlen)15);
	do_fio(&c__1, (char *)&mn7cnv_1.nfcn, (ftnlen)sizeof(integer));
	do_fio(&c__1, "   AIM=", (ftnlen)7);
	do_fio(&c__1, (char *)&aim, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "  F,A=", (ftnlen)6);
	do_fio(&c__1, (char *)&fnext, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*aopt), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (ierev > 0) {
	goto L900;
    }
    if (mn7log_1.limset && fnext <= aim) {
	goto L930;
    }
    alsb[1] = *aopt;
    ++ipt;
    mn7rpt_1.xpt[ipt - 1] = alsb[1];
    mn7rpt_1.ypt[ipt - 1] = fnext;
    *(unsigned char *)&mn7cpt_1.chpt[ipt - 1] = *(unsigned char *)&charal[ipt 
	    - 1];
    flsb[1] = fnext;
    dfda = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);
/*                   DFDA must be positive on the contour */
    if (dfda > 0.) {
	goto L460;
    }
L300:
    mnwarn_("D", "MNCROS    ", "Looking for slope of the right sign", (ftnlen)
	    1, (ftnlen)10, (ftnlen)35);
    maxlk = 15 - ipt;
    i__1 = maxlk;
    for (it = 1; it <= i__1; ++it) {
	alsb[0] = alsb[1];
	flsb[0] = flsb[1];
	*aopt = alsb[0] + (real) it * (float).2;
	mn7log_1.limset = FALSE_;
	if (*aopt > aulim) {
	    *aopt = aulim;
	    mn7log_1.limset = TRUE_;
	}
	mneval_((U_fp)fcn, aopt, &fnext, &ierev, (U_fp)futil);
/* debug printout: */
	if (ldebug) {
	    io___207.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___207);
	    do_fio(&c__1, " MNCROS: calls=", (ftnlen)15);
	    do_fio(&c__1, (char *)&mn7cnv_1.nfcn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, "   AIM=", (ftnlen)7);
	    do_fio(&c__1, (char *)&aim, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, "  F,A=", (ftnlen)6);
	    do_fio(&c__1, (char *)&fnext, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*aopt), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	if (ierev > 0) {
	    goto L900;
	}
	if (mn7log_1.limset && fnext <= aim) {
	    goto L930;
	}
	alsb[1] = *aopt;
	++ipt;
	mn7rpt_1.xpt[ipt - 1] = alsb[1];
	mn7rpt_1.ypt[ipt - 1] = fnext;
	*(unsigned char *)&mn7cpt_1.chpt[ipt - 1] = *(unsigned char *)&charal[
		ipt - 1];
	flsb[1] = fnext;
	dfda = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);
	if (dfda > 0.) {
	    goto L450;
	}
/* L400: */
    }
    mnwarn_("W", "MNCROS    ", "Cannot find slope of the right sign", (ftnlen)
	    1, (ftnlen)10, (ftnlen)35);
    goto L950;
L450:
/*                    we have two points with the right slope */
L460:
    *aopt = alsb[1] + (aim - flsb[1]) / dfda;
/* Computing MIN */
    d__3 = (d__1 = aim - flsb[0], abs(d__1)), d__4 = (d__2 = aim - flsb[1], 
	    abs(d__2));
    fdist = min(d__3,d__4);
/* Computing MIN */
    d__3 = (d__1 = *aopt - alsb[0], abs(d__1)), d__4 = (d__2 = *aopt - alsb[1]
	    , abs(d__2));
    adist = min(d__3,d__4);
    tla = .01;
    if (abs(*aopt) > 1.) {
	tla = abs(*aopt) * .01;
    }
    if (adist < tla && fdist < tlf) {
	goto L800;
    }
    if (ipt >= 15) {
	goto L950;
    }
    bmin = min(alsb[0],alsb[1]) - (float)1.;
    if (*aopt < bmin) {
	*aopt = bmin;
    }
    bmax = max(alsb[0],alsb[1]) + (float)1.;
    if (*aopt > bmax) {
	*aopt = bmax;
    }
/*                    Try a third point */
    mn7log_1.limset = FALSE_;
    if (*aopt > aulim) {
	*aopt = aulim;
	mn7log_1.limset = TRUE_;
    }
    mneval_((U_fp)fcn, aopt, &fnext, &ierev, (U_fp)futil);
/* debug printout: */
    if (ldebug) {
	io___212.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___212);
	do_fio(&c__1, " MNCROS: calls=", (ftnlen)15);
	do_fio(&c__1, (char *)&mn7cnv_1.nfcn, (ftnlen)sizeof(integer));
	do_fio(&c__1, "   AIM=", (ftnlen)7);
	do_fio(&c__1, (char *)&aim, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "  F,A=", (ftnlen)6);
	do_fio(&c__1, (char *)&fnext, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*aopt), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (ierev > 0) {
	goto L900;
    }
    if (mn7log_1.limset && fnext <= aim) {
	goto L930;
    }
    alsb[2] = *aopt;
    ++ipt;
    mn7rpt_1.xpt[ipt - 1] = alsb[2];
    mn7rpt_1.ypt[ipt - 1] = fnext;
    *(unsigned char *)&mn7cpt_1.chpt[ipt - 1] = *(unsigned char *)&charal[ipt 
	    - 1];
    flsb[2] = fnext;
    inew = 3;
/*                now we have three points, ask how many <AIM */
    ecarmn = (d__1 = fnext - aim, abs(d__1));
    ibest = 3;
    ecarmx = (float)0.;
    noless = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	ecart = (d__1 = flsb[i__ - 1] - aim, abs(d__1));
	if (ecart > ecarmx) {
	    ecarmx = ecart;
	    iworst = i__;
	}
	if (ecart < ecarmn) {
	    ecarmn = ecart;
	    ibest = i__;
	}
	if (flsb[i__ - 1] < aim) {
	    ++noless;
	}
/* L480: */
    }
    inew = ibest;
/*           if at least one on each side of AIM, fit a parabola */
    if (noless == 1 || noless == 2) {
	goto L500;
    }
/*           if all three are above AIM, third must be closest to AIM */
    if (noless == 0 && ibest != 3) {
	goto L950;
    }
/*           if all three below, and third is not best, then slope */
/*             has again gone negative, look for positive slope. */
    if (noless == 3 && ibest != 3) {
	alsb[1] = alsb[2];
	flsb[1] = flsb[2];
	goto L300;
    }
/*           in other cases, new straight line thru last two points */
    alsb[iworst - 1] = alsb[2];
    flsb[iworst - 1] = flsb[2];
    dfda = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);
    goto L460;
/*                parabola fit */
L500:
    mnpfit_(alsb, flsb, &c__3, coeff, &sdev);
    if (coeff[2] <= 0.) {
	mnwarn_("D", "MNCROS    ", "Curvature is negative near contour line.",
		 (ftnlen)1, (ftnlen)10, (ftnlen)40);
    }
/* Computing 2nd power */
    d__1 = coeff[1];
    determ = d__1 * d__1 - coeff[2] * (float)4. * (coeff[0] - aim);
    if (determ <= 0.) {
	mnwarn_("D", "MNCROS    ", "Problem 2, impossible determinant", (
		ftnlen)1, (ftnlen)10, (ftnlen)33);
	goto L950;
    }
/*                Find which root is the right one */
    rt = sqrt(determ);
    x1 = (-coeff[1] + rt) / (coeff[2] * (float)2.);
    x2 = (-coeff[1] - rt) / (coeff[2] * (float)2.);
    s1 = coeff[1] + x1 * (float)2. * coeff[2];
    s2 = coeff[1] + x2 * (float)2. * coeff[2];
    if (s1 * s2 > 0.) {
	io___229.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___229);
	do_fio(&c__1, " MNCONTour problem 1", (ftnlen)20);
	e_wsfe();
    }
    *aopt = x1;
    slope = s1;
    if (s2 > 0.) {
	*aopt = x2;
	slope = s2;
    }
/*         ask if converged */
    tla = .01;
    if (abs(*aopt) > 1.) {
	tla = abs(*aopt) * .01;
    }
    if ((d__1 = *aopt - alsb[ibest - 1], abs(d__1)) < tla && (d__2 = flsb[
	    ibest - 1] - aim, abs(d__2)) < tlf) {
	goto L800;
    }
    if (ipt >= 15) {
	goto L950;
    }
/*         see if proposed point is in acceptable zone between L and R */
/*         first find ILEFT, IRIGHT, IOUT and IBEST */
    ileft = 0;
    iright = 0;
    ibest = 1;
    ecarmx = (float)0.;
    ecarmn = (d__1 = aim - flsb[0], abs(d__1));
    for (i__ = 1; i__ <= 3; ++i__) {
	ecart = (d__1 = flsb[i__ - 1] - aim, abs(d__1));
	if (ecart < ecarmn) {
	    ecarmn = ecart;
	    ibest = i__;
	}
	if (ecart > ecarmx) {
	    ecarmx = ecart;
	}
	if (flsb[i__ - 1] > aim) {
	    if (iright == 0) {
		iright = i__;
	    } else if (flsb[i__ - 1] > flsb[iright - 1]) {
		iout = i__;
	    } else {
		iout = iright;
		iright = i__;
	    }
	} else if (ileft == 0) {
	    ileft = i__;
	} else if (flsb[i__ - 1] < flsb[ileft - 1]) {
	    iout = i__;
	} else {
	    iout = ileft;
	    ileft = i__;
	}
/* L550: */
    }
/*       avoid keeping a very bad point next time around */
    if (ecarmx > (d__1 = flsb[iout - 1] - aim, abs(d__1)) * (float)10.) {
	*aopt = *aopt * .5 + (alsb[iright - 1] + alsb[ileft - 1]) * .25;
    }
/*         knowing ILEFT and IRIGHT, get acceptable window */
    smalla = tla * (float).1;
    if (slope * smalla > tlf) {
	smalla = tlf / slope;
    }
    aleft = alsb[ileft - 1] + smalla;
    aright = alsb[iright - 1] - smalla;
/*         move proposed point AOPT into window if necessary */
    if (*aopt < aleft) {
	*aopt = aleft;
    }
    if (*aopt > aright) {
	*aopt = aright;
    }
    if (aleft > aright) {
	*aopt = (aleft + aright) * .5;
    }
/*         see if proposed point outside limits (should be impossible!) */
    mn7log_1.limset = FALSE_;
    if (*aopt > aulim) {
	*aopt = aulim;
	mn7log_1.limset = TRUE_;
    }
/*                  Evaluate function at new point AOPT */
    mneval_((U_fp)fcn, aopt, &fnext, &ierev, (U_fp)futil);
/* debug printout: */
    if (ldebug) {
	io___237.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___237);
	do_fio(&c__1, " MNCROS: calls=", (ftnlen)15);
	do_fio(&c__1, (char *)&mn7cnv_1.nfcn, (ftnlen)sizeof(integer));
	do_fio(&c__1, "   AIM=", (ftnlen)7);
	do_fio(&c__1, (char *)&aim, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "  F,A=", (ftnlen)6);
	do_fio(&c__1, (char *)&fnext, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*aopt), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (ierev > 0) {
	goto L900;
    }
    if (mn7log_1.limset && fnext <= aim) {
	goto L930;
    }
    ++ipt;
    mn7rpt_1.xpt[ipt - 1] = *aopt;
    mn7rpt_1.ypt[ipt - 1] = fnext;
    *(unsigned char *)&mn7cpt_1.chpt[ipt - 1] = *(unsigned char *)&charal[ipt 
	    - 1];
/*                Replace odd point by new one */
    alsb[iout - 1] = *aopt;
    flsb[iout - 1] = fnext;
/*          the new point may not be the best, but it is the only one */
/*          which could be good enough to pass convergence criteria */
    ibest = iout;
    goto L500;

/*       Contour has been located, return point to MNCONT OR MINOS */
L800:
    *iercr = 0;
    goto L1000;
/*                error in the minimization */
L900:
    if (ierev == 1) {
	goto L940;
    }
    goto L950;
/*                parameter up against limit */
L930:
    *iercr = 1;
    goto L1000;
/*                too many calls to FCN */
L940:
    *iercr = 2;
    goto L1000;
/*                cannot find next point */
L950:
    *iercr = 3;
/*                in any case */
L1000:
    if (ldebug) {
	itoohi = 0;
	i__1 = ipt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (mn7rpt_1.ypt[i__ - 1] > aim + mn7min_1.up) {
		mn7rpt_1.ypt[i__ - 1] = aim + mn7min_1.up;
		*(unsigned char *)&mn7cpt_1.chpt[i__ - 1] = '+';
		itoohi = 1;
	    }
/* L1100: */
	}
	s_copy(chsign, "POSI", (ftnlen)4, (ftnlen)4);
	if (mn7xcr_1.xdircr < 0.) {
	    s_copy(chsign, "NEGA", (ftnlen)4, (ftnlen)4);
	}
	if (mn7xcr_1.ke2cr == 0) {
	    io___240.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___240);
	    do_fio(&c__1, chsign, (ftnlen)4);
	    do_fio(&c__1, "TIVE MINOS ERROR, PARAMETER ", (ftnlen)28);
	    do_fio(&c__1, (char *)&mn7xcr_1.ke1cr, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	if (itoohi == 1) {
	    io___241.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___241);
	    do_fio(&c__1, "POINTS LABELLED \"+\" WERE TOO HIGH TO PLOT.", (
		    ftnlen)42);
	    e_wsfe();
	}
	if (*iercr == 1) {
	    io___242.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___242);
	    do_fio(&c__1, "RIGHTMOST POINT IS UP AGAINST LIMIT.", (ftnlen)36);
	    e_wsfe();
	}
	mnplot_(mn7rpt_1.xpt, mn7rpt_1.ypt, mn7cpt_1.chpt, &ipt, &
		mn7iou_1.isyswr, &mn7iou_1.npagwd, &mn7iou_1.npagln, (ftnlen)
		1);
    }
    return 0;
} /* mncros_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mncuve_(S_fp fcn, U_fp futil)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer i__, j;
    static doublereal dxdi;
    static integer ndex, iext;
    static doublereal wint;
    extern /* Subroutine */ int mndxdi_(doublereal *, integer *, doublereal *)
	    , mnmigr_(S_fp, U_fp), mnhess_(S_fp, U_fp), mnwarn_(char *, char *
	    , char *, ftnlen, ftnlen, ftnlen), mnwerr_();

    /* Fortran I/O blocks */
    static cilist io___243 = { 0, 0, 0, "(/A,A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Makes sure that the current point is a local */
/* C        minimum and that the error matrix exists, */
/* C        or at least something good enough for MINOS and MNCONT */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    if (mn7flg_1.isw[3] < 1) {
	io___243.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___243);
	do_fio(&c__1, " FUNCTION MUST BE MINIMIZED BEFORE CALLING ", (ftnlen)
		43);
	do_fio(&c__1, mn7tit_1.cfrom, (ftnlen)8);
	e_wsfe();
	mn7min_1.apsi = mn7min_1.epsi;
	mnmigr_((S_fp)fcn, (U_fp)futil);
    }
    if (mn7flg_1.isw[1] < 3) {
	mnhess_((S_fp)fcn, (U_fp)futil);
	if (mn7flg_1.isw[1] < 1) {
	    mnwarn_("W", mn7tit_1.cfrom, "NO ERROR MATRIX.  WILL IMPROVISE.", 
		    (ftnlen)1, (ftnlen)8, (ftnlen)33);
	    i__1 = mn7npr_1.npar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ndex = i__ * (i__ - 1) / 2;
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
		    ++ndex;
/* L554: */
		    mn7var_1.vhmat[ndex - 1] = (float)0.;
		}
		++ndex;
		if (mn7der_1.g2[i__ - 1] <= 0.) {
		    wint = mn7err_1.werr[i__ - 1];
		    iext = mn7inx_1.nexofi[i__ - 1];
		    if (mn7inx_1.nvarl[iext - 1] > 1) {
			mndxdi_(&mn7int_1.x[i__ - 1], &i__, &dxdi);
			if (abs(dxdi) < (float).001) {
			    wint = (float).01;
			} else {
			    wint /= abs(dxdi);
			}
		    }
/* Computing 2nd power */
		    d__1 = wint;
		    mn7der_1.g2[i__ - 1] = mn7min_1.up / (d__1 * d__1);
		}
		mn7var_1.vhmat[ndex - 1] = (float)2. / mn7der_1.g2[i__ - 1];
/* L555: */
	    }
	    mn7flg_1.isw[1] = 1;
	    mn7min_1.dcovar = (float)1.;
	} else {
	    mnwerr_();
	}
    }
    return 0;
} /* mncuve_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:43  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mnderi_(S_fp fcn, U_fp futil)
{
    /* Format strings */
    static char fmt_41[] = "(i4,2g11.3,5g10.2)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    doublereal d__1, d__2, d__3;
    char ch__1[48], ch__2[54];

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_wsfe(cilist *), e_wsfe();
    double sqrt(doublereal), d_sign(doublereal *, doublereal *), cos(
	    doublereal);

    /* Local variables */
    static integer i__;
    static doublereal dd, df, fs1, fs2, d1d2, xtf;
    static char cbf1[22];
    static integer icyc, ncyc, iint, iext;
    static doublereal step, dfmin;
    static integer nparx;
    static doublereal stepb4;
    static logical ldebug;
    extern /* Subroutine */ int mnamin_(S_fp, U_fp);
    static doublereal grbfor;
    extern /* Subroutine */ int mninex_(doublereal *), mnwarn_(char *, char *,
	     char *, ftnlen, ftnlen, ftnlen);
    static doublereal tlrgrd, epspri, tlrstp, optstp, stpmax, vrysml, stpmin;

    /* Fortran I/O blocks */
    static icilist io___255 = { 0, cbf1, 0, "(G12.3)", 12, 1 };
    static cilist io___256 = { 0, 0, 0, "(/'  FIRST DERIVATIVE DEBUG PRINTOU\
T.  MNDERI'/        ' PAR    DERIV     STEP      MINSTEP   OPTSTEP ',       \
        ' D1-D2    2ND DRV')", 0 };
    static cilist io___274 = { 0, 0, 0, fmt_41, 0 };
    static icilist io___275 = { 0, cbf1, 0, "(2E11.3)", 22, 1 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Calculates the first derivatives of FCN (GRD), */
/* C        either by finite differences or by transforming the user- */
/* C        supplied derivatives to internal coordinates, */
/* C        according to whether ISW(3) is zero or one. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    nparx = mn7npr_1.npar;
    ldebug = mn7flg_1.idbg[2] >= 1;
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    if (mn7flg_1.isw[2] == 1) {
	goto L100;
    }
    if (ldebug) {
/*                       make sure starting at the right place */
	mninex_(mn7int_1.x);
	nparx = mn7npr_1.npar;
	(*fcn)(&nparx, mn7der_1.gin, &fs1, mn7ext_1.u, &c__4, (U_fp)futil);
	++mn7cnv_1.nfcn;
	if (fs1 != mn7min_1.amin) {
	    df = mn7min_1.amin - fs1;
	    s_wsfi(&io___255);
	    do_fio(&c__1, (char *)&df, (ftnlen)sizeof(doublereal));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 36, a__1[0] = "function value differs from AMIN by ";
	    i__1[1] = 12, a__1[1] = cbf1;
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)48);
	    mnwarn_("D", "MNDERI", ch__1, (ftnlen)1, (ftnlen)6, (ftnlen)48);
	    mn7min_1.amin = fs1;
	}
	io___256.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___256);
	e_wsfe();
    }
    dfmin = mn7cns_1.epsma2 * (float)8. * (abs(mn7min_1.amin) + mn7min_1.up);
/* Computing 2nd power */
    d__1 = mn7cns_1.epsmac;
    vrysml = d__1 * d__1 * (float)8.;
    if (mn7cnv_1.istrat <= 0) {
	ncyc = 2;
	tlrstp = (float).5;
	tlrgrd = (float).1;
    } else if (mn7cnv_1.istrat == 1) {
	ncyc = 3;
	tlrstp = (float).3;
	tlrgrd = (float).05;
    } else {
	ncyc = 5;
	tlrstp = (float).1;
	tlrgrd = (float).02;
    }
/*                                loop over variable parameters */
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	epspri = mn7cns_1.epsma2 + (d__1 = mn7der_1.grd[i__ - 1] * 
		mn7cns_1.epsma2, abs(d__1));
/*         two-point derivatives always assumed necessary */
/*         maximum number of cycles over step size depends on strategy */
	xtf = mn7int_1.x[i__ - 1];
	stepb4 = (float)0.;
/*                               loop as little as possible here! */
	i__3 = ncyc;
	for (icyc = 1; icyc <= i__3; ++icyc) {
/*                 ........ theoretically best step */
	    optstp = sqrt(dfmin / ((d__1 = mn7der_1.g2[i__ - 1], abs(d__1)) + 
		    epspri));
/*                     step cannot decrease by more than a factor of ten */
/* Computing MAX */
	    d__2 = optstp, d__3 = (d__1 = mn7der_1.gstep[i__ - 1] * (float).1,
		     abs(d__1));
	    step = max(d__2,d__3);
/*                 but if parameter has limits, max step size = 0.5 */
	    if (mn7der_1.gstep[i__ - 1] < 0. && step > (float).5) {
		step = (float).5;
	    }
/*                 and not more than ten times the previous step */
	    stpmax = (d__1 = mn7der_1.gstep[i__ - 1], abs(d__1)) * (float)10.;
	    if (step > stpmax) {
		step = stpmax;
	    }
/*                 minimum step size allowed by machine precision */
/* Computing MAX */
	    d__2 = vrysml, d__3 = (d__1 = mn7cns_1.epsma2 * mn7int_1.x[i__ - 
		    1], abs(d__1)) * (float)8.;
	    stpmin = max(d__2,d__3);
	    if (step < stpmin) {
		step = stpmin;
	    }
/*                 end of iterations if step change less than factor 2 */
	    if ((d__1 = (step - stepb4) / step, abs(d__1)) < tlrstp) {
		goto L50;
	    }
/*         take step positive */
	    mn7der_1.gstep[i__ - 1] = d_sign(&step, &mn7der_1.gstep[i__ - 1]);
	    stepb4 = step;
	    mn7int_1.x[i__ - 1] = xtf + step;
	    mninex_(mn7int_1.x);
	    (*fcn)(&nparx, mn7der_1.gin, &fs1, mn7ext_1.u, &c__4, (U_fp)futil)
		    ;
	    ++mn7cnv_1.nfcn;
/*         take step negative */
	    mn7int_1.x[i__ - 1] = xtf - step;
	    mninex_(mn7int_1.x);
	    (*fcn)(&nparx, mn7der_1.gin, &fs2, mn7ext_1.u, &c__4, (U_fp)futil)
		    ;
	    ++mn7cnv_1.nfcn;
	    grbfor = mn7der_1.grd[i__ - 1];
	    mn7der_1.grd[i__ - 1] = (fs1 - fs2) / (step * (float)2.);
/* Computing 2nd power */
	    d__1 = step;
	    mn7der_1.g2[i__ - 1] = (fs1 + fs2 - mn7min_1.amin * (float)2.) / (
		    d__1 * d__1);
	    mn7int_1.x[i__ - 1] = xtf;
	    if (ldebug) {
		d1d2 = (fs1 + fs2 - mn7min_1.amin * (float)2.) / step;
		io___274.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___274);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&mn7der_1.grd[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&step, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&stpmin, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&optstp, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&d1d2, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&mn7der_1.g2[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
/*         see if another iteration is necessary */
	    if ((d__2 = grbfor - mn7der_1.grd[i__ - 1], abs(d__2)) / ((d__1 = 
		    mn7der_1.grd[i__ - 1], abs(d__1)) + dfmin / step) < 
		    tlrgrd) {
		goto L50;
	    }
/* L45: */
	}
/*                           end of ICYC loop. too many iterations */
	if (ncyc == 1) {
	    goto L50;
	}
	s_wsfi(&io___275);
	do_fio(&c__1, (char *)&mn7der_1.grd[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&grbfor, (ftnlen)sizeof(doublereal));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 32, a__1[0] = "First derivative not converged. ";
	i__1[1] = 22, a__1[1] = cbf1;
	s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)54);
	mnwarn_("D", "MNDERI", ch__2, (ftnlen)1, (ftnlen)6, (ftnlen)54);
L50:

/* L60: */
	;
    }
    mninex_(mn7int_1.x);
    return 0;
/*                                        .  derivatives calc by fcn */
L100:
    i__2 = mn7npr_1.npar;
    for (iint = 1; iint <= i__2; ++iint) {
	iext = mn7inx_1.nexofi[iint - 1];
	if (mn7inx_1.nvarl[iext - 1] > 1) {
	    goto L120;
	}
	mn7der_1.grd[iint - 1] = mn7der_1.gin[iext - 1];
	goto L150;
L120:
	dd = (mn7ext_1.blim[iext - 1] - mn7ext_1.alim[iext - 1]) * (float).5 *
		 cos(mn7int_1.x[iint - 1]);
	mn7der_1.grd[iint - 1] = mn7der_1.gin[iext - 1] * dd;
L150:
	;
    }
/* L200: */
    return 0;
} /* mnderi_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mndxdi_(doublereal *pint, integer *ipar, doublereal *
	dxdi)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer i__;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        calculates the transformation factor between external and */
/* C        internal parameter values.     this factor is one for */
/* C        parameters which are not limited.     called from MNEMAT. */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    i__ = mn7inx_1.nexofi[*ipar - 1];
    *dxdi = (float)1.;
    if (mn7inx_1.nvarl[i__ - 1] > 1) {
	*dxdi = (d__1 = (mn7ext_1.blim[i__ - 1] - mn7ext_1.alim[i__ - 1]) * 
		cos(*pint), abs(d__1)) * (float).5;
    }
    return 0;
} /* mndxdi_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mneig_(doublereal *a, integer *ndima, integer *n, 
	integer *mits, doublereal *work, doublereal *precis, integer *ifault)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal b, c__, f, h__;
    static integer i__, j, k, l, m;
    static doublereal r__, s;
    static integer i0, i1, j1, m1, n1;
    static doublereal hh, gl, pr, pt;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */

/*          PRECIS is the machine precision EPSMAC */
    /* Parameter adjustments */
    a_dim1 = *ndima;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --work;

    /* Function Body */
    *ifault = 1;

    i__ = *n;
    i__1 = *n;
    for (i1 = 2; i1 <= i__1; ++i1) {
	l = i__ - 2;
	f = a[i__ + (i__ - 1) * a_dim1];
	gl = 0.;

	if (l < 1) {
	    goto L25;
	}

	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
/* L20: */
/* Computing 2nd power */
	    d__1 = a[i__ + k * a_dim1];
	    gl += d__1 * d__1;
	}
L25:
/* Computing 2nd power */
	d__1 = f;
	h__ = gl + d__1 * d__1;

	if (gl > 1e-35) {
	    goto L30;
	}

	work[i__] = 0.;
	work[*n + i__] = f;
	goto L65;
L30:
	++l;

	gl = sqrt(h__);

	if (f >= 0.) {
	    gl = -gl;
	}

	work[*n + i__] = gl;
	h__ -= f * gl;
	a[i__ + (i__ - 1) * a_dim1] = f - gl;
	f = 0.;
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    a[j + i__ * a_dim1] = a[i__ + j * a_dim1] / h__;
	    gl = 0.;
	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
/* L40: */
		gl += a[j + k * a_dim1] * a[i__ + k * a_dim1];
	    }

	    if (j >= l) {
		goto L47;
	    }

	    j1 = j + 1;
	    i__3 = l;
	    for (k = j1; k <= i__3; ++k) {
/* L45: */
		gl += a[k + j * a_dim1] * a[i__ + k * a_dim1];
	    }
L47:
	    work[*n + j] = gl / h__;
	    f += gl * a[j + i__ * a_dim1];
/* L50: */
	}
	hh = f / (h__ + h__);
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = a[i__ + j * a_dim1];
	    gl = work[*n + j] - hh * f;
	    work[*n + j] = gl;
	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		a[j + k * a_dim1] = a[j + k * a_dim1] - f * work[*n + k] - gl 
			* a[i__ + k * a_dim1];
/* L60: */
	    }
	}
	work[i__] = h__;
L65:
	--i__;
/* L70: */
    }
    work[1] = 0.;
    work[*n + 1] = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = i__ - 1;

	if (work[i__] == 0. || l == 0) {
	    goto L100;
	}

	i__3 = l;
	for (j = 1; j <= i__3; ++j) {
	    gl = 0.;
	    i__2 = l;
	    for (k = 1; k <= i__2; ++k) {
/* L80: */
		gl += a[i__ + k * a_dim1] * a[k + j * a_dim1];
	    }
	    i__2 = l;
	    for (k = 1; k <= i__2; ++k) {
		a[k + j * a_dim1] -= gl * a[k + i__ * a_dim1];
/* L90: */
	    }
	}
L100:
	work[i__] = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = 1.;

	if (l == 0) {
	    goto L110;
	}

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    a[i__ + j * a_dim1] = 0.;
	    a[j + i__ * a_dim1] = 0.;
/* L105: */
	}
L110:
	;
    }


    n1 = *n - 1;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i0 = *n + i__ - 1;
/* L130: */
	work[i0] = work[i0 + 1];
    }
    work[*n + *n] = 0.;
    b = 0.;
    f = 0.;
    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
	h__ = *precis * ((d__1 = work[l], abs(d__1)) + (d__2 = work[*n + l], 
		abs(d__2)));

	if (b < h__) {
	    b = h__;
	}

	i__2 = *n;
	for (m1 = l; m1 <= i__2; ++m1) {
	    m = m1;

	    if ((d__1 = work[*n + m], abs(d__1)) <= b) {
		goto L150;
	    }

/* L140: */
	}

L150:
	if (m == l) {
	    goto L205;
	}

L160:
	if (j == *mits) {
	    return 0;
	}

	++j;
	pt = (work[l + 1] - work[l]) / (work[*n + l] * 2.);
	r__ = sqrt(pt * pt + 1.);
	pr = pt + r__;

	if (pt < 0.) {
	    pr = pt - r__;
	}

	h__ = work[l] - work[*n + l] / pr;
	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
/* L170: */
	    work[i__] -= h__;
	}
	f += h__;
	pt = work[m];
	c__ = 1.;
	s = 0.;
	m1 = m - 1;
	i__ = m;
	i__2 = m1;
	for (i1 = l; i1 <= i__2; ++i1) {
	    j = i__;
	    --i__;
	    gl = c__ * work[*n + i__];
	    h__ = c__ * pt;

	    if (abs(pt) >= (d__1 = work[*n + i__], abs(d__1))) {
		goto L180;
	    }

	    c__ = pt / work[*n + i__];
	    r__ = sqrt(c__ * c__ + 1.);
	    work[*n + j] = s * work[*n + i__] * r__;
	    s = 1. / r__;
	    c__ /= r__;
	    goto L190;
L180:
	    c__ = work[*n + i__] / pt;
	    r__ = sqrt(c__ * c__ + 1.);
	    work[*n + j] = s * pt * r__;
	    s = c__ / r__;
	    c__ = 1. / r__;
L190:
	    pt = c__ * work[i__] - s * gl;
	    work[j] = h__ + s * (c__ * gl + s * work[i__]);
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		h__ = a[k + j * a_dim1];
		a[k + j * a_dim1] = s * a[k + i__ * a_dim1] + c__ * h__;
		a[k + i__ * a_dim1] = c__ * a[k + i__ * a_dim1] - s * h__;
/* L200: */
	    }
	}
	work[*n + l] = s * pt;
	work[l] = c__ * pt;

	if ((d__1 = work[*n + l], abs(d__1)) > b) {
	    goto L160;
	}

L205:
	work[l] += f;
/* L210: */
    }
    i__1 = n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i__;
	pt = work[i__];
	i1 = i__ + 1;
	i__3 = *n;
	for (j = i1; j <= i__3; ++j) {

	    if (work[j] >= pt) {
		goto L220;
	    }

	    k = j;
	    pt = work[j];
L220:
	    ;
	}

	if (k == i__) {
	    goto L240;
	}

	work[k] = work[i__];
	work[i__] = pt;
	i__3 = *n;
	for (j = 1; j <= i__3; ++j) {
	    pt = a[j + i__ * a_dim1];
	    a[j + i__ * a_dim1] = a[j + k * a_dim1];
	    a[j + k * a_dim1] = pt;
/* L230: */
	}
L240:
	;
    }
    *ifault = 0;

    return 0;
} /* mneig_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mnemat_(doublereal *emat, integer *ndim)
{
    /* System generated locals */
    integer emat_dim1, emat_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer i__, j, k, k2, kk, iz, kga, kgb;
    static doublereal dxdi, dxdj;
    static integer npard;
    extern /* Subroutine */ int mndxdi_(doublereal *, integer *, doublereal *)
	    ;
    static integer nperln;

    /* Fortran I/O blocks */
    static cilist io___300 = { 0, 0, 0, "(/A,I4,A,I3,A,G10.3)", 0 };
    static cilist io___302 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___304 = { 0, 0, 0, "(A)", 0 };
    static cilist io___314 = { 0, 0, 0, "(1X,13E10.3)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Calculates the external error matrix from the internal */
/* C        to be called by user, who must dimension EMAT at (NDIM,NDIM) */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    /* Parameter adjustments */
    emat_dim1 = *ndim;
    emat_offset = 1 + emat_dim1 * 1;
    emat -= emat_offset;

    /* Function Body */
    if (mn7flg_1.isw[1] < 1) {
	return 0;
    }
    if (mn7flg_1.isw[4] >= 2) {
	io___300.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___300);
	do_fio(&c__1, " EXTERNAL ERROR MATRIX.    NDIM=", (ftnlen)32);
	do_fio(&c__1, (char *)&(*ndim), (ftnlen)sizeof(integer));
	do_fio(&c__1, "    NPAR=", (ftnlen)9);
	do_fio(&c__1, (char *)&mn7npr_1.npar, (ftnlen)sizeof(integer));
	do_fio(&c__1, "    ERR DEF=", (ftnlen)12);
	do_fio(&c__1, (char *)&mn7min_1.up, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*                    size of matrix to be printed */
    npard = mn7npr_1.npar;
    if (*ndim < mn7npr_1.npar) {
	npard = *ndim;
	if (mn7flg_1.isw[4] >= 0) {
	    io___302.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___302);
	    do_fio(&c__1, " USER-DIMENSIONED ", (ftnlen)18);
	    do_fio(&c__1, " ARRAY EMAT NOT BIG ENOUGH. REDUCED MATRIX CALCUL\
ATED.", (ftnlen)54);
	    e_wsfe();
	}
    }
/*                 NPERLN is the number of elements that fit on one line */
    nperln = (mn7iou_1.npagwd - 5) / 10;
    nperln = min(nperln,13);
    if (mn7flg_1.isw[4] >= 1 && npard > nperln) {
	io___304.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___304);
	do_fio(&c__1, " ELEMENTS ABOVE DIAGONAL ARE NOT PRINTED.", (ftnlen)41)
		;
	e_wsfe();
    }
/*                 I counts the rows of the matrix */
    i__1 = npard;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mndxdi_(&mn7int_1.x[i__ - 1], &i__, &dxdi);
	kga = i__ * (i__ - 1) / 2;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    mndxdi_(&mn7int_1.x[j - 1], &j, &dxdj);
	    kgb = kga + j;
	    emat[i__ + j * emat_dim1] = dxdi * mn7var_1.vhmat[kgb - 1] * dxdj 
		    * mn7min_1.up;
	    emat[j + i__ * emat_dim1] = emat[i__ + j * emat_dim1];
/* L100: */
	}
/* L110: */
    }
/*                    IZ is number of columns to be printed in row I */
    if (mn7flg_1.isw[4] >= 2) {
	i__1 = npard;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iz = npard;
	    if (npard >= nperln) {
		iz = i__;
	    }
	    i__2 = iz;
	    i__3 = nperln;
	    for (k = 1; i__3 < 0 ? k >= i__2 : k <= i__2; k += i__3) {
		k2 = k + nperln - 1;
		if (k2 > iz) {
		    k2 = iz;
		}
		io___314.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___314);
		i__4 = k2;
		for (kk = k; kk <= i__4; ++kk) {
		    do_fio(&c__1, (char *)&emat[i__ + kk * emat_dim1], (
			    ftnlen)sizeof(doublereal));
		}
		e_wsfe();
/* L150: */
	    }
/* L160: */
	}
    }
    return 0;
} /* mnemat_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mnerrs_(integer *number, doublereal *eplus, doublereal *
	eminus, doublereal *eparab, doublereal *gcc)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer iin, iex;
    static doublereal dxdi;
    static integer ndiag;
    extern /* Subroutine */ int mndxdi_(doublereal *, integer *, doublereal *)
	    ;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C    Called by user, utility routine to get MINOS errors */
/* C    If NUMBER is positive, then it is external parameter number, */
/* C                  if negative, it is -internal number. */
/* C    values returned by MNERRS: */
/* C       EPLUS, EMINUS are MINOS errors of parameter NUMBER, */
/* C       EPARAB is 'parabolic' error (from error matrix). */
/* C                 (Errors not calculated are set = 0.) */
/* C       GCC is global correlation coefficient from error matrix */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



    iex = *number;
    if (*number < 0) {
	iin = -(*number);
	if (iin > mn7npr_1.npar) {
	    goto L900;
	}
	iex = mn7inx_1.nexofi[iin - 1];
    }
    if (iex > mn7npr_1.nu || iex <= 0) {
	goto L900;
    }
    iin = mn7inx_1.niofex[iex - 1];
    if (iin <= 0) {
	goto L900;
    }
/*             IEX is external number, IIN is internal number */
    *eplus = mn7err_1.erp[iin - 1];
    if (*eplus == mn7cns_1.undefi) {
	*eplus = (float)0.;
    }
    *eminus = mn7err_1.ern[iin - 1];
    if (*eminus == mn7cns_1.undefi) {
	*eminus = (float)0.;
    }
    mndxdi_(&mn7int_1.x[iin - 1], &iin, &dxdi);
    ndiag = iin * (iin + 1) / 2;
    *eparab = (d__2 = dxdi * sqrt((d__1 = mn7min_1.up * mn7var_1.vhmat[ndiag 
	    - 1], abs(d__1))), abs(d__2));
/*              global correlation coefficient */
    *gcc = (float)0.;
    if (mn7flg_1.isw[1] < 2) {
	goto L990;
    }
    *gcc = mn7err_1.globcc[iin - 1];
    goto L990;
/*                  ERROR.  parameter number not valid */
L900:
    *eplus = (float)0.;
    *eminus = (float)0.;
    *eparab = (float)0.;
    *gcc = (float)0.;
L990:
    return 0;
} /* mnerrs_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mneval_(S_fp fcn, doublereal *anext, doublereal *fnext, 
	integer *ierev, U_fp futil)
{
    static integer nparx;
    extern /* Subroutine */ int mninex_(doublereal *), mnmigr_(S_fp, U_fp);


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C      Evaluates the function being analyzed by MNCROS, which is */
/* C      generally the minimum of FCN with respect to all remaining */
/* C      variable parameters.  Common block /MN7XCR/ contains the */
/* C      data necessary to know the values of U(KE1CR) and U(KE2CR) */
/* C      to be used, namely     U(KE1CR) = XMIDCR + ANEXT*XDIRCR */
/* C      and (if KE2CR .NE. 0)  U(KE2CR) = YMIDCR + ANEXT*YDIRCR */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/* C */
    mn7ext_1.u[mn7xcr_1.ke1cr - 1] = mn7xcr_1.xmidcr + *anext * 
	    mn7xcr_1.xdircr;
    if (mn7xcr_1.ke2cr != 0) {
	mn7ext_1.u[mn7xcr_1.ke2cr - 1] = mn7xcr_1.ymidcr + *anext * 
		mn7xcr_1.ydircr;
    }
    mninex_(mn7int_1.x);
    nparx = mn7npr_1.npar;
    (*fcn)(&nparx, mn7der_1.gin, fnext, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    *ierev = 0;
    if (mn7npr_1.npar > 0) {
	mn7cnv_1.itaur = 1;
	mn7min_1.amin = *fnext;
	mn7flg_1.isw[0] = 0;
	mnmigr_((S_fp)fcn, (U_fp)futil);
	mn7cnv_1.itaur = 0;
	*fnext = mn7min_1.amin;
	if (mn7flg_1.isw[0] >= 1) {
	    *ierev = 1;
	}
	if (mn7flg_1.isw[3] < 1) {
	    *ierev = 2;
	}
    }
    return 0;
} /* mneval_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:45  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mnexcm_(S_fp fcn, char *comand, doublereal *plist, 
	integer *llist, integer *ierflg, U_fp futil, ftnlen comand_len)
{
    /* Initialized data */

    static char clower[26+1] = "abcdefghijklmnopqrstuvwxyz";
    static char cupper[26+1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    static char cname[10*40+1] = "MINImize  SEEk      SIMplex   MIGrad    MI\
NOs     SET xxx   SHOw xxx  TOP of pagFIX       REStore   RELease   SCAn    \
  CONtour   HESse     SAVe      IMProve   CALl fcn  STAndard  END       EXIt\
      RETurn    CLEar     HELP      MNContour STOp      JUMp                \
                                                            COVARIANCEPRINTO\
UT  GRADIENT  MATOUT    ERROR DEF LIMITS    PUNCH     ";
    static integer nntot = 40;

    /* Format strings */
    static char fmt_25[] = "(\002 \002,10(\002*\002)/\002 **\002,i5,\002 *\
*\002,a,4g12.4)";
    static char fmt_3101[] = "(\002 OBSOLETE COMMAND:\002,1x,a10,5x,\002PLEA\
SE USE:\002,1x,a10)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2[3];

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfe(), s_wsfi(icilist *), e_wsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static doublereal f;
    static integer i__;
    static char c26[30];
    static integer nf, lk, it, iw, ke1, ke2, it2, kll, let, krl;
    static doublereal rno;
    static char comd[4];
    static integer icol, kcol, ierr, iint, iext;
    static doublereal step;
    static integer lnow, nptu;
    static doublereal xptu[101], yptu[101];
    static integer iflag, ierrf;
    extern /* Subroutine */ int stand_();
    static char chwhy[18];
    extern /* Subroutine */ int mnset_(S_fp, U_fp);
    static integer ilist, nparx, izero;
    extern /* Subroutine */ int mnrn15_(doublereal *, integer *);
    static logical lfreed;
    static char cvblnk[2];
    static logical lfixed;
    static integer inonde;
    static char cneway[10];
    static logical ltofix;
    extern /* Subroutine */ int mnseek_(S_fp, U_fp), mnsimp_(S_fp, U_fp), 
	    mnmigr_(S_fp, U_fp), mnwerr_();
    static integer nsuper;
    extern /* Subroutine */ int mncuve_(S_fp, U_fp), mnmnos_(S_fp, U_fp), 
	    mnrset_(integer *), mnfixp_(integer *, integer *), mnfree_(
	    integer *), mnprin_(integer *, doublereal *), mnscan_(S_fp, U_fp),
	     mncntr_(S_fp, integer *, integer *, integer *, U_fp), mnhess_(
	    S_fp, U_fp), mnmatu_(integer *), mnsave_(), mnimpr_(S_fp, U_fp);
    static integer nowprt;
    extern /* Subroutine */ int mncler_(), mnhelp_(char *, integer *, ftnlen),
	     mncont_(S_fp, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, U_fp), mninex_(doublereal *), mnamin_(
	    S_fp, U_fp);

    /* Fortran I/O blocks */
    static cilist io___330 = { 0, 0, 0, fmt_25, 0 };
    static icilist io___335 = { 0, cvblnk, 0, "(I2)", 2, 1 };
    static cilist io___337 = { 0, 0, 0, c26, 0 };
    static cilist io___338 = { 0, 0, 0, "(1H ,10(1H*))", 0 };
    static cilist io___339 = { 0, 0, 0, "(1H ,10(1H*),A,I3,A)", 0 };
    static cilist io___340 = { 0, 0, 0, "(11X,'UNKNOWN COMMAND IGNORED:',A)", 
	    0 };
    static cilist io___343 = { 0, 0, 0, "(/' TOO MANY FUNCTION CALLS. MINOS \
GIVES UP'/)", 0 };
    static cilist io___344 = { 0, 0, 0, "(1H1)", 0 };
    static cilist io___348 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___355 = { 0, 0, 0, "(A,I4,A,A)", 0 };
    static cilist io___357 = { 0, 0, 0, "(A,I4)", 0 };
    static cilist io___359 = { 0, 0, 0, "(A,I4,A)", 0 };
    static cilist io___362 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___368 = { 0, 0, 0, "(/A/)", 0 };
    static cilist io___369 = { 0, 0, 0, "(A)", 0 };
    static cilist io___378 = { 0, 0, 0, "(10X,A)", 0 };
    static cilist io___379 = { 0, 0, 0, "(A)", 0 };
    static cilist io___381 = { 0, 0, 0, fmt_3101, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Interprets a command and takes appropriate action, */
/* C        either directly by skipping to the corresponding code in */
/* C        MNEXCM, or by setting up a call to a subroutine */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/*   Cannot say DIMENSION PLIST(LLIST) since LLIST can be =0. */
/*  alphabetical order of command names! */

    /* Parameter adjustments */
    --plist;

    /* Function Body */

/*  recognized MINUIT commands: */
/*  obsolete commands: */
/*      IERFLG is now (94.5) defined the same as ICONDN in MNCOMD */
/* C            = 0: command executed normally */
/* C              1: command is blank, ignored */
/* C              2: command line unreadable, ignored */
/* C              3: unknown command, ignored */
/* C              4: abnormal termination (e.g., MIGRAD not converged) */
/* C              9: reserved */
/* C             10: END command */
/* C             11: EXIT or STOP command */
/* C             12: RETURN command */
    lk = i_len(comand, comand_len);
    if (lk > 20) {
	lk = 20;
    }
    s_copy(mn7tit_1.cword, comand, (ftnlen)20, lk);
/*              get upper case */
    i__1 = lk;
    for (icol = 1; icol <= i__1; ++icol) {
	for (let = 1; let <= 26; ++let) {
	    if (*(unsigned char *)&mn7tit_1.cword[icol - 1] == *(unsigned 
		    char *)&clower[let - 1]) {
		*(unsigned char *)&mn7tit_1.cword[icol - 1] = *(unsigned char 
			*)&cupper[let - 1];
	    }
/* L15: */
	}
/* L16: */
    }
/*           Copy the first MAXP arguments into COMMON (WORD7), making */
/*           sure that WORD7(1)=0. if LLIST=0 */
    for (iw = 1; iw <= 30; ++iw) {
	mn7arg_1.word7[iw - 1] = 0.;
	if (iw <= *llist) {
	    mn7arg_1.word7[iw - 1] = plist[iw];
	}
/* L20: */
    }
    ++mn7flg_1.icomnd;
    mn7cnv_1.nfcnlc = mn7cnv_1.nfcn;
    if (s_cmp(mn7tit_1.cword, "SET PRI", (ftnlen)7, (ftnlen)7) != 0 || 
	    mn7arg_1.word7[0] >= (float)0.) {
	if (mn7flg_1.isw[4] >= 0) {
	    lnow = *llist;
	    if (lnow > 4) {
		lnow = 4;
	    }
	    io___330.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___330);
	    do_fio(&c__1, (char *)&mn7flg_1.icomnd, (ftnlen)sizeof(integer));
	    do_fio(&c__1, mn7tit_1.cword, lk);
	    i__1 = lnow;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&plist[i__], (ftnlen)sizeof(doublereal))
			;
	    }
	    e_wsfe();
	    inonde = 0;
	    if (*llist > lnow) {
		kll = *llist;
		if (*llist > 30) {
		    inonde = 1;
		    kll = 30;
		}
		s_wsfi(&io___335);
		do_fio(&c__1, (char *)&lk, (ftnlen)sizeof(integer));
		e_wsfi();
/* Writing concatenation */
		i__2[0] = 16, a__1[0] = "(11H **********,";
		i__2[1] = 2, a__1[1] = cvblnk;
		i__2[2] = 9, a__1[2] = "X,4G12.4)";
		s_cat(c26, a__1, i__2, &c__3, (ftnlen)30);
		io___337.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___337);
		i__1 = kll;
		for (i__ = lnow + 1; i__ <= i__1; ++i__) {
		    do_fio(&c__1, (char *)&plist[i__], (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	    io___338.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___338);
	    e_wsfe();
	    if (inonde > 0) {
		io___339.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___339);
		do_fio(&c__1, "  ERROR: ABOVE CALL TO MNEXCM TRIED TO PASS M\
ORE THAN ", (ftnlen)54);
		do_fio(&c__1, (char *)&c__30, (ftnlen)sizeof(integer));
		do_fio(&c__1, " PARAMETERS.", (ftnlen)12);
		e_wsfe();
	    }
	}
    }
    mn7cnv_1.nfcnmx = (integer) mn7arg_1.word7[0];
    if (mn7cnv_1.nfcnmx <= 0) {
/* Computing 2nd power */
	i__1 = mn7npr_1.npar;
	mn7cnv_1.nfcnmx = mn7npr_1.npar * 100 + 200 + i__1 * i__1 * 5;
    }
    mn7min_1.epsi = mn7arg_1.word7[1];
    if (mn7min_1.epsi <= 0.) {
	mn7min_1.epsi = mn7min_1.up * (float).1;
    }
    mn7log_1.lnewmn = FALSE_;
    mn7log_1.lphead = TRUE_;
    mn7flg_1.isw[0] = 0;
    *ierflg = 0;
/*                look for command in list CNAME . . . . . . . . . . */
    i__1 = nntot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(mn7tit_1.cword, cname + (i__ - 1) * 10, (ftnlen)3, (ftnlen)
		3) == 0) {
	    goto L90;
	}
/* L80: */
    }
    io___340.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___340);
    do_fio(&c__1, comand, comand_len);
    e_wsfe();
    *ierflg = 3;
    goto L5000;
/*                normal case: recognized MINUIT command . . . . . . . */
L90:
    if (s_cmp(mn7tit_1.cword, "MINO", (ftnlen)4, (ftnlen)4) == 0) {
	i__ = 5;
    }
    if (i__ != 6 && i__ != 7 && i__ != 8 && i__ != 23) {
	s_copy(mn7tit_1.cfrom, cname + (i__ - 1) * 10, (ftnlen)8, (ftnlen)10);
	mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    }
/*              1    2    3    4    5    6    7    8    9   10 */
    switch (i__) {
	case 1:  goto L400;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L700;
	case 7:  goto L700;
	case 8:  goto L800;
	case 9:  goto L900;
	case 10:  goto L1000;
	case 11:  goto L1100;
	case 12:  goto L1200;
	case 13:  goto L1300;
	case 14:  goto L1400;
	case 15:  goto L1500;
	case 16:  goto L1600;
	case 17:  goto L1700;
	case 18:  goto L1800;
	case 19:  goto L1900;
	case 20:  goto L1900;
	case 21:  goto L1900;
	case 22:  goto L2200;
	case 23:  goto L2300;
	case 24:  goto L2400;
	case 25:  goto L1900;
	case 26:  goto L2600;
	case 27:  goto L3300;
	case 28:  goto L3300;
	case 29:  goto L3300;
	case 30:  goto L3300;
	case 31:  goto L3300;
	case 32:  goto L3300;
	case 33:  goto L3300;
	case 34:  goto L3400;
	case 35:  goto L3500;
	case 36:  goto L3600;
	case 37:  goto L3700;
	case 38:  goto L3800;
	case 39:  goto L3900;
	case 40:  goto L4000;
    }
/*                                        . . . . . . . . . . seek */
L200:
    mnseek_((S_fp)fcn, (U_fp)futil);
    goto L5000;
/*                                        . . . . . . . . . . simplex */
L300:
    mnsimp_((S_fp)fcn, (U_fp)futil);
    if (mn7flg_1.isw[3] < 1) {
	*ierflg = 4;
    }
    goto L5000;
/*                                        . . . . . . migrad, minimize */
L400:
    nf = mn7cnv_1.nfcn;
    mn7min_1.apsi = mn7min_1.epsi;
    mnmigr_((S_fp)fcn, (U_fp)futil);
    mnwerr_();
    if (mn7flg_1.isw[3] >= 1) {
	goto L5000;
    }
    *ierflg = 4;
    if (mn7flg_1.isw[0] == 1) {
	goto L5000;
    }
    if (s_cmp(mn7tit_1.cword, "MIG", (ftnlen)3, (ftnlen)3) == 0) {
	goto L5000;
    }
    mn7cnv_1.nfcnmx = mn7cnv_1.nfcnmx + nf - mn7cnv_1.nfcn;
    nf = mn7cnv_1.nfcn;
    mnsimp_((S_fp)fcn, (U_fp)futil);
    if (mn7flg_1.isw[0] == 1) {
	goto L5000;
    }
    mn7cnv_1.nfcnmx = mn7cnv_1.nfcnmx + nf - mn7cnv_1.nfcn;
    mnmigr_((S_fp)fcn, (U_fp)futil);
    if (mn7flg_1.isw[3] >= 1) {
	*ierflg = 0;
    }
    mnwerr_();
    goto L5000;
/*                                        . . . . . . . . . . minos */
L500:
    nsuper = mn7cnv_1.nfcn + (mn7npr_1.npar + 1 << 1) * mn7cnv_1.nfcnmx;
/*          possible loop over new minima */
    mn7min_1.epsi = mn7min_1.up * (float).1;
L510:
    mncuve_((S_fp)fcn, (U_fp)futil);
    mnmnos_((S_fp)fcn, (U_fp)futil);
    if (! mn7log_1.lnewmn) {
	goto L5000;
    }
    mnrset_(&c__0);
    mnmigr_((S_fp)fcn, (U_fp)futil);
    mnwerr_();
    if (mn7cnv_1.nfcn < nsuper) {
	goto L510;
    }
    io___343.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___343);
    e_wsfe();
    *ierflg = 4;
    goto L5000;
/*                                        . . . . . . . . . .set, show */
L700:
    mnset_((S_fp)fcn, (U_fp)futil);
    goto L5000;
/*                                        . . . . . . . . . . top of page */
L800:
    io___344.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___344);
    e_wsfe();
    goto L5000;
/*                                        . . . . . . . . . . fix */
L900:
    ltofix = TRUE_;
/*                                        . . (also release) .... */
L901:
    lfreed = FALSE_;
    lfixed = FALSE_;
    if (*llist == 0) {
	io___348.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___348);
	do_fio(&c__1, mn7tit_1.cword, (ftnlen)20);
	do_fio(&c__1, ":  NO PARAMETERS REQUESTED ", (ftnlen)27);
	e_wsfe();
	goto L5000;
    }
    i__1 = *llist;
    for (ilist = 1; ilist <= i__1; ++ilist) {
	iext = (integer) plist[ilist];
	s_copy(chwhy, " IS UNDEFINED.", (ftnlen)18, (ftnlen)14);
	if (iext <= 0) {
	    goto L930;
	}
	if (iext > mn7npr_1.nu) {
	    goto L930;
	}
	if (mn7inx_1.nvarl[iext - 1] < 0) {
	    goto L930;
	}
	s_copy(chwhy, " IS CONSTANT.  ", (ftnlen)18, (ftnlen)15);
	if (mn7inx_1.nvarl[iext - 1] == 0) {
	    goto L930;
	}
	iint = mn7inx_1.niofex[iext - 1];
	if (ltofix) {
	    s_copy(chwhy, " ALREADY FIXED.", (ftnlen)18, (ftnlen)15);
	    if (iint == 0) {
		goto L930;
	    }
	    mnfixp_(&iint, &ierr);
	    if (ierr == 0) {
		lfixed = TRUE_;
	    } else {
		*ierflg = 4;
	    }
	} else {
	    s_copy(chwhy, " ALREADY VARIABLE.", (ftnlen)18, (ftnlen)18);
	    if (iint > 0) {
		goto L930;
	    }
	    krl = -abs(iext);
	    mnfree_(&krl);
	    lfreed = TRUE_;
	}
	goto L950;
L930:
	io___355.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___355);
	do_fio(&c__1, " PARAMETER", (ftnlen)10);
	do_fio(&c__1, (char *)&iext, (ftnlen)sizeof(integer));
	do_fio(&c__1, chwhy, (ftnlen)18);
	do_fio(&c__1, " IGNORED.", (ftnlen)9);
	e_wsfe();
L950:
	;
    }
    if (lfreed || lfixed) {
	mnrset_(&c__0);
    }
    if (lfreed) {
	mn7flg_1.isw[1] = 0;
	mn7min_1.dcovar = (float)1.;
	mn7min_1.edm = mn7cns_1.bigedm;
	mn7flg_1.isw[3] = 0;
    }
    mnwerr_();
    if (mn7flg_1.isw[4] > 1) {
	mnprin_(&c__5, &mn7min_1.amin);
    }
    goto L5000;
/*                                        . . . . . . . . . . restore */
L1000:
    it = (integer) mn7arg_1.word7[0];
    if (it > 1 || it < 0) {
	goto L1005;
    }
    lfreed = mn7fx1_1.npfix > 0;
    mnfree_(&it);
    if (lfreed) {
	mnrset_(&c__0);
	mn7flg_1.isw[1] = 0;
	mn7min_1.dcovar = (float)1.;
	mn7min_1.edm = mn7cns_1.bigedm;
    }
    goto L5000;
L1005:
    io___357.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___357);
    do_fio(&c__1, " IGNORED.  UNKNOWN ARGUMENT:", (ftnlen)28);
    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
    e_wsfe();
    *ierflg = 3;
    goto L5000;
/*                                        . . . . . . . . . . release */
L1100:
    ltofix = FALSE_;
    goto L901;
/*                                       . . . . . . . . . . scan . . . */
L1200:
    iext = (integer) mn7arg_1.word7[0];
    if (iext <= 0) {
	goto L1210;
    }
    it2 = 0;
    if (iext <= mn7npr_1.nu) {
	it2 = mn7inx_1.niofex[iext - 1];
    }
    if (it2 <= 0) {
	goto L1250;
    }
L1210:
    mnscan_((S_fp)fcn, (U_fp)futil);
    goto L5000;
L1250:
    io___359.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___359);
    do_fio(&c__1, " PARAMETER", (ftnlen)10);
    do_fio(&c__1, (char *)&iext, (ftnlen)sizeof(integer));
    do_fio(&c__1, " NOT VARIABLE.", (ftnlen)14);
    e_wsfe();
    *ierflg = 3;
    goto L5000;
/*                                        . . . . . . . . . . contour */
L1300:
    ke1 = (integer) mn7arg_1.word7[0];
    ke2 = (integer) mn7arg_1.word7[1];
    if (ke1 == 0) {
	if (mn7npr_1.npar == 2) {
	    ke1 = mn7inx_1.nexofi[0];
	    ke2 = mn7inx_1.nexofi[1];
	} else {
	    io___362.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___362);
	    do_fio(&c__1, mn7tit_1.cword, (ftnlen)20);
	    do_fio(&c__1, ":  NO PARAMETERS REQUESTED ", (ftnlen)27);
	    e_wsfe();
	    *ierflg = 3;
	    goto L5000;
	}
    }
    mn7cnv_1.nfcnmx = 1000;
    mncntr_((S_fp)fcn, &ke1, &ke2, &ierrf, (U_fp)futil);
    if (ierrf > 0) {
	*ierflg = 3;
    }
    goto L5000;
/*                                        . . . . . . . . . . hesse */
L1400:
    mnhess_((S_fp)fcn, (U_fp)futil);
    mnwerr_();
    if (mn7flg_1.isw[4] >= 0) {
	mnprin_(&c__2, &mn7min_1.amin);
    }
    if (mn7flg_1.isw[4] >= 1) {
	mnmatu_(&c__1);
    }
    goto L5000;
/*                                        . . . . . . . . . . save */
L1500:
    mnsave_();
    goto L5000;
/*                                        . . . . . . . . . . improve */
L1600:
    mncuve_((S_fp)fcn, (U_fp)futil);
    mnimpr_((S_fp)fcn, (U_fp)futil);
    if (mn7log_1.lnewmn) {
	goto L400;
    }
    *ierflg = 4;
    goto L5000;
/*                                        . . . . . . . . . . call fcn */
L1700:
    iflag = (integer) mn7arg_1.word7[0];
    nparx = mn7npr_1.npar;
    f = mn7cns_1.undefi;
    (*fcn)(&nparx, mn7der_1.gin, &f, mn7ext_1.u, &iflag, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    nowprt = 0;
    if (f != mn7cns_1.undefi) {
	if (mn7min_1.amin == mn7cns_1.undefi) {
	    mn7min_1.amin = f;
	    nowprt = 1;
	} else if (f < mn7min_1.amin) {
	    mn7min_1.amin = f;
	    nowprt = 1;
	}
	if (mn7flg_1.isw[4] >= 0 && iflag <= 5 && nowprt == 1) {
	    mnprin_(&c__5, &mn7min_1.amin);
	}
	if (iflag == 3) {
	    mn7min_1.fval3 = f;
	}
    }
    if (iflag > 5) {
	mnrset_(&c__1);
    }
    goto L5000;
/*                                        . . . . . . . . . . standard */
L1800:
    stand_();
    goto L5000;
/*                                       . . . return, stop, end, exit */
L1900:
    it = (integer) mn7arg_1.word7[0];
    if (mn7min_1.fval3 != mn7min_1.amin && it == 0) {
	iflag = 3;
	if (mn7flg_1.isw[4] >= 0) {
	    io___368.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___368);
	    do_fio(&c__1, " CALL TO USER FUNCTION WITH IFLAG = 3", (ftnlen)37)
		    ;
	    e_wsfe();
	}
	nparx = mn7npr_1.npar;
	(*fcn)(&nparx, mn7der_1.gin, &f, mn7ext_1.u, &iflag, (U_fp)futil);
	++mn7cnv_1.nfcn;
	mn7min_1.fval3 = f;
    }
    *ierflg = 11;
    if (s_cmp(mn7tit_1.cword, "END", (ftnlen)3, (ftnlen)3) == 0) {
	*ierflg = 10;
    }
    if (s_cmp(mn7tit_1.cword, "RET", (ftnlen)3, (ftnlen)3) == 0) {
	*ierflg = 12;
    }
    goto L5000;
/*                                        . . . . . . . . . . clear */
L2200:
    mncler_();
    if (mn7flg_1.isw[4] >= 1) {
	io___369.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___369);
	do_fio(&c__1, " MINUIT MEMORY CLEARED. NO PARAMETERS NOW DEFINED.", (
		ftnlen)50);
	e_wsfe();
    }
    goto L5000;
/*                                        . . . . . . . . . . help */
L2300:
/* CCC      IF (INDEX(CWORD,'SHO') .GT. 0)  GO TO 700 */
/* CCC      IF (INDEX(CWORD,'SET') .GT. 0)  GO TO 700 */
    kcol = 0;
    i__1 = lk;
    for (icol = 5; icol <= i__1; ++icol) {
	if (*(unsigned char *)&mn7tit_1.cword[icol - 1] == ' ') {
	    goto L2310;
	}
	kcol = icol;
	goto L2320;
L2310:
	;
    }
L2320:
    if (kcol == 0) {
	s_copy(comd, "*   ", (ftnlen)4, (ftnlen)4);
    } else {
	s_copy(comd, mn7tit_1.cword + (kcol - 1), (ftnlen)4, lk - (kcol - 1));
    }
    mnhelp_(comd, &mn7iou_1.isyswr, (ftnlen)4);
    goto L5000;
/*                                       . . . . . . . . . . MNContour */
L2400:
    mn7min_1.epsi = mn7min_1.up * (float).05;
    ke1 = (integer) mn7arg_1.word7[0];
    ke2 = (integer) mn7arg_1.word7[1];
    if (ke1 == 0 && mn7npr_1.npar == 2) {
	ke1 = mn7inx_1.nexofi[0];
	ke2 = mn7inx_1.nexofi[1];
    }
    nptu = (integer) mn7arg_1.word7[2];
    if (nptu <= 0) {
	nptu = 20;
    }
    if (nptu > 101) {
	nptu = 101;
    }
    mn7cnv_1.nfcnmx = (nptu + 5) * 100 * (mn7npr_1.npar + 1);
    mncont_((S_fp)fcn, &ke1, &ke2, &nptu, xptu, yptu, &ierrf, (U_fp)futil);
    if (ierrf < nptu) {
	*ierflg = 4;
    }
    if (ierrf == -1) {
	*ierflg = 3;
    }
    goto L5000;
/*                                      . . . . . . . . . . jump */
L2600:
    step = mn7arg_1.word7[0];
    if (step <= 0.) {
	step = (float)2.;
    }
    rno = (float)0.;
    izero = 0;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mnrn15_(&rno, &izero);
	rno = rno * (float)2. - (float)1.;
/* L2620: */
	mn7int_1.x[i__ - 1] += rno * step * mn7err_1.werr[i__ - 1];
    }
    mninex_(mn7int_1.x);
    mnamin_((S_fp)fcn, (U_fp)futil);
    mnrset_(&c__0);
    goto L5000;
/*                                      . . . . . . . . . . blank line */
L3300:
    io___378.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___378);
    do_fio(&c__1, " BLANK COMMAND IGNORED.", (ftnlen)23);
    e_wsfe();
    *ierflg = 1;
    goto L5000;
/*  . . . . . . . . obsolete commands     . . . . . . . . . . . . . . */
/*                                      . . . . . . . . . . covariance */
L3400:
    io___379.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___379);
    do_fio(&c__1, " THE \"COVARIANCE\" COMMAND IS OSBSOLETE.", (ftnlen)39);
    do_fio(&c__1, " THE COVARIANCE MATRIX IS NOW SAVED IN A DIFFERENT FORMAT",
	     (ftnlen)57);
    do_fio(&c__1, " WITH THE \"SAVE\" COMMAND AND READ IN WITH:\"SET COVARIA\
NCE\"", (ftnlen)58);
    e_wsfe();
    *ierflg = 3;
    goto L5000;
/*                                        . . . . . . . . . . printout */
L3500:
    s_copy(cneway, "SET PRInt ", (ftnlen)10, (ftnlen)10);
    goto L3100;
/*                                        . . . . . . . . . . gradient */
L3600:
    s_copy(cneway, "SET GRAd  ", (ftnlen)10, (ftnlen)10);
    goto L3100;
/*                                        . . . . . . . . . . matout */
L3700:
    s_copy(cneway, "SHOW COVar", (ftnlen)10, (ftnlen)10);
    goto L3100;
/*                                        . . . . . . . . . error def */
L3800:
    s_copy(cneway, "SET ERRdef", (ftnlen)10, (ftnlen)10);
    goto L3100;
/*                                        . . . . . . . . . . limits */
L3900:
    s_copy(cneway, "SET LIMits", (ftnlen)10, (ftnlen)10);
    goto L3100;
/*                                        . . . . . . . . . . punch */
L4000:
    s_copy(cneway, "SAVE      ", (ftnlen)10, (ftnlen)10);
/*                                ....... come from obsolete commands */
L3100:
    io___381.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___381);
    do_fio(&c__1, mn7tit_1.cword, (ftnlen)20);
    do_fio(&c__1, cneway, (ftnlen)10);
    e_wsfe();
    s_copy(mn7tit_1.cword, cneway, (ftnlen)20, (ftnlen)10);
    if (s_cmp(mn7tit_1.cword, "SAVE      ", (ftnlen)20, (ftnlen)10) == 0) {
	goto L1500;
    }
    goto L700;
/*                                 . . . . . . . . . . . . . . . . . . */
L5000:
    return 0;
} /* mnexcm_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mnexin_(doublereal *pint)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iint, iext;
    static doublereal pinti;
    extern /* Subroutine */ int mnpint_(doublereal *, integer *, doublereal *)
	    ;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Transforms the external parameter values U to internal */
/* C        values in the dense array PINT. Subroutine MNPINT is used. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    /* Parameter adjustments */
    --pint;

    /* Function Body */
    mn7log_1.limset = FALSE_;
    i__1 = mn7npr_1.npar;
    for (iint = 1; iint <= i__1; ++iint) {
	iext = mn7inx_1.nexofi[iint - 1];
	mnpint_(&mn7ext_1.u[iext - 1], &iext, &pinti);
	pint[iint] = pinti;
/* L100: */
    }
    return 0;
} /* mnexin_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mnfixp_(integer *iint, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer i__, j, m, n, lc, ik;
    static doublereal yy[100];
    static integer kold, nold, ndex, knew, iext;
    static doublereal yyover;

    /* Fortran I/O blocks */
    static cilist io___385 = { 0, 0, 0, "(A,I4)", 0 };
    static cilist io___387 = { 0, 0, 0, "(A,I4,A,I4)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        removes parameter IINT from the internal (variable) parameter */
/* C        list, and arranges the rest of the list to fill the hole. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/*                           first see if it can be done */
    *ierr = 0;
    if (*iint > mn7npr_1.npar || *iint <= 0) {
	*ierr = 1;
	io___385.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___385);
	do_fio(&c__1, " MINUIT ERROR.  ARGUMENT TO MNFIXP=", (ftnlen)35);
	do_fio(&c__1, (char *)&(*iint), (ftnlen)sizeof(integer));
	e_wsfe();
	goto L300;
    }
    iext = mn7inx_1.nexofi[*iint - 1];
    if (mn7fx1_1.npfix >= 100) {
	*ierr = 1;
	io___387.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___387);
	do_fio(&c__1, " MINUIT CANNOT FIX PARAMETER", (ftnlen)28);
	do_fio(&c__1, (char *)&iext, (ftnlen)sizeof(integer));
	do_fio(&c__1, " MAXIMUM NUMBER THAT CAN BE FIXED IS", (ftnlen)36);
	do_fio(&c__1, (char *)&c__100, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L300;
    }
/*                           reduce number of variable parameters by one */
    mn7inx_1.niofex[iext - 1] = 0;
    nold = mn7npr_1.npar;
    --mn7npr_1.npar;
/*                       save values in case parameter is later restored */
    ++mn7fx1_1.npfix;
    mn7fx1_1.ipfix[mn7fx1_1.npfix - 1] = iext;
    lc = *iint;
    mn7fx2_1.xs[mn7fx1_1.npfix - 1] = mn7int_1.x[lc - 1];
    mn7fx2_1.xts[mn7fx1_1.npfix - 1] = mn7int_1.xt[lc - 1];
    mn7fx2_1.dirins[mn7fx1_1.npfix - 1] = mn7err_1.werr[lc - 1];
    mn7fx3_1.grds[mn7fx1_1.npfix - 1] = mn7der_1.grd[lc - 1];
    mn7fx3_1.g2s[mn7fx1_1.npfix - 1] = mn7der_1.g2[lc - 1];
    mn7fx3_1.gsteps[mn7fx1_1.npfix - 1] = mn7der_1.gstep[lc - 1];
/*                        shift values for other parameters to fill hole */
    i__1 = mn7npr_1.nu;
    for (ik = iext + 1; ik <= i__1; ++ik) {
	if (mn7inx_1.niofex[ik - 1] > 0) {
	    lc = mn7inx_1.niofex[ik - 1] - 1;
	    mn7inx_1.niofex[ik - 1] = lc;
	    mn7inx_1.nexofi[lc - 1] = ik;
	    mn7int_1.x[lc - 1] = mn7int_1.x[lc];
	    mn7int_1.xt[lc - 1] = mn7int_1.xt[lc];
	    mn7int_1.dirin[lc - 1] = mn7int_1.dirin[lc];
	    mn7err_1.werr[lc - 1] = mn7err_1.werr[lc];
	    mn7der_1.grd[lc - 1] = mn7der_1.grd[lc];
	    mn7der_1.g2[lc - 1] = mn7der_1.g2[lc];
	    mn7der_1.gstep[lc - 1] = mn7der_1.gstep[lc];
	}
/* L100: */
    }
    if (mn7flg_1.isw[1] <= 0) {
	goto L300;
    }
/*                    remove one row and one column from variance matrix */
    if (mn7npr_1.npar <= 0) {
	goto L300;
    }
    i__1 = nold;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = max(i__,*iint);
	n = min(i__,*iint);
	ndex = m * (m - 1) / 2 + n;
/* L260: */
	yy[i__ - 1] = mn7var_1.vhmat[ndex - 1];
    }
    yyover = (float)1. / yy[*iint - 1];
    knew = 0;
    kold = 0;
    i__1 = nold;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ++kold;
	    if (j == *iint || i__ == *iint) {
		goto L292;
	    }
	    ++knew;
	    mn7var_1.vhmat[knew - 1] = mn7var_1.vhmat[kold - 1] - yy[j - 1] * 
		    yy[i__ - 1] * yyover;
L292:
	    ;
	}
/* L294: */
    }
L300:
    return 0;
} /* mnfixp_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mnfree_(integer *k)
{
    /* Format strings */
    static char fmt_510[] = "(\002 CALL TO MNFREE IGNORED.  ARGUMENT GREATER\
 THAN ONE\002/)";
    static char fmt_500[] = "(\002 CALL TO MNFREE IGNORED.  THERE ARE NO FIX\
ED PA\002,\002RAMETERS\002/)";
    static char fmt_540[] = "(\002 IGNORED.  PARAMETER SPECIFIED IS ALREADY \
VARIABLE.\002)";
    static char fmt_530[] = "(\002 PARAMETER\002,i4,\002 NOT FIXED.  CANNOT \
BE RELEASED.\002)";
    static char fmt_520[] = "(20x,\002PARAMETER\002,i4,\002, \002,a10,\002 R\
ESTORED TO VARIABLE.\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, ka, lc, ik, iq, ir, is;
    static doublereal xv, g2v, xtv, grdv;
    static integer ipsav;
    static doublereal dirinv, gstepv;
    extern /* Subroutine */ int mnexin_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___400 = { 0, 0, 0, fmt_510, 0 };
    static cilist io___401 = { 0, 0, 0, fmt_500, 0 };
    static cilist io___403 = { 0, 0, 0, fmt_540, 0 };
    static cilist io___405 = { 0, 0, 0, fmt_530, 0 };
    static cilist io___418 = { 0, 0, 0, fmt_520, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Restores one or more fixed parameter(s) to variable status */
/* C        by inserting it into the internal parameter list at the */
/* C        appropriate place. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/* --       K = 0 means restore all parameters */
/* --       K = 1 means restore the last parameter fixed */
/* --       K = -I means restore external parameter I (if possible) */
/* --       IQ = fix-location where internal parameters were stored */
/* --       IR = external number of parameter being restored */
/* --       IS = internal number of parameter being restored */
    if (*k > 1) {
	io___400.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___400);
	e_wsfe();
    }
    if (mn7fx1_1.npfix < 1) {
	io___401.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___401);
	e_wsfe();
    }
    if (*k == 1 || *k == 0) {
	goto L40;
    }
/*                   release parameter with specified external number */
    ka = abs(*k);
    if (mn7inx_1.niofex[ka - 1] == 0) {
	goto L15;
    }
    io___403.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___403);
    e_wsfe();
    return 0;
L15:
    if (mn7fx1_1.npfix < 1) {
	goto L21;
    }
    i__1 = mn7fx1_1.npfix;
    for (ik = 1; ik <= i__1; ++ik) {
	if (mn7fx1_1.ipfix[ik - 1] == ka) {
	    goto L24;
	}
/* L20: */
    }
L21:
    io___405.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___405);
    do_fio(&c__1, (char *)&ka, (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;
L24:
    if (ik == mn7fx1_1.npfix) {
	goto L40;
    }
/*                   move specified parameter to end of list */
    ipsav = ka;
    xv = mn7fx2_1.xs[ik - 1];
    xtv = mn7fx2_1.xts[ik - 1];
    dirinv = mn7fx2_1.dirins[ik - 1];
    grdv = mn7fx3_1.grds[ik - 1];
    g2v = mn7fx3_1.g2s[ik - 1];
    gstepv = mn7fx3_1.gsteps[ik - 1];
    i__1 = mn7fx1_1.npfix;
    for (i__ = ik + 1; i__ <= i__1; ++i__) {
	mn7fx1_1.ipfix[i__ - 2] = mn7fx1_1.ipfix[i__ - 1];
	mn7fx2_1.xs[i__ - 2] = mn7fx2_1.xs[i__ - 1];
	mn7fx2_1.xts[i__ - 2] = mn7fx2_1.xts[i__ - 1];
	mn7fx2_1.dirins[i__ - 2] = mn7fx2_1.dirins[i__ - 1];
	mn7fx3_1.grds[i__ - 2] = mn7fx3_1.grds[i__ - 1];
	mn7fx3_1.g2s[i__ - 2] = mn7fx3_1.g2s[i__ - 1];
	mn7fx3_1.gsteps[i__ - 2] = mn7fx3_1.gsteps[i__ - 1];
/* L30: */
    }
    mn7fx1_1.ipfix[mn7fx1_1.npfix - 1] = ipsav;
    mn7fx2_1.xs[mn7fx1_1.npfix - 1] = xv;
    mn7fx2_1.xts[mn7fx1_1.npfix - 1] = xtv;
    mn7fx2_1.dirins[mn7fx1_1.npfix - 1] = dirinv;
    mn7fx3_1.grds[mn7fx1_1.npfix - 1] = grdv;
    mn7fx3_1.g2s[mn7fx1_1.npfix - 1] = g2v;
    mn7fx3_1.gsteps[mn7fx1_1.npfix - 1] = gstepv;
/*                restore last parameter in fixed list  -- IPFIX(NPFIX) */
L40:
    if (mn7fx1_1.npfix < 1) {
	goto L300;
    }
    ir = mn7fx1_1.ipfix[mn7fx1_1.npfix - 1];
    is = 0;
    i__1 = ir;
    for (ik = mn7npr_1.nu; ik >= i__1; --ik) {
	if (mn7inx_1.niofex[ik - 1] > 0) {
	    lc = mn7inx_1.niofex[ik - 1] + 1;
	    is = lc - 1;
	    mn7inx_1.niofex[ik - 1] = lc;
	    mn7inx_1.nexofi[lc - 1] = ik;
	    mn7int_1.x[lc - 1] = mn7int_1.x[lc - 2];
	    mn7int_1.xt[lc - 1] = mn7int_1.xt[lc - 2];
	    mn7int_1.dirin[lc - 1] = mn7int_1.dirin[lc - 2];
	    mn7err_1.werr[lc - 1] = mn7err_1.werr[lc - 2];
	    mn7der_1.grd[lc - 1] = mn7der_1.grd[lc - 2];
	    mn7der_1.g2[lc - 1] = mn7der_1.g2[lc - 2];
	    mn7der_1.gstep[lc - 1] = mn7der_1.gstep[lc - 2];
	}
/* L100: */
    }
    ++mn7npr_1.npar;
    if (is == 0) {
	is = mn7npr_1.npar;
    }
    mn7inx_1.niofex[ir - 1] = is;
    mn7inx_1.nexofi[is - 1] = ir;
    iq = mn7fx1_1.npfix;
    mn7int_1.x[is - 1] = mn7fx2_1.xs[iq - 1];
    mn7int_1.xt[is - 1] = mn7fx2_1.xts[iq - 1];
    mn7int_1.dirin[is - 1] = mn7fx2_1.dirins[iq - 1];
    mn7err_1.werr[is - 1] = mn7fx2_1.dirins[iq - 1];
    mn7der_1.grd[is - 1] = mn7fx3_1.grds[iq - 1];
    mn7der_1.g2[is - 1] = mn7fx3_1.g2s[iq - 1];
    mn7der_1.gstep[is - 1] = mn7fx3_1.gsteps[iq - 1];
    --mn7fx1_1.npfix;
    mn7flg_1.isw[1] = 0;
    mn7min_1.dcovar = (float)1.;
    if (mn7flg_1.isw[4] - mn7cnv_1.itaur >= 1) {
	io___418.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___418);
	do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
	do_fio(&c__1, mn7nam_1.cpnam + (ir - 1) * 10, (ftnlen)10);
	e_wsfe();
    }
    if (*k == 0) {
	goto L40;
    }
L300:
/*         if different from internal, external values are taken */
    mnexin_(mn7int_1.x);
/* L400: */
    return 0;
} /* mnfree_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:19  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mngrad_(S_fp fcn, U_fp futil)
{
    /* Format strings */
    static char fmt_51[] = "(/\002 CHECK OF GRADIENT CALCULATION IN FCN\002/\
12x,\002PARAMETER\002,6x,\002G(IN FCN)\002,3x,\002G(MINUIT)\002,2x,\002DG(MI\
NUIT)\002,3x,\002AGREEMENT\002)";
    static char fmt_99[] = "(7x,i5,2x,a10,3e12.4,4x,a4)";
    static char fmt_1003[] = "(/\002 MINUIT DOES NOT ACCEPT DERIVATIVE CALCU\
LATIONS BY FCN\002/\002 TO FORCE ACCEPTANCE, ENTER \"SET GRAD    1\"\002/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), do_fio(integer *, char *, 
	    ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal gf[100];
    static integer lc;
    static char cwd[4];
    static doublereal err;
    static logical lnone;
    static integer nparx;
    static doublereal fzero;
    extern /* Subroutine */ int mnhes1_(S_fp, U_fp), mnderi_(S_fp, U_fp), 
	    mninex_(doublereal *);
    static integer istsav;

    /* Fortran I/O blocks */
    static cilist io___424 = { 0, 0, 0, fmt_51, 0 };
    static cilist io___429 = { 0, 0, 0, fmt_99, 0 };
    static cilist io___430 = { 0, 0, 0, "(A)", 0 };
    static cilist io___431 = { 0, 0, 0, fmt_1003, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       Called from MNSET */
/* C       Interprets the SET GRAD command, which informs MINUIT whether */
/* C       the first derivatives of FCN will be calculated by the user */
/* C       inside FCN.  It can check the user's derivative calculation */
/* C       by comparing it with a finite difference approximation. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */




    mn7flg_1.isw[2] = 1;
    nparx = mn7npr_1.npar;
    if (mn7arg_1.word7[0] > 0.) {
	goto L2000;
    }
/*                  get user-calculated first derivatives from FCN */
    i__1 = mn7npr_1.nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	mn7der_1.gin[i__ - 1] = mn7cns_1.undefi;
    }
    mninex_(mn7int_1.x);
    (*fcn)(&nparx, mn7der_1.gin, &fzero, mn7ext_1.u, &c__2, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    mnderi_((S_fp)fcn, (U_fp)futil);
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	gf[i__ - 1] = mn7der_1.grd[i__ - 1];
    }
/*                    get MINUIT-calculated first derivatives */
    mn7flg_1.isw[2] = 0;
    istsav = mn7cnv_1.istrat;
    mn7cnv_1.istrat = 2;
    mnhes1_((S_fp)fcn, (U_fp)futil);
    mn7cnv_1.istrat = istsav;
    io___424.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___424);
    e_wsfe();
    mn7flg_1.isw[2] = 1;
    lnone = FALSE_;
    i__1 = mn7npr_1.npar;
    for (lc = 1; lc <= i__1; ++lc) {
	i__ = mn7inx_1.nexofi[lc - 1];
	s_copy(cwd, "GOOD", (ftnlen)4, (ftnlen)4);
	err = mn7der_1.dgrd[lc - 1];
	if ((d__1 = gf[lc - 1] - mn7der_1.grd[lc - 1], abs(d__1)) > err) {
	    s_copy(cwd, " BAD", (ftnlen)4, (ftnlen)4);
	}
	if (mn7der_1.gin[i__ - 1] == mn7cns_1.undefi) {
	    s_copy(cwd, "NONE", (ftnlen)4, (ftnlen)4);
	    lnone = TRUE_;
	    gf[lc - 1] = (float)0.;
	}
	if (s_cmp(cwd, "GOOD", (ftnlen)4, (ftnlen)4) != 0) {
	    mn7flg_1.isw[2] = 0;
	}
	io___429.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___429);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, mn7nam_1.cpnam + (i__ - 1) * 10, (ftnlen)10);
	do_fio(&c__1, (char *)&gf[lc - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&mn7der_1.grd[lc - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&err, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, cwd, (ftnlen)4);
	e_wsfe();
/* L100: */
    }
    if (lnone) {
	io___430.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___430);
	do_fio(&c__1, "  AGREEMENT=NONE  MEANS FCN DID NOT CALCULATE THE DER\
IVATIVE", (ftnlen)60);
	e_wsfe();
    }
    if (mn7flg_1.isw[2] == 0) {
	io___431.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___431);
	e_wsfe();
    }

L2000:
    return 0;
} /* mngrad_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1999/09/03 09:17:47  couet */
/* - \Cind{} removed in the help of minuit. This was a Tex directive which very */
/*   likely has been forgotten during a Tex to f77 translation. This didn't */
/*   compile on RH6. */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnhelp_(char *comd, integer *lout, ftnlen comd_len)
{
    /* Format strings */
    static char fmt_10000[] = "(\002   ==>List of MINUIT Interactive command\
s:\002,/,\002 CLEar     Reset all parameter names and values undefined\002,/,\
\002 CONtour   Make contour map of the user function\002,/,\002 EXIT      Ex\
it from Interactive Minuit\002,/,\002 FIX       Cause parameter(s) to remain\
 constant\002,/,\002 HESse     Calculate the Hessian or error matrix.\002,/\
,\002 IMPROVE   Search for a new minimum around current minimum\002,/,\002 M\
IGrad    Minimize by the method of Migrad\002,/,\002 MINImize  MIGRAD + SIMP\
LEX method if Migrad fails\002,/,\002 MINOs     Exact (non-linear) parameter\
 error analysis\002)";
    static char fmt_10001[] = "(\002 MNContour Calculate one MINOS function \
contour\002,/,\002 PARameter Define or redefine new parameters and values\
\002,/,\002 RELease   Make previously FIXed parameters variable again\002,/\
,\002 REStore   Release last parameter fixed\002,/,\002 SAVe      Save curre\
nt parameter values on a file\002,/,\002 SCAn      Scan the user function by\
 varying parameters\002,/,\002 SEEk      Minimize by the method of Monte Car\
lo\002,/,\002 SET       Set various MINUIT constants or conditions\002,/,\
\002 SHOw      Show values of current constants or conditions\002,/,\002 SIM\
plex   Minimize by the method of Simplex\002)";
    static char fmt_10100[] = "(\002 ***>CLEAR\002,/,\002 Resets all paramet\
er names and values to undefined.\002,/,\002 Must normally be followed by a \
PARameters command or \002,/,\002 equivalent, in order to define parameter v\
alues.\002)";
    static char fmt_10200[] = "(\002 ***>CONTOUR <par1>  <par2>  [devs]  [ng\
rid]\002,/,\002 Instructs Minuit to trace contour lines of the user functio\
n\002,/,\002 with respect to the two parameters whose external numbers\002,/,\
\002 are <par1> and <par2>.\002,/,\002 Other variable parameters of the func\
tion, if any, will have\002,/,\002 their values fixed at the current values \
during the contour\002,/,\002 tracing. The optional parameter [devs] (defaul\
t value 2.)\002,/,\002 gives the number of standard deviations in each param\
eter\002,/,\002 which should lie entirely within the plotting area.\002,/\
,\002 Optional parameter [ngrid] (default value 25 unless page\002,/,\002 si\
ze is too small) determines the resolution of the plot,\002,/,\002 i.e. the \
number of rows and columns of the grid at which the\002,/,\002 function will\
 be evaluated. [See also MNContour.]\002)";
    static char fmt_10300[] = "(\002 ***>END\002,/,\002 Signals the end of a\
 data block (i.e., the end of a fit),\002,/,\002 and implies that execution \
should continue, because another\002,/,\002 Data Block follows. A Data Block\
 is a set of Minuit data\002,/,\002 consisting of\002,/,\002     (1) A Title,\
\002,/,\002     (2) One or more Parameter Definitions,\002,/,\002     (3) A \
blank line, and\002,/,\002     (4) A set of Minuit Commands.\002,/,\002 The \
END command is used when more than one Data Block is to\002,/,\002 be used w\
ith the same FCN function. It first causes Minuit\002,/,\002 to issue a CALL\
 FCN with IFLAG=3, in order to allow FCN to\002,/,\002 perform any calculati\
ons associated with the final fitted\002,/,\002 parameter values, unless a C\
ALL FCN 3 command has already\002,/,\002 been executed at the current FCN va\
lue.\002)";
    static char fmt_10400[] = "(\002 ***>EXIT\002,/,\002 Signals the end of \
execution.\002,/,\002 The EXIT command first causes Minuit to issue a CALL F\
CN\002,/,\002 with IFLAG=3, to allow FCN to perform any calculations\002,/\
,\002 associated with the final fitted parameter values, unless a\002,/,\002\
 CALL FCN 3 command has already been executed.\002)";
    static char fmt_10500[] = "(\002 ***>FIX} <parno> [parno] ... [parno]\
\002,/,\002 Causes parameter(s) <parno> to be removed from the list of\002,/,\
\002 variable parameters, and their value(s) will remain constant\002,/,\002\
 during subsequent minimizations, etc., until another command\002,/,\002 cha\
nges their value(s) or status.\002)";
    static char fmt_10600[] = "(\002 ***>HESse  [maxcalls]\002,/,\002 Calcul\
ate, by finite differences, the Hessian or error matrix.\002,/,\002  That is\
, it calculates the full matrix of second derivatives\002,/,\002 of the func\
tion with respect to the currently variable\002,/,\002 parameters, and inver\
ts it, printing out the resulting error\002,/,\002 matrix. The optional argu\
ment [maxcalls] specifies the\002,/,\002 (approximate) maximum number of fun\
ction calls after which\002,/,\002 the calculation will be stopped.\002)";
    static char fmt_10700[] = "(\002 ***>IMPROVE  [maxcalls]\002,/,\002 If a\
 previous minimization has converged, and the current\002,/,\002 values of t\
he parameters therefore correspond to a local\002,/,\002 minimum of the func\
tion, this command requests a search for\002,/,\002 additional distinct loca\
l minima.\002,/,\002 The optional argument [maxcalls] specifies the (approxi\
mate)\002,/,\002 maximum number of function calls after which the calculation\
\002,/,\002 will be stopped.\002)";
    static char fmt_10800[] = "(\002 ***>MIGrad  [maxcalls]  [tolerance]\002\
,/,\002 Causes minimization of the function by the method of Migrad,\002,/\
,\002 the most efficient and complete single method, recommended\002,/,\002 \
for general functions (see also MINImize).\002,/,\002 The minimization produ\
ces as a by-product the error matrix\002,/,\002 of the parameters, which is \
usually reliable unless warning\002,/,\002 messages are produced.\002,/,\002\
 The optional argument [maxcalls] specifies the (approximate)\002,/,\002 max\
imum number of function calls after which the calculation\002,/,\002 will be\
 stopped even if it has not yet converged.\002,/,\002 The optional argument \
[tolerance] specifies required tolerance\002,/,\002 on the function value at\
 the minimum.\002,/,\002 The default tolerance is 0.1, and the minimization \
will stop\002,/,\002 when the estimated vertical distance to the minimum (ED\
M) is\002,/,\002 less than 0.001*[tolerance]*UP (see [SET ERRordef]).\002)";
    static char fmt_10900[] = "(\002 ***>MINImize  [maxcalls] [tolerance]\
\002,/,\002 Causes minimization of the function by the method of Migrad,\002\
,/,\002 as does the MIGrad command, but switches to the SIMplex method\002,/,\
\002 if Migrad fails to converge. Arguments are as for MIGrad.\002,/,\002 No\
te that command requires four characters to be unambiguous.\002)";
    static char fmt_11000[] = "(\002 ***>MINOs  [maxcalls]  [parno] [parno] \
...\002,/,\002 Causes a Minos error analysis to be performed on the paramete\
rs\002,/,\002 whose numbers [parno] are specified. If none are specified,\
\002,/,\002 Minos errors are calculated for all variable parameters.\002,/\
,\002 Minos errors may be expensive to calculate, but are very\002,/,\002 re\
liable since they take account of non-linearities in the\002,/,\002 problem \
as well as parameter correlations, and are in general\002,/\002 asymmetric\
.\002,/,\002 The optional argument [maxcalls] specifies the (approximate)\
\002,/,\002 maximum number of function calls per parameter requested,\002,/\
,\002 after which the calculation will stop for that parameter.\002)";
    static char fmt_11100[] = "(\002 ***>MNContour  <par1> <par2> [npts]\002\
,/,\002 Calculates one function contour of FCN with respect to\002,/,\002 pa\
rameters par1 and par2, with FCN minimized always with\002,/,\002 respect to\
 all other NPAR-2 variable parameters (if any).\002,/,\002 Minuit will try t\
o find npts points on the contour (default 20)\002,/,\002 If only two parame\
ters are variable at the time, it is not\002,/,\002 necessary to specify the\
ir numbers. To calculate more than\002,/,\002 one contour, it is necessary t\
o SET ERRordef to the appropriate\002,/,\002 value and issue the MNContour c\
ommand for each contour.\002)";
    static char fmt_11150[] = "(\002 ***>PARameters\002,/,\002 followed by o\
ne or more parameter definitions.\002,/,\002 Parameter definitions are of th\
e form:\002,/,\002   <number>  'name'  <value>  <step>  [lolim] [uplim] \002\
,/,\002 for example:\002,/,\002  3  'K width'  1.2   0.1\002,/,\002 the last\
 definition is followed by a blank line or a zero.\002)";
    static char fmt_11200[] = "(\002 ***>RELease  <parno> [parno] ... [par\
no]\002,/,\002 If <parno> is the number of a previously variable paramete\
r\002,/,\002 which has been fixed by a command: FIX <parno>, then that\002,/,\
\002 parameter will return to variable status.  Otherwise a warning\002,/\
,\002 message is printed and the command is ignored.\002,/,\002 Note that th\
is command operates only on parameters which were\002,/\002 at one time vari\
able and have been FIXed. It cannot make\002,/,\002 constant parameters vari\
able; that must be done by redefining\002,/\002 the parameter with a PARamet\
ers command.\002)";
    static char fmt_11300[] = "(\002 ***>REStore  [code]\002,/,\002 If no [c\
ode] is specified, this command restores all previously\002,/,\002 FIXed par\
ameters to variable status. If [code]=1, then only\002,/,\002 the last param\
eter FIXed is restored to variable status.\002,/,\002 If code is neither zer\
o nor one, the command is ignored.\002)";
    static char fmt_11400[] = "(\002 ***>RETURN\002,/,\002 Signals the end o\
f a data block, and instructs Minuit to return\002,/,\002 to the program whi\
ch called it. The RETurn command first\002,/,\002 causes Minuit to CALL FCN \
with IFLAG=3, in order to allow FCN\002,/,\002 to perform any calculations a\
ssociated with the final fitted\002,/,\002 parameter values, unless a CALL F\
CN 3 command has already been\002,/,\002 executed at the current FCN value\
.\002)";
    static char fmt_11500[] = "(\002 ***>SAVe\002,/,\002 Causes the current \
parameter values to be saved on a file in\002,/,\002 such a format that they\
 can be read in again as Minuit\002,/,\002 parameter definitions. If the cov\
ariance matrix exists, it is\002,/,\002 also output in such a format. The un\
it number is by default 7,\002,/,\002 or that specified by the user in his c\
all to MINTIO or\002,/,\002 MNINIT. The user is responsible for opening the \
file previous\002,/,\002 to issuing the [SAVe] command (except where this ca\
n be done\002,/,\002 interactively).\002)";
    static char fmt_11600[] = "(\002 ***>SCAn  [parno]  [numpts] [from]  [\
to]\002,/,\002 Scans the value of the user function by varying parameter\002\
,/,\002 number [parno], leaving all other parameters fixed at the\002,/,\002\
 current value. If [parno] is not specified, all variable\002,/,\002 paramet\
ers are scanned in sequence.\002,/,\002 The number of points [numpts] in the\
 scan is 40 by default,\002,/,\002 and cannot exceed 100. The range of the s\
can is by default\002,/,\002 2 standard deviations on each side of the curre\
nt best value,\002,/,\002 but can be specified as from [from] to [to].\002,/,\
\002 After each scan, if a new minimum is found, the best parameter\002,/\
,\002 values are retained as start values for future scans or\002,/,\002 min\
imizations. The curve resulting from each scan is plotted\002,/,\002 on the \
output unit in order to show the approximate behaviour\002,/,\002 of the fun\
ction.\002,/,\002 This command is not intended for minimization, but is some\
times\002,/,\002 useful for debugging the user function or finding a\002,/\
,\002 reasonable starting point.\002)";
    static char fmt_11700[] = "(\002 ***>SEEk  [maxcalls]  [devs]\002,/,\002\
 Causes a Monte Carlo minimization of the function, by choosing\002,/,\002 r\
andom values of the variable parameters, chosen uniformly\002,/,\002 over a \
hypercube centered at the current best value.\002,/,\002 The region size is \
by default 3 standard deviations on each\002,/,\002 side, but can be changed\
 by specifying the value of [devs].\002)";
    static char fmt_11800[] = "(\002 ***>SET <option_name>\002,/,/,\002  SET\
 BATch\002,/,\002    Informs Minuit that it is running in batch mode.\002,//,\
\002  SET EPSmachine  <accuracy>\002,/,\002    Informs Minuit that the relat\
ive floating point arithmetic\002,/\002    precision is <accuracy>. Minuit d\
etermines the nominal\002,/,\002    precision itself, but the SET EPSmachine\
 command can be\002,/,\002    used to override Minuit own determination, whe\
n the user\002,/,\002    knows that the FCN function value is not calculated\
 to\002,/,\002    the nominal machine accuracy. Typical values of <accuracy\
>\002,/\002    are between 10**-5 and 10**-14.\002)";
    static char fmt_11801[] = "(/,\002  SET ERRordef  <up>\002,/,\002    Set\
s the value of UP (default value= 1.), defining parameter\002,/,\002    erro\
rs. Minuit defines parameter errors as the change\002,/,\002    in parameter\
 value required to change the function value\002,/,\002    by UP. Normally, \
for chisquared fits UP=1, and for negative\002,/,\002    log likelihood, UP=\
0.5.\002)";
    static char fmt_11802[] = "(/,\002   SET GRAdient  [force]\002,/,\002   \
 Informs Minuit that the user function is prepared to\002,/,\002    calculat\
e its own first derivatives and return their values\002,/,\002    in the arr\
ay GRAD when IFLAG=2 (see specs of FCN).\002,/,\002    If [force] is not spe\
cified, Minuit will calculate\002,/,\002    the FCN derivatives by finite di\
fferences at the current\002,/,\002    point and compare with the user calcu\
lation at that point,\002,/,\002    accepting the user values only if they a\
gree.\002,/,\002    If [force]=1, Minuit does not do its own derivative\002,\
/,\002    calculation, and uses the derivatives calculated in FCN.\002)";
    static char fmt_11803[] = "(/,\002   SET INPut  [unitno]  [filename]\002\
,/,\002    Causes Minuit, in data-driven mode only, to read subsequent\002,/,\
\002    commands (or parameter definitions) from a different input\002,/,\
\002    file. If no [unitno] is specified, reading reverts to the\002,/,\002\
    previous input file, assuming that there was one.\002,/,\002    If [unit\
no] is specified, and that unit has not been opened,\002,/,\002    then Minu\
it attempts to open the file [filename]} if a\002,/,\002    name is specifie\
d. If running in interactive mode and\002,/,\002    [filename] is not specif\
ied and [unitno] is not opened,\002,/,\002    Minuit prompts the user to ent\
er a file name.\002,/,\002    If the word REWIND is added to the command (no\
te:no blanks\002,/\002    between INPUT and REWIND), the file is rewound bef\
ore\002,/,\002    reading. Note that this command is implemented in standar\
d\002,/\002    Fortran 77 and the results may depend on the  system;\002,/\
,\002    for example, if a filename is given under VM/CMS, it must\002,/,\
\002    be preceeded by a slash.\002)";
    static char fmt_11804[] = "(/,\002   SET INTeractive\002,/,\002    Infor\
ms Minuit that it is running interactively.\002)";
    static char fmt_11805[] = "(/,\002   SET LIMits  [parno]  [lolim]  [upli\
m]\002,/,\002    Allows the user to change the limits on one or all\002,/\
,\002    parameters. If no arguments are specified, all limits are\002,/,\
\002    removed from all parameters. If [parno] alone is specified,\002,/\
,\002    limits are removed from parameter [parno].\002,/,\002    If all arg\
uments are specified, then parameter [parno] will\002,/,\002    be bounded b\
etween [lolim] and [uplim].\002,/,\002    Limits can be specified in either \
order, Minuit will take\002,/,\002    the smaller as [lolim] and the larger \
as [uplim].\002,/,\002    However, if [lolim] is equal to [uplim], an error \
condition\002,/,\002    results.\002)";
    static char fmt_11806[] = "(/,\002   SET LINesperpage\002,/,\002     Set\
s the number of lines for one page of output.\002,/,\002     Default value i\
s 24 for interactive mode\002)";
    static char fmt_11807[] = "(/,\002   SET NOGradient\002,/,\002    The in\
verse of SET GRAdient, instructs Minuit not to\002,/,\002    use the first d\
erivatives calculated by the user in FCN.\002)";
    static char fmt_11808[] = "(/,\002   SET NOWarnings\002,/,\002    Supres\
ses Minuit warning messages.\002)";
    static char fmt_11809[] = "(/,\002   SET OUTputfile  <unitno>\002,/,\002\
    Instructs Minuit to write further output to unit <unitno>.\002)";
    static char fmt_11810[] = "(/,\002   SET PAGethrow  <integer>\002,/,\002\
    Sets the carriage control character for ``new page' to\002,/,\002    <in\
teger>. Thus the value 1 produces a new page, and 0\002,/,\002    produces a\
 blank line, on some devices (see TOPofpage)\002)";
    static char fmt_11811[] = "(/,\002   SET PARameter  <parno>  <value>\002\
,/,\002    Sets the value of parameter <parno> to <value>.\002,/,\002    The\
 parameter in question may be variable, fixed, or\002,/,\002    constant, bu\
t must be defined.\002)";
    static char fmt_11812[] = "(/,\002   SET PRIntout  <level>\002,/,\002   \
 Sets the print level, determining how much output will be\002,/,\002    pro\
duced. Allowed values and their meanings are displayed\002,/,\002    after a\
 SHOw PRInt command, and are currently <level>=:\002,/,\002      [-1]  no ou\
tput except from SHOW commands\002,/,\002       [0]  minimum output\002,/\
,\002       [1]  default value, normal output\002,/,\002       [2]  addition\
al output giving intermediate results.\002,/,\002       [3]  maximum output,\
 showing progress of minimizations.\002,/\002    Note: See also the SET WARn\
ings command.\002)";
    static char fmt_11813[] = "(/,\002   SET RANdomgenerator  <seed>\002,/\
,\002    Sets the seed of the random number generator used in SEEk.\002,/\
\002    This can be any integer between 10000 and 900000000, for\002,/,\002 \
   example one which was output from a SHOw RANdom command of\002,/\002    a\
 previous run.\002)";
    static char fmt_11814[] = "(/,\002   SET STRategy  <level>\002,/,\002   \
 Sets the strategy to be used in calculating first and second\002,/,\002    \
derivatives and in certain minimization methods.\002,/,\002    In general, l\
ow values of <level> mean fewer function calls\002,/,\002    and high values\
 mean more reliable minimization.\002,/,\002    Currently allowed values are\
 0, 1 (default), and 2.\002)";
    static char fmt_11815[] = "(/,\002   SET TITle\002,/,\002    Informs Min\
uit that the next input line is to be considered\002,/,\002    the (new) tit\
le for this task or sub-task.  This is for\002,/,\002    the convenience of \
the user in reading his output.\002)";
    static char fmt_11816[] = "(/,\002   SET WARnings\002,/,\002    Instruct\
s Minuit to output warning messages when suspicious\002,/,\002    conditions\
 arise which may indicate unreliable results.\002,/\002    This is the defau\
lt.\002)";
    static char fmt_11817[] = "(/,\002    SET WIDthpage\002,/,\002    Inform\
s Minuit of the output page width.\002,/,\002    Default values are 80 for i\
nteractive jobs\002)";
    static char fmt_11900[] = "(\002 ***>SHOw  <option_name>\002,/,\002  All\
 SET XXXX commands have a corresponding SHOw XXXX command.\002,/,\002  In ad\
dition, the SHOw commands listed starting here have no\002,/,\002  correspon\
ding SET command for obvious reasons.\002)";
    static char fmt_11901[] = "(/,\002   SHOw CORrelations\002,/,\002    Cal\
culates and prints the parameter correlations from the\002,/,\002    error m\
atrix.\002)";
    static char fmt_11902[] = "(/,\002   SHOw COVariance\002,/,\002    Print\
s the (external) covariance (error) matrix.\002)";
    static char fmt_11903[] = "(/,\002   SHOw EIGenvalues\002,/,\002    Calc\
ulates and prints the eigenvalues of the covariance\002,/,\002    matrix.\
\002)";
    static char fmt_11904[] = "(/,\002   SHOw FCNvalue\002,/,\002    Prints \
the current value of FCN.\002)";
    static char fmt_12000[] = "(\002 ***>SIMplex  [maxcalls]  [tolerance]\
\002,/,\002 Performs a function minimization using the simplex method of\002\
,/\002 Nelder and Mead. Minimization terminates either when the\002,/,\002 f\
unction has been called (approximately) [maxcalls] times,\002,/,\002 or when\
 the estimated vertical distance to minimum (EDM) is\002,/,\002 less than [t\
olerance].\002,/,\002 The default value of [tolerance] is 0.1*UP(see SET ERR\
ordef).\002)";
    static char fmt_12100[] = "(\002 ***>STAndard\002,/,\002 Causes Minuit t\
o execute the Fortran instruction CALL STAND\002,/,\002 where STAND is a sub\
routine supplied by the user.\002)";
    static char fmt_12200[] = "(\002 ***>STOP\002,/,\002 Same as EXIT.\002)";
    static char fmt_12300[] = "(\002 ***>TOPofpage\002,/,\002 Causes Minuit \
to write the character specified in a\002,/,\002 SET PAGethrow command (defa\
ult = 1) to column 1 of the output\002,/,\002 file, which may or may not pos\
ition your output medium to\002,/,\002 the top of a page depending on the de\
vice and system.\002)";
    static char fmt_13000[] = "(\002 Unknown MINUIT command. Type HELP for l\
ist of commands.\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char cmd3[3];

    /* Fortran I/O blocks */
    static cilist io___432 = { 0, 0, 0, fmt_10000, 0 };
    static cilist io___433 = { 0, 0, 0, fmt_10001, 0 };
    static cilist io___435 = { 0, 0, 0, fmt_10100, 0 };
    static cilist io___436 = { 0, 0, 0, fmt_10200, 0 };
    static cilist io___437 = { 0, 0, 0, fmt_10300, 0 };
    static cilist io___438 = { 0, 0, 0, fmt_10400, 0 };
    static cilist io___439 = { 0, 0, 0, fmt_10500, 0 };
    static cilist io___440 = { 0, 0, 0, fmt_10600, 0 };
    static cilist io___441 = { 0, 0, 0, fmt_10700, 0 };
    static cilist io___442 = { 0, 0, 0, fmt_10800, 0 };
    static cilist io___443 = { 0, 0, 0, fmt_10900, 0 };
    static cilist io___444 = { 0, 0, 0, fmt_11000, 0 };
    static cilist io___445 = { 0, 0, 0, fmt_11100, 0 };
    static cilist io___446 = { 0, 0, 0, fmt_11150, 0 };
    static cilist io___447 = { 0, 0, 0, fmt_11200, 0 };
    static cilist io___448 = { 0, 0, 0, fmt_11300, 0 };
    static cilist io___449 = { 0, 0, 0, fmt_11400, 0 };
    static cilist io___450 = { 0, 0, 0, fmt_11500, 0 };
    static cilist io___451 = { 0, 0, 0, fmt_11600, 0 };
    static cilist io___452 = { 0, 0, 0, fmt_11700, 0 };
    static cilist io___453 = { 0, 0, 0, fmt_11800, 0 };
    static cilist io___454 = { 0, 0, 0, fmt_11801, 0 };
    static cilist io___455 = { 0, 0, 0, fmt_11802, 0 };
    static cilist io___456 = { 0, 0, 0, fmt_11803, 0 };
    static cilist io___457 = { 0, 0, 0, fmt_11804, 0 };
    static cilist io___458 = { 0, 0, 0, fmt_11805, 0 };
    static cilist io___459 = { 0, 0, 0, fmt_11806, 0 };
    static cilist io___460 = { 0, 0, 0, fmt_11807, 0 };
    static cilist io___461 = { 0, 0, 0, fmt_11808, 0 };
    static cilist io___462 = { 0, 0, 0, fmt_11809, 0 };
    static cilist io___463 = { 0, 0, 0, fmt_11810, 0 };
    static cilist io___464 = { 0, 0, 0, fmt_11811, 0 };
    static cilist io___465 = { 0, 0, 0, fmt_11812, 0 };
    static cilist io___466 = { 0, 0, 0, fmt_11813, 0 };
    static cilist io___467 = { 0, 0, 0, fmt_11814, 0 };
    static cilist io___468 = { 0, 0, 0, fmt_11815, 0 };
    static cilist io___469 = { 0, 0, 0, fmt_11816, 0 };
    static cilist io___470 = { 0, 0, 0, fmt_11817, 0 };
    static cilist io___471 = { 0, 0, 0, fmt_11900, 0 };
    static cilist io___472 = { 0, 0, 0, fmt_11901, 0 };
    static cilist io___473 = { 0, 0, 0, fmt_11902, 0 };
    static cilist io___474 = { 0, 0, 0, fmt_11903, 0 };
    static cilist io___475 = { 0, 0, 0, fmt_11904, 0 };
    static cilist io___476 = { 0, 0, 0, fmt_12000, 0 };
    static cilist io___477 = { 0, 0, 0, fmt_12100, 0 };
    static cilist io___478 = { 0, 0, 0, fmt_12200, 0 };
    static cilist io___479 = { 0, 0, 0, fmt_12300, 0 };
    static cilist io___480 = { 0, 0, 0, fmt_13000, 0 };


/* . */
/* .         HELP routine for MINUIT interactive commands. */
/* . */
/* .      COMD ='*   '  prints a global help for all commands */
/* .      COMD =Command_name: print detailed help for one command. */
/* .          Note that at least 3 characters must be given for the command name. */
/* . */
/* .     Author: Rene Brun */
/*             comments extracted from the MINUIT documentation file. */
/* . */
/* . */
/* -- command name ASSUMED to be in upper case */
/* __________________________________________________________________ */
/* -- */
/* --  Global HELP: Summary of all commands */
/* --  ==================================== */
/* -- */
    if (*(unsigned char *)comd == '*') {
	io___432.ciunit = *lout;
	s_wsfe(&io___432);
	e_wsfe();
	io___433.ciunit = *lout;
	s_wsfe(&io___433);
	e_wsfe();
	goto L99;
    }

    s_copy(cmd3, comd, (ftnlen)3, (ftnlen)3);
/* __________________________________________________________________ */
/* -- */
/* --  Command CLEAR */
/* --  ============= */
/* . */
    if (s_cmp(cmd3, "CLE", (ftnlen)3, (ftnlen)3) == 0) {
	io___435.ciunit = *lout;
	s_wsfe(&io___435);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command CONTOUR */
/* --  =============== */
/* . */
    if (s_cmp(cmd3, "CON", (ftnlen)3, (ftnlen)3) == 0) {
	io___436.ciunit = *lout;
	s_wsfe(&io___436);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command END */
/* --  =========== */
/* . */
    if (s_cmp(cmd3, "END", (ftnlen)3, (ftnlen)3) == 0) {
	io___437.ciunit = *lout;
	s_wsfe(&io___437);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* . */
/* -- */
/* --  Command EXIT */
/* --  ============ */
    if (s_cmp(cmd3, "EXI", (ftnlen)3, (ftnlen)3) == 0) {
	io___438.ciunit = *lout;
	s_wsfe(&io___438);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command FIX */
/* --  =========== */
/* . */
    if (s_cmp(cmd3, "FIX", (ftnlen)3, (ftnlen)3) == 0) {
	io___439.ciunit = *lout;
	s_wsfe(&io___439);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command HESSE */
/* --  ============= */
/* . */
    if (s_cmp(cmd3, "HES", (ftnlen)3, (ftnlen)3) == 0) {
	io___440.ciunit = *lout;
	s_wsfe(&io___440);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command IMPROVE */
/* --  =============== */
/* . */
    if (s_cmp(cmd3, "IMP", (ftnlen)3, (ftnlen)3) == 0) {
	io___441.ciunit = *lout;
	s_wsfe(&io___441);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command MIGRAD */
/* --  ============== */
/* . */
    if (s_cmp(cmd3, "MIG", (ftnlen)3, (ftnlen)3) == 0) {
	io___442.ciunit = *lout;
	s_wsfe(&io___442);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command MINIMIZE */
/* --  ================ */
/* . */
    if (s_cmp(comd, "MINI", (ftnlen)4, (ftnlen)4) == 0) {
	io___443.ciunit = *lout;
	s_wsfe(&io___443);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command MINOS */
/* --  ============= */
/* . */
    if (s_cmp(comd, "MINO", (ftnlen)4, (ftnlen)4) == 0) {
	io___444.ciunit = *lout;
	s_wsfe(&io___444);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command MNCONTOUR */
/* --  ================= */
/* . */
    if (s_cmp(cmd3, "MNC", (ftnlen)3, (ftnlen)3) == 0) {
	io___445.ciunit = *lout;
	s_wsfe(&io___445);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command PARAMETER */
/* --  ================= */
/* . */
    if (s_cmp(cmd3, "PAR", (ftnlen)3, (ftnlen)3) == 0) {
	io___446.ciunit = *lout;
	s_wsfe(&io___446);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command RELEASE */
/* --  =============== */
/* . */
    if (s_cmp(cmd3, "REL", (ftnlen)3, (ftnlen)3) == 0) {
	io___447.ciunit = *lout;
	s_wsfe(&io___447);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command RESTORE */
/* --  =============== */
/* . */
    if (s_cmp(cmd3, "RES", (ftnlen)3, (ftnlen)3) == 0) {
	io___448.ciunit = *lout;
	s_wsfe(&io___448);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command RETURN */
/* --  ============== */
/* . */
    if (s_cmp(cmd3, "RET", (ftnlen)3, (ftnlen)3) == 0) {
	io___449.ciunit = *lout;
	s_wsfe(&io___449);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command SAVE */
/* --  ============ */
/* . */
    if (s_cmp(cmd3, "SAV", (ftnlen)3, (ftnlen)3) == 0) {
	io___450.ciunit = *lout;
	s_wsfe(&io___450);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command SCAN */
/* --  ============ */
/* . */
    if (s_cmp(cmd3, "SCA", (ftnlen)3, (ftnlen)3) == 0) {
	io___451.ciunit = *lout;
	s_wsfe(&io___451);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command SEEK */
/* --  ============ */
/* . */
    if (s_cmp(cmd3, "SEE", (ftnlen)3, (ftnlen)3) == 0) {
	io___452.ciunit = *lout;
	s_wsfe(&io___452);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command SET */
/* --  =========== */
/* . */
    if (s_cmp(cmd3, "SET", (ftnlen)3, (ftnlen)3) == 0) {
	io___453.ciunit = *lout;
	s_wsfe(&io___453);
	e_wsfe();
	io___454.ciunit = *lout;
	s_wsfe(&io___454);
	e_wsfe();
	io___455.ciunit = *lout;
	s_wsfe(&io___455);
	e_wsfe();
	io___456.ciunit = *lout;
	s_wsfe(&io___456);
	e_wsfe();
	io___457.ciunit = *lout;
	s_wsfe(&io___457);
	e_wsfe();
	io___458.ciunit = *lout;
	s_wsfe(&io___458);
	e_wsfe();
	io___459.ciunit = *lout;
	s_wsfe(&io___459);
	e_wsfe();
	io___460.ciunit = *lout;
	s_wsfe(&io___460);
	e_wsfe();
	io___461.ciunit = *lout;
	s_wsfe(&io___461);
	e_wsfe();
	io___462.ciunit = *lout;
	s_wsfe(&io___462);
	e_wsfe();
	io___463.ciunit = *lout;
	s_wsfe(&io___463);
	e_wsfe();
	io___464.ciunit = *lout;
	s_wsfe(&io___464);
	e_wsfe();
	io___465.ciunit = *lout;
	s_wsfe(&io___465);
	e_wsfe();
	io___466.ciunit = *lout;
	s_wsfe(&io___466);
	e_wsfe();
	io___467.ciunit = *lout;
	s_wsfe(&io___467);
	e_wsfe();
	io___468.ciunit = *lout;
	s_wsfe(&io___468);
	e_wsfe();
	io___469.ciunit = *lout;
	s_wsfe(&io___469);
	e_wsfe();
	io___470.ciunit = *lout;
	s_wsfe(&io___470);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command SHOW */
/* --  ============ */
/* . */
    if (s_cmp(cmd3, "SHO", (ftnlen)3, (ftnlen)3) == 0) {
	io___471.ciunit = *lout;
	s_wsfe(&io___471);
	e_wsfe();
	io___472.ciunit = *lout;
	s_wsfe(&io___472);
	e_wsfe();
	io___473.ciunit = *lout;
	s_wsfe(&io___473);
	e_wsfe();
	io___474.ciunit = *lout;
	s_wsfe(&io___474);
	e_wsfe();
	io___475.ciunit = *lout;
	s_wsfe(&io___475);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command SIMPLEX */
/* --  =============== */
/* . */
    if (s_cmp(cmd3, "SIM", (ftnlen)3, (ftnlen)3) == 0) {
	io___476.ciunit = *lout;
	s_wsfe(&io___476);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command STANDARD */
/* --  ================ */
/* . */
    if (s_cmp(cmd3, "STA", (ftnlen)3, (ftnlen)3) == 0) {
	io___477.ciunit = *lout;
	s_wsfe(&io___477);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command STOP */
/* --  ============ */
/* . */
    if (s_cmp(cmd3, "STO", (ftnlen)3, (ftnlen)3) == 0) {
	io___478.ciunit = *lout;
	s_wsfe(&io___478);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */
/* -- */
/* --  Command TOPOFPAGE */
/* --  ================= */
/* . */
    if (s_cmp(cmd3, "TOP", (ftnlen)3, (ftnlen)3) == 0) {
	io___479.ciunit = *lout;
	s_wsfe(&io___479);
	e_wsfe();
	goto L99;
    }
/* __________________________________________________________________ */

    io___480.ciunit = *lout;
    s_wsfe(&io___480);
    e_wsfe();

L99:
    return 0;
} /* mnhelp_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnhes1_(S_fp fcn, U_fp futil)
{
    /* Format strings */
    static char fmt_11[] = "(i4,i2,6g12.5)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    doublereal d__1, d__2, d__3;
    char ch__1[48];

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double d_sign(doublereal *, doublereal *);
    integer s_wsfi(icilist *), e_wsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static doublereal d__;
    static integer i__;
    static doublereal fs1, fs2, sag, xtf;
    static char cbf1[22];
    static doublereal dmin__;
    static integer icyc, ncyc, idrv;
    static doublereal dfmin, dgmin;
    static integer nparx;
    static doublereal change, chgold;
    static logical ldebug;
    static doublereal grdold, epspri;
    extern /* Subroutine */ int mninex_(doublereal *);
    static doublereal grdnew;
    extern /* Subroutine */ int mnwarn_(char *, char *, char *, ftnlen, 
	    ftnlen, ftnlen);
    static doublereal optstp;

    /* Fortran I/O blocks */
    static cilist io___500 = { 0, 0, 0, fmt_11, 0 };
    static icilist io___503 = { 0, cbf1, 0, "(2G11.3)", 22, 1 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C      Called from MNHESS and MNGRAD */
/* C      Calculate first derivatives (GRD) and uncertainties (DGRD) */
/* C         and appropriate step sizes GSTEP */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    ldebug = mn7flg_1.idbg[5] >= 1;
    if (mn7cnv_1.istrat <= 0) {
	ncyc = 1;
    }
    if (mn7cnv_1.istrat == 1) {
	ncyc = 2;
    }
    if (mn7cnv_1.istrat > 1) {
	ncyc = 6;
    }
    idrv = 1;
    nparx = mn7npr_1.npar;
    dfmin = mn7cns_1.epsma2 * (float)4. * (abs(mn7min_1.amin) + mn7min_1.up);
/*                                     main loop over parameters */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xtf = mn7int_1.x[i__ - 1];
	dmin__ = mn7cns_1.epsma2 * (float)4. * abs(xtf);
	epspri = mn7cns_1.epsma2 + (d__1 = mn7der_1.grd[i__ - 1] * 
		mn7cns_1.epsma2, abs(d__1));
	optstp = sqrt(dfmin / ((d__1 = mn7der_1.g2[i__ - 1], abs(d__1)) + 
		epspri));
	d__ = (d__1 = mn7der_1.gstep[i__ - 1], abs(d__1)) * (float).2;
	if (d__ > optstp) {
	    d__ = optstp;
	}
	if (d__ < dmin__) {
	    d__ = dmin__;
	}
	chgold = (float)1e4;
/*                                       iterate reducing step size */
	i__2 = ncyc;
	for (icyc = 1; icyc <= i__2; ++icyc) {
	    mn7int_1.x[i__ - 1] = xtf + d__;
	    mninex_(mn7int_1.x);
	    (*fcn)(&nparx, mn7der_1.gin, &fs1, mn7ext_1.u, &c__4, (U_fp)futil)
		    ;
	    ++mn7cnv_1.nfcn;
	    mn7int_1.x[i__ - 1] = xtf - d__;
	    mninex_(mn7int_1.x);
	    (*fcn)(&nparx, mn7der_1.gin, &fs2, mn7ext_1.u, &c__4, (U_fp)futil)
		    ;
	    ++mn7cnv_1.nfcn;
	    mn7int_1.x[i__ - 1] = xtf;
/*                                       check if step sizes appropriate */
	    sag = (fs1 + fs2 - mn7min_1.amin * (float)2.) * (float).5;
	    grdold = mn7der_1.grd[i__ - 1];
	    grdnew = (fs1 - fs2) / (d__ * (float)2.);
	    dgmin = mn7cns_1.epsmac * (abs(fs1) + abs(fs2)) / d__;
	    if (ldebug) {
		io___500.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___500);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&idrv, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&mn7der_1.gstep[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&d__, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&mn7der_1.g2[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&grdnew, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&sag, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (grdnew == 0.) {
		goto L60;
	    }
	    change = (d__1 = (grdold - grdnew) / grdnew, abs(d__1));
	    if (change > chgold && icyc > 1) {
		goto L60;
	    }
	    chgold = change;
	    mn7der_1.grd[i__ - 1] = grdnew;
	    mn7der_1.gstep[i__ - 1] = d_sign(&d__, &mn7der_1.gstep[i__ - 1]);
/*                  decrease step until first derivative changes by <5% */
	    if (change < (float).05) {
		goto L60;
	    }
	    if ((d__1 = grdold - grdnew, abs(d__1)) < dgmin) {
		goto L60;
	    }
	    if (d__ < dmin__) {
		mnwarn_("D", "MNHES1", "Step size too small for 1st drv.", (
			ftnlen)1, (ftnlen)6, (ftnlen)32);
		goto L60;
	    }
	    d__ *= (float).2;
/* L50: */
	}
/*                                       loop satisfied = too many iter */
	s_wsfi(&io___503);
	do_fio(&c__1, (char *)&grdold, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&grdnew, (ftnlen)sizeof(doublereal));
	e_wsfi();
/* Writing concatenation */
	i__3[0] = 26, a__1[0] = "Too many iterations on D1.";
	i__3[1] = 22, a__1[1] = cbf1;
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)48);
	mnwarn_("D", "MNHES1", ch__1, (ftnlen)1, (ftnlen)6, (ftnlen)48);
L60:
/* Computing MAX */
	d__2 = dgmin, d__3 = (d__1 = grdold - grdnew, abs(d__1));
	mn7der_1.dgrd[i__ - 1] = max(d__2,d__3);
/* L100: */
    }
/*                                        end of first deriv. loop */
    mninex_(mn7int_1.x);
    return 0;
} /* mnhes1_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnhess_(S_fp fcn, U_fp futil)
{
    /* Format strings */
    static char fmt_31[] = "(i4,i2,6g12.5)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    doublereal d__1, d__2;
    char ch__1[48], ch__2[41], ch__3[40], ch__4[45];

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), e_wsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal d__;
    static integer i__, j;
    static doublereal df;
    static integer id;
    static doublereal yy[100], g2i, fs1, fs2, sag, xtf, xti, xtj;
    static char cbf1[22];
    static doublereal dmin__, dxdi;
    static integer icyc;
    static doublereal elem;
    static integer ncyc, ndex, idrv, iext;
    static doublereal wint;
    static integer npar2;
    static doublereal tlrg2;
    static integer ifail, npard;
    static doublereal dlast;
    static integer nparx;
    static doublereal ztemp, g2bfor;
    extern /* Subroutine */ int mnhes1_(S_fp, U_fp);
    static doublereal aimsag;
    static logical ldebug;
    extern /* Subroutine */ int mnamin_(S_fp, U_fp), mndxdi_(doublereal *, 
	    integer *, doublereal *), mninex_(doublereal *), mnwarn_(char *, 
	    char *, char *, ftnlen, ftnlen, ftnlen);
    static doublereal stpinm, tlrstp;
    static integer multpy;
    extern /* Subroutine */ int mnpsdf_(), mnvert_(doublereal *, integer *, 
	    integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___508 = { 0, 0, 0, "(A)", 0 };
    static icilist io___514 = { 0, cbf1, 0, "(G12.3)", 12, 1 };
    static cilist io___515 = { 0, 0, 0, "(A,A)", 0 };
    static icilist io___522 = { 0, cbf1, 0, "(I4)", 4, 1 };
    static icilist io___532 = { 0, cbf1, 0, "(I4)", 4, 1 };
    static cilist io___534 = { 0, 0, 0, fmt_31, 0 };
    static icilist io___538 = { 0, cbf1, 0, "(I2,2E10.2)", 22, 1 };
    static cilist io___546 = { 0, 0, 0, "(A)", 0 };
    static cilist io___547 = { 0, 0, 0, "(A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Calculates the full second-derivative matrix of FCN */
/* C        by taking finite differences. When calculating diagonal */
/* C        elements, it may iterate so that step size is nearly that */
/* C        which gives function change= UP/10. The first derivatives */
/* C        of course come as a free side effect, but with a smaller */
/* C        step size in order to obtain a known accuracy. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



    ldebug = mn7flg_1.idbg[3] >= 1;
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    if (mn7cnv_1.istrat <= 0) {
	ncyc = 3;
	tlrstp = (float).5;
	tlrg2 = (float).1;
    } else if (mn7cnv_1.istrat == 1) {
	ncyc = 5;
	tlrstp = (float).3;
	tlrg2 = (float).05;
    } else {
	ncyc = 7;
	tlrstp = (float).1;
	tlrg2 = (float).02;
    }
    if (mn7flg_1.isw[4] >= 2 || ldebug) {
	io___508.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___508);
	do_fio(&c__1, "   START COVARIANCE MATRIX CALCULATION.", (ftnlen)39);
	e_wsfe();
    }
    s_copy(mn7tit_1.cfrom, "HESSE   ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "OK        ", (ftnlen)10, (ftnlen)10);
    npard = mn7npr_1.npar;
/*                 make sure starting at the right place */
    mninex_(mn7int_1.x);
    nparx = mn7npr_1.npar;
    (*fcn)(&nparx, mn7der_1.gin, &fs1, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    if (fs1 != mn7min_1.amin) {
	df = mn7min_1.amin - fs1;
	s_wsfi(&io___514);
	do_fio(&c__1, (char *)&df, (ftnlen)sizeof(doublereal));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 36, a__1[0] = "function value differs from AMIN by ";
	i__1[1] = 12, a__1[1] = cbf1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)48);
	mnwarn_("D", "MNHESS", ch__1, (ftnlen)1, (ftnlen)6, (ftnlen)48);
    }
    mn7min_1.amin = fs1;
    if (ldebug) {
	io___515.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___515);
	do_fio(&c__1, " PAR D   GSTEP          ", (ftnlen)24);
	do_fio(&c__1, " D          G2         GRD         SAG    ", (ftnlen)
		42);
	e_wsfe();
    }
/*                                        . . . . . . diagonal elements . */
/*         ISW(2) = 1 if approx, 2 if not posdef, 3 if ok */
/*         AIMSAG is the sagitta we are aiming for in second deriv calc. */
    aimsag = sqrt(mn7cns_1.epsma2) * (abs(mn7min_1.amin) + mn7min_1.up);
/*         Zero the second derivative matrix */
    npar2 = mn7npr_1.npar * (mn7npr_1.npar + 1) / 2;
    i__2 = npar2;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L10: */
	mn7var_1.vhmat[i__ - 1] = (float)0.;
    }

/*         Loop over variable parameters for second derivatives */
    idrv = 2;
    i__2 = npard;
    for (id = 1; id <= i__2; ++id) {
	i__ = id + mn7npr_1.npar - npard;
	iext = mn7inx_1.nexofi[i__ - 1];
	if (mn7der_1.g2[i__ - 1] == 0.) {
	    s_wsfi(&io___522);
	    do_fio(&c__1, (char *)&iext, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 37, a__1[0] = "Second derivative enters zero, param ";
	    i__1[1] = 4, a__1[1] = cbf1;
	    s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)41);
	    mnwarn_("W", "HESSE", ch__2, (ftnlen)1, (ftnlen)5, (ftnlen)41);
	    wint = mn7err_1.werr[i__ - 1];
	    if (mn7inx_1.nvarl[iext - 1] > 1) {
		mndxdi_(&mn7int_1.x[i__ - 1], &i__, &dxdi);
		if (abs(dxdi) < (float).001) {
		    wint = (float).01;
		} else {
		    wint /= abs(dxdi);
		}
	    }
/* Computing 2nd power */
	    d__1 = wint;
	    mn7der_1.g2[i__ - 1] = mn7min_1.up / (d__1 * d__1);
	}
	xtf = mn7int_1.x[i__ - 1];
	dmin__ = mn7cns_1.epsma2 * (float)8. * abs(xtf);

/*                               find step which gives sagitta = AIMSAG */
	d__ = (d__1 = mn7der_1.gstep[i__ - 1], abs(d__1));
	i__3 = ncyc;
	for (icyc = 1; icyc <= i__3; ++icyc) {
/*                               loop here only if SAG=0. */
	    for (multpy = 1; multpy <= 5; ++multpy) {
/*           take two steps */
		mn7int_1.x[i__ - 1] = xtf + d__;
		mninex_(mn7int_1.x);
		nparx = mn7npr_1.npar;
		(*fcn)(&nparx, mn7der_1.gin, &fs1, mn7ext_1.u, &c__4, (U_fp)
			futil);
		++mn7cnv_1.nfcn;
		mn7int_1.x[i__ - 1] = xtf - d__;
		mninex_(mn7int_1.x);
		(*fcn)(&nparx, mn7der_1.gin, &fs2, mn7ext_1.u, &c__4, (U_fp)
			futil);
		++mn7cnv_1.nfcn;
		mn7int_1.x[i__ - 1] = xtf;
		sag = (fs1 + fs2 - mn7min_1.amin * (float)2.) * (float).5;
		if (sag != 0.) {
		    goto L30;
		}
		if (mn7der_1.gstep[i__ - 1] < 0.) {
		    if (d__ >= (float).5) {
			goto L26;
		    }
		    d__ *= (float)10.;
		    if (d__ > (float).5) {
			d__ = (float).51;
		    }
		    goto L25;
		}
		d__ *= (float)10.;
L25:
		;
	    }
L26:
	    s_wsfi(&io___532);
	    do_fio(&c__1, (char *)&iext, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 36, a__1[0] = "Second derivative zero for parameter";
	    i__1[1] = 4, a__1[1] = cbf1;
	    s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)40);
	    mnwarn_("W", "HESSE", ch__3, (ftnlen)1, (ftnlen)5, (ftnlen)40);
	    goto L390;
/*                             SAG is not zero */
L30:
	    g2bfor = mn7der_1.g2[i__ - 1];
/* Computing 2nd power */
	    d__1 = d__;
	    mn7der_1.g2[i__ - 1] = sag * (float)2. / (d__1 * d__1);
	    mn7der_1.grd[i__ - 1] = (fs1 - fs2) / (d__ * (float)2.);
	    if (ldebug) {
		io___534.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___534);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&idrv, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&mn7der_1.gstep[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&d__, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&mn7der_1.g2[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&mn7der_1.grd[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&sag, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    mn7der_1.gstep[i__ - 1] = d_sign(&d__, &mn7der_1.gstep[i__ - 1]);
	    mn7int_1.dirin[i__ - 1] = d__;
	    yy[i__ - 1] = fs1;
	    dlast = d__;
	    d__ = sqrt(aimsag * (float)2. / (d__1 = mn7der_1.g2[i__ - 1], abs(
		    d__1)));
/*         if parameter has limits, max int step size = 0.5 */
	    stpinm = (float).5;
	    if (mn7der_1.gstep[i__ - 1] < 0.) {
		d__ = min(d__,stpinm);
	    }
	    if (d__ < dmin__) {
		d__ = dmin__;
	    }
/*           see if converged */
	    if ((d__1 = (d__ - dlast) / d__, abs(d__1)) < tlrstp) {
		goto L50;
	    }
	    if ((d__1 = (mn7der_1.g2[i__ - 1] - g2bfor) / mn7der_1.g2[i__ - 1]
		    , abs(d__1)) < tlrg2) {
		goto L50;
	    }
/* Computing MIN */
	    d__1 = d__, d__2 = dlast * (float)10.;
	    d__ = min(d__1,d__2);
/* Computing MAX */
	    d__1 = d__, d__2 = dlast * (float).1;
	    d__ = max(d__1,d__2);
/* L40: */
	}
/*                       end of step size loop */
	s_wsfi(&io___538);
	do_fio(&c__1, (char *)&iext, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&sag, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&aimsag, (ftnlen)sizeof(doublereal));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 23, a__1[0] = "Second Deriv. SAG,AIM= ";
	i__1[1] = 22, a__1[1] = cbf1;
	s_cat(ch__4, a__1, i__1, &c__2, (ftnlen)45);
	mnwarn_("D", "MNHESS", ch__4, (ftnlen)1, (ftnlen)6, (ftnlen)45);

L50:
	ndex = i__ * (i__ + 1) / 2;
	mn7var_1.vhmat[ndex - 1] = mn7der_1.g2[i__ - 1];
/* L100: */
    }
/*                              end of diagonal second derivative loop */
    mninex_(mn7int_1.x);
/*                                     refine the first derivatives */
    if (mn7cnv_1.istrat > 0) {
	mnhes1_((S_fp)fcn, (U_fp)futil);
    }
    mn7flg_1.isw[1] = 3;
    mn7min_1.dcovar = (float)0.;
/*                                        . . . .  off-diagonal elements */
    if (mn7npr_1.npar == 1) {
	goto L214;
    }
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__3 = i__ - 1;
	for (j = 1; j <= i__3; ++j) {
	    xti = mn7int_1.x[i__ - 1];
	    xtj = mn7int_1.x[j - 1];
	    mn7int_1.x[i__ - 1] = xti + mn7int_1.dirin[i__ - 1];
	    mn7int_1.x[j - 1] = xtj + mn7int_1.dirin[j - 1];
	    mninex_(mn7int_1.x);
	    (*fcn)(&nparx, mn7der_1.gin, &fs1, mn7ext_1.u, &c__4, (U_fp)futil)
		    ;
	    ++mn7cnv_1.nfcn;
	    mn7int_1.x[i__ - 1] = xti;
	    mn7int_1.x[j - 1] = xtj;
	    elem = (fs1 + mn7min_1.amin - yy[i__ - 1] - yy[j - 1]) / (
		    mn7int_1.dirin[i__ - 1] * mn7int_1.dirin[j - 1]);
	    ndex = i__ * (i__ - 1) / 2 + j;
	    mn7var_1.vhmat[ndex - 1] = elem;
/* L180: */
	}
/* L200: */
    }
L214:
    mninex_(mn7int_1.x);
/*                  verify matrix positive-definite */
    mnpsdf_();
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__3 = i__;
	for (j = 1; j <= i__3; ++j) {
	    ndex = i__ * (i__ - 1) / 2 + j;
	    mn7sim_1.p[i__ + j * 100 - 101] = mn7var_1.vhmat[ndex - 1];
/* L219: */
	    mn7sim_1.p[j + i__ * 100 - 101] = mn7sim_1.p[i__ + j * 100 - 101];
	}
/* L220: */
    }
    mnvert_(mn7sim_1.p, &mn7npr_1.maxint, &mn7npr_1.maxint, &mn7npr_1.npar, &
	    ifail);
    if (ifail > 0) {
	mnwarn_("W", "HESSE", "Matrix inversion fails.", (ftnlen)1, (ftnlen)5,
		 (ftnlen)23);
	goto L390;
    }
/*                                        . . . . . . .  calculate  e d m */
    mn7min_1.edm = (float)0.;
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*                              off-diagonal elements */
	ndex = i__ * (i__ - 1) / 2;
	i__3 = i__ - 1;
	for (j = 1; j <= i__3; ++j) {
	    ++ndex;
	    ztemp = mn7sim_1.p[i__ + j * 100 - 101] * (float)2.;
	    mn7min_1.edm += mn7der_1.grd[i__ - 1] * ztemp * mn7der_1.grd[j - 
		    1];
/* L225: */
	    mn7var_1.vhmat[ndex - 1] = ztemp;
	}
/*                              diagonal elements */
	++ndex;
	mn7var_1.vhmat[ndex - 1] = mn7sim_1.p[i__ + i__ * 100 - 101] * (float)
		2.;
/* Computing 2nd power */
	d__1 = mn7der_1.grd[i__ - 1];
	mn7min_1.edm += mn7sim_1.p[i__ + i__ * 100 - 101] * (d__1 * d__1);
/* L230: */
    }
    if (mn7flg_1.isw[4] >= 1 && mn7flg_1.isw[1] == 3 && mn7cnv_1.itaur == 0) {
	io___546.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___546);
	do_fio(&c__1, " COVARIANCE MATRIX CALCULATED SUCCESSFULLY", (ftnlen)
		42);
	e_wsfe();
    }
    goto L900;
/*                              failure to invert 2nd deriv matrix */
L390:
    mn7flg_1.isw[1] = 1;
    mn7min_1.dcovar = (float)1.;
    s_copy(mn7tit_1.cstatu, "FAILED    ", (ftnlen)10, (ftnlen)10);
    if (mn7flg_1.isw[4] >= 0) {
	io___547.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___547);
	do_fio(&c__1, "  MNHESS FAILS AND WILL RETURN DIAGONAL MATRIX. ", (
		ftnlen)48);
	e_wsfe();
    }
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ndex = i__ * (i__ - 1) / 2;
	i__3 = i__ - 1;
	for (j = 1; j <= i__3; ++j) {
	    ++ndex;
/* L394: */
	    mn7var_1.vhmat[ndex - 1] = (float)0.;
	}
	++ndex;
	g2i = mn7der_1.g2[i__ - 1];
	if (g2i <= 0.) {
	    g2i = (float)1.;
	}
/* L395: */
	mn7var_1.vhmat[ndex - 1] = (float)2. / g2i;
    }
L900:
    return 0;
} /* mnhess_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnimpr_(S_fp fcn, U_fp futil)
{
    /* Initialized data */

    static doublereal rnum = 0.;

    /* Format strings */
    static char fmt_1040[] = "(/\002 START ATTEMPT NO.\002,i2,\002 TO FIND N\
EW MINIMUM\002)";
    static char fmt_1000[] = "(\002 AN IMPROVEMENT ON THE PREVIOUS MINIMUM H\
AS BEEN FOUND\002)";
    static char fmt_1030[] = "(/\002 IMPROVE HAS FOUND A TRULY NEW MINIMU\
M\002/\002 \002,37(\002*\002)/)";
    static char fmt_1020[] = "(/\002 COVARIANCE MATRIX WAS NOT POSITIVE-DEFI\
NITE\002)";
    static char fmt_1010[] = "(\002 IMPROVE HAS RETURNED TO REGION OF ORIGIN\
AL MINIMUM\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer i__, j;
    static doublereal y[101];
    static integer jh;
    static doublereal pb, ep;
    static integer jl;
    static doublereal wg, xi, reg, sig2, amax, dsav[100];
    static integer npfn, ndex, loop, ifail, iseed;
    static doublereal ycalf;
    static integer jhold, nloop, nparx;
    extern /* Subroutine */ int mnrn15_(doublereal *, integer *);
    static doublereal ystar, ystst;
    static integer nparp1;
    extern /* Subroutine */ int mncalf_(S_fp, doublereal *, doublereal *, 
	    U_fp), mnamin_(S_fp, U_fp);
    static doublereal sigsav;
    extern /* Subroutine */ int mnvert_(doublereal *, integer *, integer *, 
	    integer *, integer *), mnrazz_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), mninex_(doublereal *), 
	    mnsimp_(S_fp, U_fp), mnrset_(integer *), mnprin_(integer *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___564 = { 0, 0, 0, fmt_1040, 0 };
    static cilist io___577 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___578 = { 0, 0, 0, fmt_1030, 0 };
    static cilist io___579 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___580 = { 0, 0, 0, fmt_1010, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Attempts to improve on a good local minimum by finding a */
/* C        better one.   The quadratic part of FCN is removed by MNCALF */
/* C        and this transformed function is minimized using the simplex */
/* C        method from several random starting points. */
/* C        ref. -- Goldstein and Price, Math.Comp. 25, 569 (1971) */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    if (mn7npr_1.npar <= 0) {
	return 0;
    }
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    s_copy(mn7tit_1.cstatu, "UNCHANGED ", (ftnlen)10, (ftnlen)10);
    mn7cnv_1.itaur = 1;
    mn7min_1.epsi = mn7min_1.up * (float).1;
    npfn = mn7cnv_1.nfcn;
    nloop = (integer) mn7arg_1.word7[1];
    if (nloop <= 0) {
	nloop = mn7npr_1.npar + 4;
    }
    nparx = mn7npr_1.npar;
    nparp1 = mn7npr_1.npar + 1;
    wg = (float)1. / (real) mn7npr_1.npar;
    sigsav = mn7min_1.edm;
    mn7min_1.apsi = mn7min_1.amin;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7int_1.xt[i__ - 1] = mn7int_1.x[i__ - 1];
	dsav[i__ - 1] = mn7err_1.werr[i__ - 1];
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ndex = i__ * (i__ - 1) / 2 + j;
	    mn7sim_1.p[i__ + j * 100 - 101] = mn7var_1.vhmat[ndex - 1];
/* L2: */
	    mn7sim_1.p[j + i__ * 100 - 101] = mn7sim_1.p[i__ + j * 100 - 101];
	}
    }
    mnvert_(mn7sim_1.p, &mn7npr_1.maxint, &mn7npr_1.maxint, &mn7npr_1.npar, &
	    ifail);
    if (ifail >= 1) {
	goto L280;
    }
/*               Save inverted matrix in VT */
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ndex = i__ * (i__ - 1) / 2;
	i__1 = i__;
	for (j = 1; j <= i__1; ++j) {
	    ++ndex;
/* L12: */
	    mn7vat_1.vthmat[ndex - 1] = mn7sim_1.p[i__ + j * 100 - 101];
	}
    }
    loop = 0;

L20:
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7int_1.dirin[i__ - 1] = dsav[i__ - 1] * (float)2.;
	mnrn15_(&rnum, &iseed);
/* L25: */
	mn7int_1.x[i__ - 1] = mn7int_1.xt[i__ - 1] + mn7int_1.dirin[i__ - 1] *
		 (float)2. * (rnum - (float).5);
    }
    ++loop;
    reg = (float)2.;
    if (mn7flg_1.isw[4] >= 0) {
	io___564.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___564);
	do_fio(&c__1, (char *)&loop, (ftnlen)sizeof(integer));
	e_wsfe();
    }
L30:
    mncalf_((S_fp)fcn, mn7int_1.x, &ycalf, (U_fp)futil);
    mn7min_1.amin = ycalf;
/*                                        . . . . set up  random simplex */
    jl = nparp1;
    jh = nparp1;
    y[nparp1 - 1] = mn7min_1.amin;
    amax = mn7min_1.amin;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = mn7int_1.x[i__ - 1];
	mnrn15_(&rnum, &iseed);
	mn7int_1.x[i__ - 1] = xi - mn7int_1.dirin[i__ - 1] * (rnum - (float)
		.5);
	mncalf_((S_fp)fcn, mn7int_1.x, &ycalf, (U_fp)futil);
	y[i__ - 1] = ycalf;
	if (y[i__ - 1] < mn7min_1.amin) {
	    mn7min_1.amin = y[i__ - 1];
	    jl = i__;
	} else if (y[i__ - 1] > amax) {
	    amax = y[i__ - 1];
	    jh = i__;
	}
	i__2 = mn7npr_1.npar;
	for (j = 1; j <= i__2; ++j) {
/* L40: */
	    mn7sim_1.p[j + i__ * 100 - 101] = mn7int_1.x[j - 1];
	}
	mn7sim_1.p[i__ + nparp1 * 100 - 101] = xi;
	mn7int_1.x[i__ - 1] = xi;
/* L45: */
    }

    mn7min_1.edm = mn7min_1.amin;
    sig2 = mn7min_1.edm;
/*                                        . . . . . . .  start main loop */
L50:
    if (mn7min_1.amin < 0.) {
	goto L95;
    }
    if (mn7flg_1.isw[1] <= 2) {
	goto L280;
    }
    ep = mn7min_1.amin * (float).1;
    if (sig2 < ep && mn7min_1.edm < ep) {
	goto L100;
    }
    sig2 = mn7min_1.edm;
    if (mn7cnv_1.nfcn - npfn > mn7cnv_1.nfcnmx) {
	goto L300;
    }
/*         calculate new point * by reflection */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pb = (float)0.;
	i__2 = nparp1;
	for (j = 1; j <= i__2; ++j) {
/* L59: */
	    pb += wg * mn7sim_1.p[i__ + j * 100 - 101];
	}
	mn7sim_1.pbar[i__ - 1] = pb - wg * mn7sim_1.p[i__ + jh * 100 - 101];
/* L60: */
	mn7sim_1.pstar[i__ - 1] = mn7sim_1.pbar[i__ - 1] * 2. - mn7sim_1.p[
		i__ + jh * 100 - 101] * 1.;
    }
    mncalf_((S_fp)fcn, mn7sim_1.pstar, &ycalf, (U_fp)futil);
    ystar = ycalf;
    if (ystar >= mn7min_1.amin) {
	goto L70;
    }
/*         point * better than jl, calculate new point ** */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L61: */
	mn7sim_1.pstst[i__ - 1] = mn7sim_1.pstar[i__ - 1] * 2. + 
		mn7sim_1.pbar[i__ - 1] * -1.;
    }
    mncalf_((S_fp)fcn, mn7sim_1.pstst, &ycalf, (U_fp)futil);
    ystst = ycalf;
/* L66: */
    if (ystst < y[jl - 1]) {
	goto L67;
    }
    mnrazz_(&ystar, mn7sim_1.pstar, y, &jh, &jl);
    goto L50;
L67:
    mnrazz_(&ystst, mn7sim_1.pstst, y, &jh, &jl);
    goto L50;
/*         point * is not as good as jl */
L70:
    if (ystar >= y[jh - 1]) {
	goto L73;
    }
    jhold = jh;
    mnrazz_(&ystar, mn7sim_1.pstar, y, &jh, &jl);
    if (jhold != jh) {
	goto L50;
    }
/*         calculate new point ** */
L73:
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L74: */
	mn7sim_1.pstst[i__ - 1] = mn7sim_1.p[i__ + jh * 100 - 101] * .5 + 
		mn7sim_1.pbar[i__ - 1] * .5;
    }
    mncalf_((S_fp)fcn, mn7sim_1.pstst, &ycalf, (U_fp)futil);
    ystst = ycalf;
    if (ystst > y[jh - 1]) {
	goto L30;
    }
/*     point ** is better than jh */
    if (ystst < mn7min_1.amin) {
	goto L67;
    }
    mnrazz_(&ystst, mn7sim_1.pstst, y, &jh, &jl);
    goto L50;
/*                                        . . . . . .  end main loop */
L95:
    if (mn7flg_1.isw[4] >= 0) {
	io___577.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___577);
	e_wsfe();
    }
    reg = (float).1;
/*                                        . . . . . ask if point is new */
L100:
    mninex_(mn7int_1.x);
    (*fcn)(&nparx, mn7der_1.gin, &mn7min_1.amin, mn7ext_1.u, &c__4, (U_fp)
	    futil);
    ++mn7cnv_1.nfcn;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7int_1.dirin[i__ - 1] = reg * dsav[i__ - 1];
	if ((d__1 = mn7int_1.x[i__ - 1] - mn7int_1.xt[i__ - 1], abs(d__1)) > 
		mn7int_1.dirin[i__ - 1]) {
	    goto L150;
	}
/* L120: */
    }
    goto L230;
L150:
    mn7cnv_1.nfcnmx = mn7cnv_1.nfcnmx + npfn - mn7cnv_1.nfcn;
    npfn = mn7cnv_1.nfcn;
    mnsimp_((S_fp)fcn, (U_fp)futil);
    if (mn7min_1.amin >= mn7min_1.apsi) {
	goto L325;
    }
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7int_1.dirin[i__ - 1] = dsav[i__ - 1] * (float).1;
	if ((d__1 = mn7int_1.x[i__ - 1] - mn7int_1.xt[i__ - 1], abs(d__1)) > 
		mn7int_1.dirin[i__ - 1]) {
	    goto L250;
	}
/* L220: */
    }
L230:
    if (mn7min_1.amin < mn7min_1.apsi) {
	goto L350;
    }
    goto L325;
/*                                        . . . . . . truly new minimum */
L250:
    mn7log_1.lnewmn = TRUE_;
    if (mn7flg_1.isw[1] >= 1) {
	mn7flg_1.isw[1] = 1;
	mn7min_1.dcovar = max(mn7min_1.dcovar,.5);
    } else {
	mn7min_1.dcovar = (float)1.;
    }
    mn7cnv_1.itaur = 0;
    mn7cnv_1.nfcnmx = mn7cnv_1.nfcnmx + npfn - mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "NEW MINIMU", (ftnlen)10, (ftnlen)10);
    if (mn7flg_1.isw[4] >= 0) {
	io___578.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___578);
	e_wsfe();
    }
    return 0;
/*                                        . . . return to previous region */
L280:
    if (mn7flg_1.isw[4] > 0) {
	io___579.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___579);
	e_wsfe();
    }
    goto L325;
L300:
    mn7flg_1.isw[0] = 1;
L325:
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7int_1.dirin[i__ - 1] = dsav[i__ - 1] * (float).01;
/* L330: */
	mn7int_1.x[i__ - 1] = mn7int_1.xt[i__ - 1];
    }
    mn7min_1.amin = mn7min_1.apsi;
    mn7min_1.edm = sigsav;
L350:
    mninex_(mn7int_1.x);
    if (mn7flg_1.isw[4] > 0) {
	io___580.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___580);
	e_wsfe();
    }
    s_copy(mn7tit_1.cstatu, "UNCHANGED ", (ftnlen)10, (ftnlen)10);
    mnrset_(&c__0);
    if (mn7flg_1.isw[1] < 2) {
	goto L380;
    }
    if (loop < nloop && mn7flg_1.isw[0] < 1) {
	goto L20;
    }
L380:
    mnprin_(&c__5, &mn7min_1.amin);
    mn7cnv_1.itaur = 0;
    return 0;
} /* mnimpr_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mninex_(doublereal *pint)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    static integer i__, j;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Transforms from internal coordinates (PINT) to external */
/* C        parameters (U).   The minimizing routines which work in */
/* C        internal coordinates call this routine before calling FCN. */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    /* Parameter adjustments */
    --pint;

    /* Function Body */
    i__1 = mn7npr_1.npar;
    for (j = 1; j <= i__1; ++j) {
	i__ = mn7inx_1.nexofi[j - 1];
	if (mn7inx_1.nvarl[i__ - 1] == 1) {
	    mn7ext_1.u[i__ - 1] = pint[j];
	} else {
	    mn7ext_1.u[i__ - 1] = mn7ext_1.alim[i__ - 1] + (sin(pint[j]) + (
		    float)1.) * (float).5 * (mn7ext_1.blim[i__ - 1] - 
		    mn7ext_1.alim[i__ - 1]);
	}
/* L100: */
    }
    return 0;
} /* mninex_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.4  1997/09/02 15:16:08  mclareni */
/* WINNT corrections */

/* Revision 1.3  1997/03/14 17:18:00  mclareni */
/* WNT mods */

/* Revision 1.2.2.1  1997/01/21 11:33:28  mclareni */
/* All mods for Winnt 96a on winnt branch */

/* Revision 1.2  1996/03/15 18:02:47  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mninit_(integer *i1, integer *i2, integer *i3)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double sqrt(doublereal), atan(doublereal);

    /* Local variables */
    static integer i__, idb;
    static doublereal piby2, epsp1, dummy, epsbak;
    extern logical intrac_(doublereal *);
    extern /* Subroutine */ int mncler_(), mntiny_(doublereal *, doublereal *)
	    ;
    static doublereal epstry, distnn;

    /* Fortran I/O blocks */
    static cilist io___589 = { 0, 0, 0, "(A,A,E10.2)", 0 };
    static cilist io___592 = { 0, 0, 0, "(3A,I3,A,I3,A,E10.2)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        This is the main initialization subroutine for MINUIT */
/* C     It initializes some constants in common */
/* C                (including the logical I/O unit nos.), */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



/*            I/O unit numbers */
    mn7iou_1.isysrd = *i1;
    mn7iou_1.isyswr = *i2;
    mn7io2_1.istkwr[0] = mn7iou_1.isyswr;
    mn7io2_1.nstkwr = 1;
    mn7iou_1.isyssa = *i3;
    mn7io2_1.nstkrd = 0;
/*               version identifier */
    s_copy(mn7tit_1.cvrsn, "96.03 ", (ftnlen)6, (ftnlen)6);
/*               some CONSTANT constants in COMMON */
    mn7npr_1.maxint = 100;
    mn7npr_1.maxext = 100;
    mn7cns_1.undefi = (float)-54321.;
    mn7cns_1.bigedm = (float)123456.;
    s_copy(mn7tit_1.cundef, ")UNDEFINED", (ftnlen)10, (ftnlen)10);
    s_copy(mn7tit_1.covmes, "NO ERROR MATRIX       ", (ftnlen)22, (ftnlen)22);
    s_copy(mn7tit_1.covmes + 22, "ERR MATRIX APPROXIMATE", (ftnlen)22, (
	    ftnlen)22);
    s_copy(mn7tit_1.covmes + 44, "ERR MATRIX NOT POS-DEF", (ftnlen)22, (
	    ftnlen)22);
    s_copy(mn7tit_1.covmes + 66, "ERROR MATRIX ACCURATE ", (ftnlen)22, (
	    ftnlen)22);
/*                some starting values in COMMON */
    mn7flg_1.nblock = 0;
    mn7flg_1.icomnd = 0;
    s_copy(mn7tit_1.ctitl, mn7tit_1.cundef, (ftnlen)50, (ftnlen)10);
    s_copy(mn7tit_1.cfrom, "INPUT   ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "INITIALIZE", (ftnlen)10, (ftnlen)10);
    mn7flg_1.isw[2] = 0;
    mn7flg_1.isw[3] = 0;
    mn7flg_1.isw[4] = 1;
/*         ISW(6)=0 for batch jobs,  =1 for interactive jobs */
/*                      =-1 for originally interactive temporarily batch */
    mn7flg_1.isw[5] = 0;
    if (intrac_(&dummy)) {
	mn7flg_1.isw[5] = 1;
    }
/*        DEBUG options set to default values */
    for (idb = 0; idb <= 10; ++idb) {
/* L10: */
	mn7flg_1.idbg[idb] = 0;
    }
    mn7log_1.lrepor = FALSE_;
    mn7log_1.lwarn = TRUE_;
    mn7log_1.limset = FALSE_;
    mn7log_1.lnewmn = FALSE_;
    mn7cnv_1.istrat = 1;
    mn7cnv_1.itaur = 0;
/*        default page dimensions and 'new page' carriage control integer */
    mn7iou_1.npagwd = 120;
    mn7iou_1.npagln = 56;
    mn7iou_1.newpag = 1;
    if (mn7flg_1.isw[5] > 0) {
	mn7iou_1.npagwd = 80;
	mn7iou_1.npagln = 30;
	mn7iou_1.newpag = 0;
    }
    mn7min_1.up = (float)1.;
    mn7cns_1.updflt = mn7min_1.up;
/*                   determine machine accuracy epsmac */
    epstry = (float).5;
    for (i__ = 1; i__ <= 100; ++i__) {
	epstry *= (float).5;
	epsp1 = epstry + 1.;
	mntiny_(&epsp1, &epsbak);
	if (epsbak < epstry) {
	    goto L35;
	}
/* L33: */
    }
    epstry = (float)1e-7;
    mn7cns_1.epsmac = epstry * (float)4.;
    io___589.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___589);
    do_fio(&c__1, " MNINIT UNABLE TO DETERMINE", (ftnlen)27);
    do_fio(&c__1, " ARITHMETIC PRECISION. WILL ASSUME:", (ftnlen)35);
    do_fio(&c__1, (char *)&mn7cns_1.epsmac, (ftnlen)sizeof(doublereal));
    e_wsfe();
L35:
    mn7cns_1.epsmac = epstry * (float)8.;
    mn7cns_1.epsma2 = sqrt(mn7cns_1.epsmac) * (float)2.;
/*                 the vlims are a non-negligible distance from pi/2 */
/*         used by MNPINT to set variables "near" the physical limits */
    piby2 = atan((float)1.) * (float)2.;
    distnn = sqrt(mn7cns_1.epsma2) * (float)8.;
    mn7cns_1.vlimhi = piby2 - distnn;
    mn7cns_1.vlimlo = -piby2 + distnn;
    mncler_();
    io___592.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___592);
    do_fio(&c__1, "  MINUIT RELEASE ", (ftnlen)17);
    do_fio(&c__1, mn7tit_1.cvrsn, (ftnlen)6);
    do_fio(&c__1, " INITIALIZED.   DIMENSIONS ", (ftnlen)27);
    do_fio(&c__1, (char *)&c__100, (ftnlen)sizeof(integer));
    do_fio(&c__1, "/", (ftnlen)1);
    do_fio(&c__1, (char *)&c__100, (ftnlen)sizeof(integer));
    do_fio(&c__1, "  EPSMAC=", (ftnlen)9);
    do_fio(&c__1, (char *)&mn7cns_1.epsmac, (ftnlen)sizeof(doublereal));
    e_wsfe();
    return 0;
} /* mninit_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mninpu_(integer *iunit, integer *ierr)
{
    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Fortran I/O blocks */
    static cilist io___593 = { 0, 0, 0, "(A)", 0 };
    static cilist io___594 = { 0, 0, 0, "(A)", 0 };
    static cilist io___595 = { 0, 0, 0, "(A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C      called by the user to SET INPUT to IUNIT, */
/* C      an alternative to MNSTIN where the user can specify just */
/* C      a logical unit number and he is not interrogated about */
/* C      open files and rewinding, all that is the responsibility */
/* C      of the user and cannot be fixed interactively. */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



    *ierr = 0;
/*                              IUNIT = 0, revert to previous input file */
    if (*iunit == 0) {
	if (mn7io2_1.nstkrd == 0) {
	    io___593.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___593);
	    do_fio(&c__1, " CALL TO MNINPU(0) IGNORED", (ftnlen)26);
	    e_wsfe();
	    io___594.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___594);
	    do_fio(&c__1, " ALREADY READING FROM PRIMARY INPUT", (ftnlen)35);
	    e_wsfe();
	} else {
	    mn7iou_1.isysrd = mn7io2_1.istkrd[mn7io2_1.nstkrd - 1];
	    --mn7io2_1.nstkrd;
	}

/*                               new input file */
    } else {
	if (mn7io2_1.nstkrd >= 10) {
	    io___595.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___595);
	    do_fio(&c__1, " INPUT FILE STACK SIZE EXCEEDED.", (ftnlen)32);
	    e_wsfe();
	    goto L800;
	}
	++mn7io2_1.nstkrd;
	mn7io2_1.istkrd[mn7io2_1.nstkrd - 1] = mn7iou_1.isysrd;
	mn7iou_1.isysrd = *iunit;
    }

    return 0;
L800:
    *ierr = 1;
    return 0;
} /* mninpu_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnintr_(S_fp fcn, U_fp futil)
{
    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    extern /* Subroutine */ int mnread_(S_fp, integer *, integer *, U_fp);
    static integer iflgin, iflgut;

    /* Fortran I/O blocks */
    static cilist io___598 = { 0, 0, 0, "(2A/)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       Called by user. Interfaces to MNREAD to allow user to change */
/* C       easily from Fortran-callable to interactive mode. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    iflgin = 3;
    mnread_((S_fp)fcn, &iflgin, &iflgut, (U_fp)futil);
    io___598.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___598);
    do_fio(&c__1, " END OF MINUIT COMMAND INPUT. ", (ftnlen)30);
    do_fio(&c__1, "   RETURN TO USER PROGRAM.", (ftnlen)26);
    e_wsfe();
    return 0;
} /* mnintr_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:48  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnlims_()
{
    /* Format strings */
    static char fmt_134[] = "(\002 LIMITS REMOVED FROM PARAMETER\002,i4)";
    static char fmt_237[] = "(\002 PARAMETER\002,i3,\002 LIMITS SET TO\002,2\
g15.5)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i2, ifx, inu;
    static doublereal dxdi;
    static integer kint;
    static doublereal snew;
    static integer newcod;
    extern /* Subroutine */ int mndxdi_(doublereal *, integer *, doublereal *)
	    , mnexin_(doublereal *), mnrset_(integer *);

    /* Fortran I/O blocks */
    static cilist io___603 = { 0, 0, 0, "(11X,A,I3)", 0 };
    static cilist io___604 = { 0, 0, 0, fmt_134, 0 };
    static cilist io___607 = { 0, 0, 0, fmt_237, 0 };
    static cilist io___608 = { 0, 0, 0, "(A,I3,A)", 0 };
    static cilist io___609 = { 0, 0, 0, "(A,I3)", 0 };
    static cilist io___611 = { 0, 0, 0, "(A)", 0 };
    static cilist io___612 = { 0, 0, 0, fmt_134, 0 };
    static cilist io___613 = { 0, 0, 0, "(A,I3)", 0 };
    static cilist io___614 = { 0, 0, 0, fmt_237, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       Called from MNSET */
/* C       Interprets the SET LIM command, to reset the parameter limits */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



    s_copy(mn7tit_1.cfrom, "SET LIM ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "NO CHANGE ", (ftnlen)10, (ftnlen)10);
    i2 = (integer) mn7arg_1.word7[0];
    if (i2 > mn7npr_1.maxext || i2 < 0) {
	goto L900;
    }
    if (i2 > 0) {
	goto L30;
    }
/*                                     set limits on all parameters */
    newcod = 4;
    if (mn7arg_1.word7[1] == mn7arg_1.word7[2]) {
	newcod = 1;
    }
    i__1 = mn7npr_1.nu;
    for (inu = 1; inu <= i__1; ++inu) {
	if (mn7inx_1.nvarl[inu - 1] <= 0) {
	    goto L20;
	}
	if (mn7inx_1.nvarl[inu - 1] == 1 && newcod == 1) {
	    goto L20;
	}
	kint = mn7inx_1.niofex[inu - 1];
/*             see if parameter has been fixed */
	if (kint <= 0) {
	    if (mn7flg_1.isw[4] >= 0) {
		io___603.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___603);
		do_fio(&c__1, " LIMITS NOT CHANGED FOR FIXED PARAMETER:", (
			ftnlen)40);
		do_fio(&c__1, (char *)&inu, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    goto L20;
	}
	if (newcod == 1) {
/*            remove limits from parameter */
	    if (mn7flg_1.isw[4] > 0) {
		io___604.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___604);
		do_fio(&c__1, (char *)&inu, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    s_copy(mn7tit_1.cstatu, "NEW LIMITS", (ftnlen)10, (ftnlen)10);
	    mndxdi_(&mn7int_1.x[kint - 1], &kint, &dxdi);
	    snew = mn7der_1.gstep[kint - 1] * dxdi;
	    mn7der_1.gstep[kint - 1] = abs(snew);
	    mn7inx_1.nvarl[inu - 1] = 1;
	} else {
/*             put limits on parameter */
	    mn7ext_1.alim[inu - 1] = min(mn7arg_1.word7[1],mn7arg_1.word7[2]);
	    mn7ext_1.blim[inu - 1] = max(mn7arg_1.word7[1],mn7arg_1.word7[2]);
	    if (mn7flg_1.isw[4] > 0) {
		io___607.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___607);
		do_fio(&c__1, (char *)&inu, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&mn7ext_1.alim[inu - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&mn7ext_1.blim[inu - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    mn7inx_1.nvarl[inu - 1] = 4;
	    s_copy(mn7tit_1.cstatu, "NEW LIMITS", (ftnlen)10, (ftnlen)10);
	    mn7der_1.gstep[kint - 1] = (float)-.1;
	}
L20:
	;
    }
    goto L900;
/*                                       set limits on one parameter */
L30:
    if (mn7inx_1.nvarl[i2 - 1] <= 0) {
	io___608.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___608);
	do_fio(&c__1, " PARAMETER ", (ftnlen)11);
	do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	do_fio(&c__1, " IS NOT VARIABLE.", (ftnlen)17);
	e_wsfe();
	goto L900;
    }
    kint = mn7inx_1.niofex[i2 - 1];
/*                                       see if parameter was fixed */
    if (kint == 0) {
	io___609.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___609);
	do_fio(&c__1, " REQUEST TO CHANGE LIMITS ON FIXED PARAMETER:", (
		ftnlen)45);
	do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	e_wsfe();
	i__1 = mn7fx1_1.npfix;
	for (ifx = 1; ifx <= i__1; ++ifx) {
	    if (i2 == mn7fx1_1.ipfix[ifx - 1]) {
		goto L92;
	    }
/* L82: */
	}
	io___611.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___611);
	do_fio(&c__1, " MINUIT BUG IN MNLIMS. SEE F. JAMES", (ftnlen)35);
	e_wsfe();
L92:
	;
    }
    if (mn7arg_1.word7[1] != mn7arg_1.word7[2]) {
	goto L235;
    }
/*                                       remove limits */
    if (mn7inx_1.nvarl[i2 - 1] != 1) {
	if (mn7flg_1.isw[4] > 0) {
	    io___612.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___612);
	    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	s_copy(mn7tit_1.cstatu, "NEW LIMITS", (ftnlen)10, (ftnlen)10);
	if (kint <= 0) {
	    mn7fx3_1.gsteps[ifx - 1] = (d__1 = mn7fx3_1.gsteps[ifx - 1], abs(
		    d__1));
	} else {
	    mndxdi_(&mn7int_1.x[kint - 1], &kint, &dxdi);
	    if (abs(dxdi) < (float).01) {
		dxdi = (float).01;
	    }
	    mn7der_1.gstep[kint - 1] = (d__1 = mn7der_1.gstep[kint - 1] * 
		    dxdi, abs(d__1));
	    mn7der_1.grd[kint - 1] *= dxdi;
	}
	mn7inx_1.nvarl[i2 - 1] = 1;
    } else {
	io___613.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___613);
	do_fio(&c__1, " NO LIMITS SPECIFIED.  PARAMETER ", (ftnlen)33);
	do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	do_fio(&c__1, " IS ALREADY UNLIMITED.  NO CHANGE.", (ftnlen)34);
	e_wsfe();
    }
    goto L900;
/*                                        put on limits */
L235:
    mn7ext_1.alim[i2 - 1] = min(mn7arg_1.word7[1],mn7arg_1.word7[2]);
    mn7ext_1.blim[i2 - 1] = max(mn7arg_1.word7[1],mn7arg_1.word7[2]);
    mn7inx_1.nvarl[i2 - 1] = 4;
    if (mn7flg_1.isw[4] > 0) {
	io___614.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___614);
	do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&mn7ext_1.alim[i2 - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&mn7ext_1.blim[i2 - 1], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
    }
    s_copy(mn7tit_1.cstatu, "NEW LIMITS", (ftnlen)10, (ftnlen)10);
    if (kint <= 0) {
	mn7fx3_1.gsteps[ifx - 1] = (float)-.1;
    } else {
	mn7der_1.gstep[kint - 1] = (float)-.1;
    }

L900:
    if (s_cmp(mn7tit_1.cstatu, "NO CHANGE ", (ftnlen)10, (ftnlen)10) != 0) {
	mnexin_(mn7int_1.x);
	mnrset_(&c__1);
    }
    return 0;
} /* mnlims_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnline_(S_fp fcn, doublereal *start, doublereal *fstart, 
	doublereal *step, doublereal *slope, doublereal *toler, U_fp futil)
{
    /* Initialized data */

    static char charal[26+1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal f1, f2, f3;
    static integer kk, ipt;
    static doublereal xpq[12], ypq[12];
    static char chpq[1*12];
    static doublereal slam, sdev, coeff[3], denom, flast;
    static char cmess[60];
    static doublereal fvals[3], xvals[3];
    static integer nparx;
    static doublereal fvmin, xvmin, ratio, fvmax;
    static integer nvmax, nxypt;
    static doublereal toler8, toler9;
    static logical ldebug;
    static doublereal slamin, undral, overal;
    extern /* Subroutine */ int mninex_(doublereal *);
    static doublereal slamax;
    extern /* Subroutine */ int mnpfit_(doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *);
    static doublereal slopem;
    extern /* Subroutine */ int mnwarn_(char *, char *, char *, ftnlen, 
	    ftnlen, ftnlen), mnplot_(doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___621 = { 0, 0, 0, "(A/2E14.5/2X,10F10.5)", 0 };
    static cilist io___649 = { 0, 0, 0, "(A/(2X,6G12.4))", 0 };
    static cilist io___650 = { 0, 0, 0, "(' AFTER',I3,' POINTS,',A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Perform a line search from position START */
/* C        along direction STEP, where the length of vector STEP */
/* C                   gives the expected position of minimum. */
/* C        FSTART is value of function at START */
/* C        SLOPE (if non-zero) is df/dx along STEP at START */
/* C        TOLER is initial tolerance of minimum in direction STEP */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/* SLAMBG and ALPHA control the maximum individual steps allowed. */
/* The first step is always =1. The max length of second step is SLAMBG. */
/* The max size of subsequent steps is the maximum previous successful */
/*   step multiplied by ALPHA + the size of most recent successful step, */
/*   but cannot be smaller than SLAMBG. */
    /* Parameter adjustments */
    --step;
    --start;

    /* Function Body */
    ldebug = mn7flg_1.idbg[1] >= 1;
/*                  starting values for overall limits on total step SLAM */
    overal = (float)1e3;
    undral = (float)-100.;
/*                              debug check if start is ok */
    if (ldebug) {
	mninex_(&start[1]);
	(*fcn)(&nparx, mn7der_1.gin, &f1, mn7ext_1.u, &c__4, (U_fp)futil);
	++mn7cnv_1.nfcn;
	if (f1 != *fstart) {
	    io___621.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___621);
	    do_fio(&c__1, " MNLINE start point not consistent, F values, par\
ameters=", (ftnlen)57);
	    i__1 = mn7npr_1.npar;
	    for (kk = 1; kk <= i__1; ++kk) {
		do_fio(&c__1, (char *)&mn7int_1.x[kk - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	}
    }
/*                                      . set up linear search along STEP */
    fvmin = *fstart;
    xvmin = 0.;
    nxypt = 1;
    *(unsigned char *)&chpq[0] = *(unsigned char *)&charal[0];
    xpq[0] = (float)0.;
    ypq[0] = *fstart;
/*               SLAMIN = smallest possible value of ABS(SLAM) */
    slamin = 0.;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (step[i__] == 0.) {
	    goto L20;
	}
	ratio = (d__1 = start[i__] / step[i__], abs(d__1));
	if (slamin == 0.) {
	    slamin = ratio;
	}
	if (ratio < slamin) {
	    slamin = ratio;
	}
L20:
	mn7int_1.x[i__ - 1] = start[i__] + step[i__];
    }
    if (slamin == 0.) {
	slamin = mn7cns_1.epsmac;
    }
    slamin *= mn7cns_1.epsma2;
    nparx = mn7npr_1.npar;

    mninex_(mn7int_1.x);
    (*fcn)(&nparx, mn7der_1.gin, &f1, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    ++nxypt;
    *(unsigned char *)&chpq[nxypt - 1] = *(unsigned char *)&charal[nxypt - 1];
    xpq[nxypt - 1] = (float)1.;
    ypq[nxypt - 1] = f1;
    if (f1 < *fstart) {
	fvmin = f1;
	xvmin = (float)1.;
    }
/*                         . quadr interp using slope GDEL and two points */
    slam = (float)1.;
    toler8 = *toler;
    slamax = 5.;
    flast = f1;
/*                         can iterate on two-points (cut) if no imprvmnt */
L25:
/* Computing 2nd power */
    d__1 = slam;
    denom = (flast - *fstart - *slope * slam) * (float)2. / (d__1 * d__1);
/*     IF (DENOM .EQ. ZERO)  DENOM = -0.1*SLOPE */
    slam = (float)1.;
    if (denom != 0.) {
	slam = -(*slope) / denom;
    }
    if (slam < 0.) {
	slam = slamax;
    }
    if (slam > slamax) {
	slam = slamax;
    }
    if (slam < toler8) {
	slam = toler8;
    }
    if (slam < slamin) {
	goto L80;
    }
    if ((d__1 = slam - (float)1., abs(d__1)) < toler8 && f1 < *fstart) {
	goto L70;
    }
    if ((d__1 = slam - (float)1., abs(d__1)) < toler8) {
	slam = toler8 + (float)1.;
    }
    if (nxypt >= 12) {
	goto L65;
    }
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	mn7int_1.x[i__ - 1] = start[i__] + slam * step[i__];
    }
    mninex_(mn7int_1.x);
    (*fcn)(&mn7npr_1.npar, mn7der_1.gin, &f2, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    ++nxypt;
    *(unsigned char *)&chpq[nxypt - 1] = *(unsigned char *)&charal[nxypt - 1];
    xpq[nxypt - 1] = slam;
    ypq[nxypt - 1] = f2;
    if (f2 < fvmin) {
	fvmin = f2;
	xvmin = slam;
    }
    if (*fstart == fvmin) {
	flast = f2;
	toler8 = *toler * slam;
	overal = slam - toler8;
	slamax = overal;
	goto L25;
    }
/*                                        . quadr interp using 3 points */
    xvals[0] = xpq[0];
    fvals[0] = ypq[0];
    xvals[1] = xpq[nxypt - 2];
    fvals[1] = ypq[nxypt - 2];
    xvals[2] = xpq[nxypt - 1];
    fvals[2] = ypq[nxypt - 1];
/*                             begin iteration, calculate desired step */
L50:
/* Computing MAX */
    d__1 = slamax, d__2 = abs(xvmin) * 2.;
    slamax = max(d__1,d__2);
    mnpfit_(xvals, fvals, &c__3, coeff, &sdev);
    if (coeff[2] <= 0.) {
	slopem = coeff[2] * (float)2. * xvmin + coeff[1];
	if (slopem <= 0.) {
	    slam = xvmin + slamax;
	} else {
	    slam = xvmin - slamax;
	}
    } else {
	slam = -coeff[1] / (coeff[2] * (float)2.);
	if (slam > xvmin + slamax) {
	    slam = xvmin + slamax;
	}
	if (slam < xvmin - slamax) {
	    slam = xvmin - slamax;
	}
    }
    if (slam > 0.) {
	if (slam > overal) {
	    slam = overal;
	}
    } else {
	if (slam < undral) {
	    slam = undral;
	}
    }
/*               come here if step was cut below */
L52:
/* Computing MAX */
    d__2 = toler8, d__3 = (d__1 = toler8 * slam, abs(d__1));
    toler9 = max(d__2,d__3);
    for (ipt = 1; ipt <= 3; ++ipt) {
	if ((d__1 = slam - xvals[ipt - 1], abs(d__1)) < toler9) {
	    goto L70;
	}
/* L55: */
    }
/*                take the step */
    if (nxypt >= 12) {
	goto L65;
    }
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	mn7int_1.x[i__ - 1] = start[i__] + slam * step[i__];
    }
    mninex_(mn7int_1.x);
    (*fcn)(&nparx, mn7der_1.gin, &f3, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    ++nxypt;
    *(unsigned char *)&chpq[nxypt - 1] = *(unsigned char *)&charal[nxypt - 1];
    xpq[nxypt - 1] = slam;
    ypq[nxypt - 1] = f3;
/*             find worst previous point out of three */
    fvmax = fvals[0];
    nvmax = 1;
    if (fvals[1] > fvmax) {
	fvmax = fvals[1];
	nvmax = 2;
    }
    if (fvals[2] > fvmax) {
	fvmax = fvals[2];
	nvmax = 3;
    }
/*              if latest point worse than all three previous, cut step */
    if (f3 >= fvmax) {
	if (nxypt >= 12) {
	    goto L65;
	}
	if (slam > xvmin) {
/* Computing MIN */
	    d__1 = overal, d__2 = slam - toler8;
	    overal = min(d__1,d__2);
	}
	if (slam < xvmin) {
/* Computing MAX */
	    d__1 = undral, d__2 = slam + toler8;
	    undral = max(d__1,d__2);
	}
	slam = (slam + xvmin) * (float).5;
	goto L52;
    }
/*              prepare another iteration, replace worst previous point */
    xvals[nvmax - 1] = slam;
    fvals[nvmax - 1] = f3;
    if (f3 < fvmin) {
	fvmin = f3;
	xvmin = slam;
    } else {
	if (slam > xvmin) {
/* Computing MIN */
	    d__1 = overal, d__2 = slam - toler8;
	    overal = min(d__1,d__2);
	}
	if (slam < xvmin) {
/* Computing MAX */
	    d__1 = undral, d__2 = slam + toler8;
	    undral = max(d__1,d__2);
	}
    }
    if (nxypt < 12) {
	goto L50;
    }
/*                                            . . end of iteration . . . */
/*            stop because too many iterations */
L65:
    s_copy(cmess, " LINE SEARCH HAS EXHAUSTED THE LIMIT OF FUNCTION CALLS ", (
	    ftnlen)60, (ftnlen)55);
    if (ldebug) {
	io___649.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___649);
	do_fio(&c__1, " MNLINE DEBUG: steps=", (ftnlen)21);
	i__1 = mn7npr_1.npar;
	for (kk = 1; kk <= i__1; ++kk) {
	    do_fio(&c__1, (char *)&step[kk], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    goto L100;
/*            stop because within tolerance */
L70:
    s_copy(cmess, " LINE SEARCH HAS ATTAINED TOLERANCE ", (ftnlen)60, (ftnlen)
	    36);
    goto L100;
L80:
    s_copy(cmess, " STEP SIZE AT ARITHMETICALLY ALLOWED MINIMUM", (ftnlen)60, 
	    (ftnlen)44);
L100:
    mn7min_1.amin = fvmin;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7int_1.dirin[i__ - 1] = step[i__] * xvmin;
/* L120: */
	mn7int_1.x[i__ - 1] = start[i__] + mn7int_1.dirin[i__ - 1];
    }
    mninex_(mn7int_1.x);
    if (xvmin < (float)0.) {
	mnwarn_("D", "MNLINE", " LINE MINIMUM IN BACKWARDS DIRECTION", (
		ftnlen)1, (ftnlen)6, (ftnlen)36);
    }
    if (fvmin == *fstart) {
	mnwarn_("D", "MNLINE", " LINE SEARCH FINDS NO IMPROVEMENT ", (ftnlen)
		1, (ftnlen)6, (ftnlen)34);
    }
    if (ldebug) {
	io___650.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___650);
	do_fio(&c__1, (char *)&nxypt, (ftnlen)sizeof(integer));
	do_fio(&c__1, cmess, (ftnlen)60);
	e_wsfe();
	mnplot_(xpq, ypq, chpq, &nxypt, &mn7iou_1.isyswr, &mn7iou_1.npagwd, &
		mn7iou_1.npagln, (ftnlen)1);
    }
    return 0;
} /* mnline_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnmatu_(integer *kode)
{
    /* Format strings */
    static char fmt_150[] = "(/\002 PARAMETER  CORRELATION COEFFICIENTS\002\
/\002       NO.  GLOBAL\002,20i6)";
    static char fmt_171[] = "(6x,i3,2x,f7.5,1x,20f6.3)";
    static char fmt_181[] = "(19x,20f6.3)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, m, n, id, it, ix, ndi, ndj, iso, isw2, isw5, ndex, 
	    ncoef;
    static doublereal vline[100];
    static integer nparm;
    extern /* Subroutine */ int mnemat_(doublereal *, integer *), mnwerr_();
    static integer nsofar;

    /* Fortran I/O blocks */
    static cilist io___652 = { 0, 0, 0, "(1X,A)", 0 };
    static cilist io___653 = { 0, 0, 0, "(' MNMATU: NPAR=0')", 0 };
    static cilist io___655 = { 0, 0, 0, "(1X,A)", 0 };
    static cilist io___658 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___669 = { 0, 0, 0, fmt_171, 0 };
    static cilist io___673 = { 0, 0, 0, fmt_181, 0 };
    static cilist io___674 = { 0, 0, 0, "(1X,A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        prints the covariance matrix v when KODE=1. */
/* C        always prints the global correlations, and */
/* C        calculates and prints the individual correlation coefficients */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    isw2 = mn7flg_1.isw[1];
    if (isw2 < 1) {
	io___652.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___652);
	do_fio(&c__1, mn7tit_1.covmes + isw2 * 22, (ftnlen)22);
	e_wsfe();
	goto L500;
    }
    if (mn7npr_1.npar == 0) {
	io___653.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___653);
	e_wsfe();
	goto L500;
    }
/*                                       . . . . .external error matrix */
    if (*kode == 1) {
	isw5 = mn7flg_1.isw[4];
	mn7flg_1.isw[4] = 2;
	mnemat_(mn7sim_1.p, &mn7npr_1.maxint);
	if (isw2 < 3) {
	    io___655.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___655);
	    do_fio(&c__1, mn7tit_1.covmes + isw2 * 22, (ftnlen)22);
	    e_wsfe();
	}
	mn7flg_1.isw[4] = isw5;
    }
/*                                       . . . . . correlation coeffs. . */
    if (mn7npr_1.npar <= 1) {
	goto L500;
    }
    mnwerr_();
/*     NCOEF is number of coeff. that fit on one line, not to exceed 20 */
    ncoef = (mn7iou_1.npagwd - 19) / 6;
    ncoef = min(ncoef,20);
    nparm = min(mn7npr_1.npar,ncoef);
    io___658.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___658);
    i__1 = nparm;
    for (id = 1; id <= i__1; ++id) {
	do_fio(&c__1, (char *)&mn7inx_1.nexofi[id - 1], (ftnlen)sizeof(
		integer));
    }
    e_wsfe();
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ix = mn7inx_1.nexofi[i__ - 1];
	ndi = i__ * (i__ + 1) / 2;
	i__2 = mn7npr_1.npar;
	for (j = 1; j <= i__2; ++j) {
	    m = max(i__,j);
	    n = min(i__,j);
	    ndex = m * (m - 1) / 2 + n;
	    ndj = j * (j + 1) / 2;
/* L170: */
	    vline[j - 1] = mn7var_1.vhmat[ndex - 1] / sqrt((d__1 = 
		    mn7var_1.vhmat[ndi - 1] * mn7var_1.vhmat[ndj - 1], abs(
		    d__1)));
	}
	nparm = min(mn7npr_1.npar,ncoef);
	io___669.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___669);
	do_fio(&c__1, (char *)&ix, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&mn7err_1.globcc[i__ - 1], (ftnlen)sizeof(
		doublereal));
	i__2 = nparm;
	for (it = 1; it <= i__2; ++it) {
	    do_fio(&c__1, (char *)&vline[it - 1], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	if (i__ <= nparm) {
	    goto L200;
	}
	for (iso = 1; iso <= 10; ++iso) {
	    nsofar = nparm;
/* Computing MIN */
	    i__2 = mn7npr_1.npar, i__3 = nsofar + ncoef;
	    nparm = min(i__2,i__3);
	    io___673.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___673);
	    i__2 = nparm;
	    for (it = nsofar + 1; it <= i__2; ++it) {
		do_fio(&c__1, (char *)&vline[it - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    if (i__ <= nparm) {
		goto L192;
	    }
/* L190: */
	}
L192:
L200:
	;
    }
    if (isw2 < 3) {
	io___674.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___674);
	do_fio(&c__1, mn7tit_1.covmes + isw2 * 22, (ftnlen)22);
	e_wsfe();
    }
L500:
    return 0;
} /* mnmatu_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:49  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnmigr_(S_fp fcn, U_fp futil)
{
    /* Format strings */
    static char fmt_470[] = "(\002 START MIGRAD MINIMIZATION.  STRATEGY\002,\
i2,\002.  CONVERGENCE WHEN EDM .LT.\002,e9.2)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, m, n, kk;
    static doublereal fs, gs[100], ri, vg[100], gvg, vgi, xxs[100], gdel, 
	    gami;
    static integer npfn;
    static doublereal flnu[100];
    static integer lenv, ndex, iext, iter;
    static doublereal step[100], dsum, gssq, vsum;
    static integer npsdf, nparx;
    static doublereal fzero;
    static integer iswtr, lined2;
    static doublereal delgam;
    static logical ldebug;
    static integer nfcnmg;
    extern /* Subroutine */ int mnamin_(S_fp, U_fp), mnderi_(S_fp, U_fp), 
	    mnline_(S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, U_fp);
    static doublereal rhotol;
    extern /* Subroutine */ int mnhess_(S_fp, U_fp), mnwerr_(), mninex_(
	    doublereal *), mnwarn_(char *, char *, char *, ftnlen, ftnlen, 
	    ftnlen), mnprin_(integer *, doublereal *);
    static integer nrstrt;
    extern /* Subroutine */ int mnmatu_(integer *), mnpsdf_();
    static doublereal gdgssq;

    /* Fortran I/O blocks */
    static cilist io___685 = { 0, 0, 0, fmt_470, 0 };
    static cilist io___693 = { 0, 0, 0, "(A,I3,2G13.3)", 0 };
    static cilist io___696 = { 0, 0, 0, "(A,A/(1X,10G10.2))", 0 };
    static cilist io___710 = { 0, 0, 0, "(A,(1X,10G10.3))", 0 };
    static cilist io___714 = { 0, 0, 0, "(A,F5.1,A)", 0 };
    static cilist io___715 = { 0, 0, 0, "(A,(1X,10G10.3))", 0 };
    static cilist io___717 = { 0, 0, 0, "(A)", 0 };
    static cilist io___718 = { 0, 0, 0, "(A)", 0 };
    static cilist io___719 = { 0, 0, 0, "(A)", 0 };
    static cilist io___720 = { 0, 0, 0, "(A)", 0 };
    static cilist io___721 = { 0, 0, 0, "(A)", 0 };
    static cilist io___722 = { 0, 0, 0, "(/A)", 0 };
    static cilist io___723 = { 0, 0, 0, "(/A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Performs a local function minimization using basically the */
/* C        method of Davidon-Fletcher-Powell as modified by Fletcher */
/* C        ref. -- Fletcher, Comp.J. 13,317 (1970)   "switching method" */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    if (mn7npr_1.npar <= 0) {
	return 0;
    }
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    ldebug = mn7flg_1.idbg[4] >= 1;
    s_copy(mn7tit_1.cfrom, "MIGRAD  ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    nfcnmg = mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "INITIATE  ", (ftnlen)10, (ftnlen)10);
    iswtr = mn7flg_1.isw[4] - (mn7cnv_1.itaur << 1);
    npfn = mn7cnv_1.nfcn;
    nparx = mn7npr_1.npar;
    lenv = mn7npr_1.npar * (mn7npr_1.npar + 1) / 2;
    nrstrt = 0;
    npsdf = 0;
    lined2 = 0;
    mn7flg_1.isw[3] = -1;
    rhotol = mn7min_1.apsi * (float).001;
    if (iswtr >= 1) {
	io___685.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___685);
	do_fio(&c__1, (char *)&mn7cnv_1.istrat, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rhotol, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*                                           initialization strategy */
    if (mn7cnv_1.istrat < 2 || mn7flg_1.isw[1] >= 3) {
	goto L2;
    }
/*                                come (back) here to restart completely */
L1:
    if (nrstrt > mn7cnv_1.istrat) {
	s_copy(mn7tit_1.cstatu, "FAILED    ", (ftnlen)10, (ftnlen)10);
	mn7flg_1.isw[3] = -1;
	goto L230;
    }
/*                                      . get full covariance and gradient */
    mnhess_((S_fp)fcn, (U_fp)futil);
    mnwerr_();
    npsdf = 0;
    if (mn7flg_1.isw[1] >= 1) {
	goto L10;
    }
/*                                        . get gradient at start point */
L2:
    mninex_(mn7int_1.x);
    if (mn7flg_1.isw[2] == 1) {
	(*fcn)(&nparx, mn7der_1.gin, &fzero, mn7ext_1.u, &c__2, (U_fp)futil);
	++mn7cnv_1.nfcn;
    }
    mnderi_((S_fp)fcn, (U_fp)futil);
    if (mn7flg_1.isw[1] >= 1) {
	goto L10;
    }
/*                                   sometimes start with diagonal matrix */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xxs[i__ - 1] = mn7int_1.x[i__ - 1];
	step[i__ - 1] = 0.;
/* L3: */
    }
/*                           do line search if second derivative negative */
    ++lined2;
    if (lined2 < (mn7cnv_1.istrat + 1) * mn7npr_1.npar) {
	i__1 = mn7npr_1.npar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (mn7der_1.g2[i__ - 1] > 0.) {
		goto L5;
	    }
	    step[i__ - 1] = -d_sign(&mn7der_1.gstep[i__ - 1], &mn7der_1.grd[
		    i__ - 1]);
	    gdel = step[i__ - 1] * mn7der_1.grd[i__ - 1];
	    fs = mn7min_1.amin;
	    mnline_((S_fp)fcn, xxs, &fs, step, &gdel, &c_b1209, (U_fp)futil);
	    mnwarn_("D", "MNMIGR", "Negative G2 line search", (ftnlen)1, (
		    ftnlen)6, (ftnlen)23);
	    iext = mn7inx_1.nexofi[i__ - 1];
	    if (ldebug) {
		io___693.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___693);
		do_fio(&c__1, " Negative G2 line search, param ", (ftnlen)32);
		do_fio(&c__1, (char *)&iext, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fs, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&mn7min_1.amin, (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    goto L2;
L5:
	    ;
	}
    }
/*                           make diagonal error matrix */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ndex = i__ * (i__ - 1) / 2;
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    ++ndex;
/* L7: */
	    mn7var_1.vhmat[ndex - 1] = (float)0.;
	}
	++ndex;
	if (mn7der_1.g2[i__ - 1] <= 0.) {
	    mn7der_1.g2[i__ - 1] = (float)1.;
	}
	mn7var_1.vhmat[ndex - 1] = (float)2. / mn7der_1.g2[i__ - 1];
/* L8: */
    }
    mn7min_1.dcovar = (float)1.;
    if (ldebug) {
	io___696.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___696);
	do_fio(&c__1, " DEBUG MNMIGR,", (ftnlen)14);
	do_fio(&c__1, " STARTING MATRIX DIAGONAL,  VHMAT=", (ftnlen)34);
	i__1 = lenv;
	for (kk = 1; kk <= i__1; ++kk) {
	    do_fio(&c__1, (char *)&mn7var_1.vhmat[kk - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
    }
/*                                         ready to start first iteration */
L10:
    ++nrstrt;
    if (nrstrt > mn7cnv_1.istrat + 1) {
	s_copy(mn7tit_1.cstatu, "FAILED    ", (ftnlen)10, (ftnlen)10);
	goto L230;
    }
    fs = mn7min_1.amin;
/*                                        . . . get EDM and set up loop */
    mn7min_1.edm = (float)0.;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gs[i__ - 1] = mn7der_1.grd[i__ - 1];
	xxs[i__ - 1] = mn7int_1.x[i__ - 1];
	ndex = i__ * (i__ - 1) / 2;
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    ++ndex;
/* L17: */
	    mn7min_1.edm += gs[i__ - 1] * mn7var_1.vhmat[ndex - 1] * gs[j - 1]
		    ;
	}
	++ndex;
/* L18: */
/* Computing 2nd power */
	d__1 = gs[i__ - 1];
	mn7min_1.edm += d__1 * d__1 * (float).5 * mn7var_1.vhmat[ndex - 1];
    }
    mn7min_1.edm = mn7min_1.edm * (float).5 * (mn7min_1.dcovar * (float)3. + (
	    float)1.);
    if (mn7min_1.edm < 0.) {
	mnwarn_("W", "MIGRAD", "STARTING MATRIX NOT POS-DEFINITE.", (ftnlen)1,
		 (ftnlen)6, (ftnlen)33);
	mn7flg_1.isw[1] = 0;
	mn7min_1.dcovar = (float)1.;
	goto L2;
    }
    if (mn7flg_1.isw[1] == 0) {
	mn7min_1.edm = mn7cns_1.bigedm;
    }
    iter = 0;
    mninex_(mn7int_1.x);
    mnwerr_();
    if (iswtr >= 1) {
	mnprin_(&c__3, &mn7min_1.amin);
    }
    if (iswtr >= 2) {
	mnmatu_(&c__0);
    }
/*                                        . . . . .  start main loop */
L24:
    if (mn7cnv_1.nfcn - npfn >= mn7cnv_1.nfcnmx) {
	goto L190;
    }
    gdel = (float)0.;
    gssq = (float)0.;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ri = (float)0.;
/* Computing 2nd power */
	d__1 = gs[i__ - 1];
	gssq += d__1 * d__1;
	i__2 = mn7npr_1.npar;
	for (j = 1; j <= i__2; ++j) {
	    m = max(i__,j);
	    n = min(i__,j);
	    ndex = m * (m - 1) / 2 + n;
/* L25: */
	    ri += mn7var_1.vhmat[ndex - 1] * gs[j - 1];
	}
	step[i__ - 1] = ri * (float)-.5;
/* L30: */
	gdel += step[i__ - 1] * gs[i__ - 1];
    }
    if (gssq == 0.) {
	mnwarn_("D", "MIGRAD", " FIRST DERIVATIVES OF FCN ARE ALL ZERO", (
		ftnlen)1, (ftnlen)6, (ftnlen)38);
	goto L300;
    }
/*                 if gdel positive, V not posdef */
    if (gdel >= 0.) {
	mnwarn_("D", "MIGRAD", " NEWTON STEP NOT DESCENT.", (ftnlen)1, (
		ftnlen)6, (ftnlen)25);
	if (npsdf == 1) {
	    goto L1;
	}
	mnpsdf_();
	npsdf = 1;
	goto L24;
    }
/*                                        . . . . do line search */
    mnline_((S_fp)fcn, xxs, &fs, step, &gdel, &c_b1209, (U_fp)futil);
    if (mn7min_1.amin == fs) {
	goto L200;
    }
    s_copy(mn7tit_1.cfrom, "MIGRAD  ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = nfcnmg;
    s_copy(mn7tit_1.cstatu, "PROGRESS  ", (ftnlen)10, (ftnlen)10);
/*                                        . get gradient at new point */
    mninex_(mn7int_1.x);
    if (mn7flg_1.isw[2] == 1) {
	(*fcn)(&nparx, mn7der_1.gin, &fzero, mn7ext_1.u, &c__2, (U_fp)futil);
	++mn7cnv_1.nfcn;
    }
    mnderi_((S_fp)fcn, (U_fp)futil);
/*                                         . calculate new EDM */
    npsdf = 0;
L81:
    mn7min_1.edm = (float)0.;
    gvg = (float)0.;
    delgam = (float)0.;
    gdgssq = (float)0.;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ri = (float)0.;
	vgi = (float)0.;
	i__2 = mn7npr_1.npar;
	for (j = 1; j <= i__2; ++j) {
	    m = max(i__,j);
	    n = min(i__,j);
	    ndex = m * (m - 1) / 2 + n;
	    vgi += mn7var_1.vhmat[ndex - 1] * (mn7der_1.grd[j - 1] - gs[j - 1]
		    );
/* L90: */
	    ri += mn7var_1.vhmat[ndex - 1] * mn7der_1.grd[j - 1];
	}
	vg[i__ - 1] = vgi * (float).5;
	gami = mn7der_1.grd[i__ - 1] - gs[i__ - 1];
/* Computing 2nd power */
	d__1 = gami;
	gdgssq += d__1 * d__1;
	gvg += gami * vg[i__ - 1];
	delgam += mn7int_1.dirin[i__ - 1] * gami;
/* L100: */
	mn7min_1.edm += mn7der_1.grd[i__ - 1] * ri * (float).5;
    }
    mn7min_1.edm = mn7min_1.edm * (float).5 * (mn7min_1.dcovar * (float)3. + (
	    float)1.);
/*                          . if EDM negative,  not positive-definite */
    if (mn7min_1.edm < 0. || gvg <= 0.) {
	mnwarn_("D", "MIGRAD", "NOT POS-DEF. EDM OR GVG NEGATIVE.", (ftnlen)1,
		 (ftnlen)6, (ftnlen)33);
	s_copy(mn7tit_1.cstatu, "NOT POSDEF", (ftnlen)10, (ftnlen)10);
	if (npsdf == 1) {
	    goto L230;
	}
	mnpsdf_();
	npsdf = 1;
	goto L81;
    }
/*                            print information about this iteration */
    ++iter;
    if (iswtr >= 3 || iswtr == 2 && iter % 10 == 1) {
	mnwerr_();
	mnprin_(&c__3, &mn7min_1.amin);
    }
    if (gdgssq == 0.) {
	mnwarn_("D", "MIGRAD", "NO CHANGE IN FIRST DERIVATIVES OVER LAST STEP"
		, (ftnlen)1, (ftnlen)6, (ftnlen)45);
    }
    if (delgam < 0.) {
	mnwarn_("D", "MIGRAD", "FIRST DERIVATIVES INCREASING ALONG SEARCH LI\
NE", (ftnlen)1, (ftnlen)6, (ftnlen)46);
    }
/*                                        .  update covariance matrix */
    s_copy(mn7tit_1.cstatu, "IMPROVEMNT", (ftnlen)10, (ftnlen)10);
    if (ldebug) {
	io___710.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___710);
	do_fio(&c__1, " VHMAT 1 =", (ftnlen)10);
	for (kk = 1; kk <= 10; ++kk) {
	    do_fio(&c__1, (char *)&mn7var_1.vhmat[kk - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
    }
    dsum = (float)0.;
    vsum = (float)0.;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    d__ = mn7int_1.dirin[i__ - 1] * mn7int_1.dirin[j - 1] / delgam - 
		    vg[i__ - 1] * vg[j - 1] / gvg;
	    dsum += abs(d__);
	    ndex = i__ * (i__ - 1) / 2 + j;
	    mn7var_1.vhmat[ndex - 1] += d__ * (float)2.;
	    vsum += (d__1 = mn7var_1.vhmat[ndex - 1], abs(d__1));
/* L120: */
	}
    }
/*                smooth local fluctuations by averaging DCOVAR */
    mn7min_1.dcovar = (mn7min_1.dcovar + dsum / vsum) * (float).5;
    if (iswtr >= 3 || ldebug) {
	io___714.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___714);
	do_fio(&c__1, " RELATIVE CHANGE IN COV. MATRIX=", (ftnlen)32);
	d__1 = mn7min_1.dcovar * (float)100.;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "%", (ftnlen)1);
	e_wsfe();
    }
    if (ldebug) {
	io___715.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___715);
	do_fio(&c__1, " VHMAT 2 =", (ftnlen)10);
	for (kk = 1; kk <= 10; ++kk) {
	    do_fio(&c__1, (char *)&mn7var_1.vhmat[kk - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
    }
    if (delgam <= gvg) {
	goto L135;
    }
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L125: */
	flnu[i__ - 1] = mn7int_1.dirin[i__ - 1] / delgam - vg[i__ - 1] / gvg;
    }
    i__2 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = i__;
	for (j = 1; j <= i__1; ++j) {
	    ndex = i__ * (i__ - 1) / 2 + j;
/* L130: */
	    mn7var_1.vhmat[ndex - 1] += gvg * (float)2. * flnu[i__ - 1] * 
		    flnu[j - 1];
	}
    }
L135:
/*                                              and see if converged */
    if (mn7min_1.edm < rhotol * (float).1) {
	goto L300;
    }
/*                                    if not, prepare next iteration */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xxs[i__ - 1] = mn7int_1.x[i__ - 1];
	gs[i__ - 1] = mn7der_1.grd[i__ - 1];
/* L140: */
    }
    fs = mn7min_1.amin;
    if (mn7flg_1.isw[1] == 0 && mn7min_1.dcovar < (float).5) {
	mn7flg_1.isw[1] = 1;
    }
    if (mn7flg_1.isw[1] == 3 && mn7min_1.dcovar > (float).1) {
	mn7flg_1.isw[1] = 1;
    }
    if (mn7flg_1.isw[1] == 1 && mn7min_1.dcovar < (float).05) {
	mn7flg_1.isw[1] = 3;
    }
    goto L24;
/*                                        . . . . .  end main loop */
/*                                         . . call limit in MNMIGR */
L190:
    mn7flg_1.isw[0] = 1;
    if (mn7flg_1.isw[4] >= 0) {
	io___717.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___717);
	do_fio(&c__1, " CALL LIMIT EXCEEDED IN MIGRAD.", (ftnlen)31);
	e_wsfe();
    }
    s_copy(mn7tit_1.cstatu, "CALL LIMIT", (ftnlen)10, (ftnlen)10);
    goto L230;
/*                                         . . fails to improve . . */
L200:
    if (iswtr >= 1) {
	io___718.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___718);
	do_fio(&c__1, " MIGRAD FAILS TO FIND IMPROVEMENT", (ftnlen)33);
	e_wsfe();
    }
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L210: */
	mn7int_1.x[i__ - 1] = xxs[i__ - 1];
    }
    if (mn7min_1.edm < rhotol) {
	goto L300;
    }
    if (mn7min_1.edm < (d__1 = mn7cns_1.epsma2 * mn7min_1.amin, abs(d__1))) {
	if (iswtr >= 0) {
	    io___719.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___719);
	    do_fio(&c__1, " MACHINE ACCURACY LIMITS FURTHER IMPROVEMENT.", (
		    ftnlen)45);
	    e_wsfe();
	}
	goto L300;
    }
    if (mn7cnv_1.istrat < 1) {
	if (mn7flg_1.isw[4] >= 0) {
	    io___720.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___720);
	    do_fio(&c__1, " MIGRAD FAILS WITH STRATEGY=0.   WILL TRY WITH ST\
RATEGY=1.", (ftnlen)58);
	    e_wsfe();
	}
	mn7cnv_1.istrat = 1;
    }
    goto L1;
/*                                         . . fails to converge */
L230:
    if (iswtr >= 0) {
	io___721.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___721);
	do_fio(&c__1, " MIGRAD TERMINATED WITHOUT CONVERGENCE.", (ftnlen)39);
	e_wsfe();
    }
    if (mn7flg_1.isw[1] == 3) {
	mn7flg_1.isw[1] = 1;
    }
    mn7flg_1.isw[3] = -1;
    goto L400;
/*                                         . . apparent convergence */
L300:
    if (iswtr >= 0) {
	io___722.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___722);
	do_fio(&c__1, " MIGRAD MINIMIZATION HAS CONVERGED.", (ftnlen)35);
	e_wsfe();
    }
    if (mn7cnv_1.itaur == 0) {
	if (mn7cnv_1.istrat >= 2 || mn7cnv_1.istrat == 1 && mn7flg_1.isw[1] < 
		3) {
	    if (mn7flg_1.isw[4] >= 0) {
		io___723.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___723);
		do_fio(&c__1, " MIGRAD WILL VERIFY CONVERGENCE AND ERROR MAT\
RIX.", (ftnlen)49);
		e_wsfe();
	    }
	    mnhess_((S_fp)fcn, (U_fp)futil);
	    mnwerr_();
	    npsdf = 0;
	    if (mn7min_1.edm > rhotol) {
		goto L10;
	    }
	}
    }
    s_copy(mn7tit_1.cstatu, "CONVERGED ", (ftnlen)10, (ftnlen)10);
    mn7flg_1.isw[3] = 1;
/*                                           come here in any case */
L400:
    s_copy(mn7tit_1.cfrom, "MIGRAD  ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = nfcnmg;
    mninex_(mn7int_1.x);
    mnwerr_();
    if (iswtr >= 0) {
	mnprin_(&c__3, &mn7min_1.amin);
    }
    if (iswtr >= 1) {
	mnmatu_(&c__1);
    }
    return 0;
} /* mnmigr_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnmnos_(S_fp fcn, U_fp futil)
{
    /* Format strings */
    static char fmt_564[] = "(\002 PARAMETER NUMBER \002,i5,\002 NOT VARIABL\
E. IGNORED.\002)";
    static char fmt_675[] = "(/\002 NEW MINIMUM FOUND.  GO BACK TO MINIMIZAT\
ION STEP.\002/\002 \002,60(\002=\002)/60x,\002V\002/60x,\002V\002/60x,\002\
V\002/57x,\002VVVVVVV\002/58x,\002VVVVV\002/59x,\002VVV\002/60x,\002V\002//)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer iin, knt, nbad, ilax, ilax2, ngood;
    static doublereal val2mi, val2pl;
    static integer nfcnmi;
    extern /* Subroutine */ int mnmnot_(S_fp, integer *, integer *, 
	    doublereal *, doublereal *, U_fp), mnprin_(integer *, doublereal *
	    ), mnmatu_(integer *);

    /* Fortran I/O blocks */
    static cilist io___729 = { 0, 0, 0, fmt_564, 0 };
    static cilist io___734 = { 0, 0, 0, fmt_675, 0 };
    static cilist io___735 = { 0, 0, 0, "(A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Performs a MINOS error analysis on those parameters for */
/* C        which it is requested on the MINOS command by calling */
/* C        MNMNOT for each parameter requested. */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    if (mn7npr_1.npar <= 0) {
	goto L700;
    }
    ngood = 0;
    nbad = 0;
    nfcnmi = mn7cnv_1.nfcn;
/*                                      . loop over parameters requested */
    i__1 = mn7npr_1.npar;
    for (knt = 1; knt <= i__1; ++knt) {
	if ((integer) mn7arg_1.word7[1] == 0) {
	    ilax = mn7inx_1.nexofi[knt - 1];
	} else {
	    if (knt >= 7) {
		goto L580;
	    }
	    ilax = (integer) mn7arg_1.word7[knt];
	    if (ilax == 0) {
		goto L580;
	    }
	    if (ilax > 0 && ilax <= mn7npr_1.nu) {
		if (mn7inx_1.niofex[ilax - 1] > 0) {
		    goto L565;
		}
	    }
	    io___729.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___729);
	    do_fio(&c__1, (char *)&ilax, (ftnlen)sizeof(integer));
	    e_wsfe();
	    goto L570;
	}
L565:
/*                                         calculate one pair of M E's */
	ilax2 = 0;
	mnmnot_((S_fp)fcn, &ilax, &ilax2, &val2pl, &val2mi, (U_fp)futil);
	if (mn7log_1.lnewmn) {
	    goto L650;
	}
/*                                          update NGOOD and NBAD */
	iin = mn7inx_1.niofex[ilax - 1];
	if (mn7err_1.erp[iin - 1] > 0.) {
	    ++ngood;
	} else {
	    ++nbad;
	}
	if (mn7err_1.ern[iin - 1] < 0.) {
	    ++ngood;
	} else {
	    ++nbad;
	}
L570:
	;
    }
/*                                           end of loop . . . . . . . */
L580:
/*                                        . . . . printout final values . */
    s_copy(mn7tit_1.cfrom, "MINOS   ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = nfcnmi;
    s_copy(mn7tit_1.cstatu, "UNCHANGED ", (ftnlen)10, (ftnlen)10);
    if (ngood == 0 && nbad == 0) {
	goto L700;
    }
    if (ngood > 0 && nbad == 0) {
	s_copy(mn7tit_1.cstatu, "SUCCESSFUL", (ftnlen)10, (ftnlen)10);
    }
    if (ngood == 0 && nbad > 0) {
	s_copy(mn7tit_1.cstatu, "FAILURE   ", (ftnlen)10, (ftnlen)10);
    }
    if (ngood > 0 && nbad > 0) {
	s_copy(mn7tit_1.cstatu, "PROBLEMS  ", (ftnlen)10, (ftnlen)10);
    }
    if (mn7flg_1.isw[4] >= 0) {
	mnprin_(&c__4, &mn7min_1.amin);
    }
    if (mn7flg_1.isw[4] >= 2) {
	mnmatu_(&c__0);
    }
    goto L900;
/*                                        . . . new minimum found . . . . */
L650:
    s_copy(mn7tit_1.cfrom, "MINOS   ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = nfcnmi;
    s_copy(mn7tit_1.cstatu, "NEW MINIMU", (ftnlen)10, (ftnlen)10);
    if (mn7flg_1.isw[4] >= 0) {
	mnprin_(&c__4, &mn7min_1.amin);
    }
    io___734.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___734);
    e_wsfe();
    goto L900;
L700:
    io___735.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___735);
    do_fio(&c__1, " THERE ARE NO MINOS ERRORS TO CALCULATE.", (ftnlen)40);
    e_wsfe();
L900:
    return 0;
} /* mnmnos_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int mnmnot_(S_fp fcn, integer *ilax, integer *ilax2, 
	doublereal *val2pl, doublereal *val2mi, U_fp futil)
{
    /* Format strings */
    static char fmt_806[] = "(/\002 DETERMINATION OF \002,a4,\002TIVE MINOS \
ERROR FOR PARAMETER\002,i3,2x,a)";
    static char fmt_801[] = "(/\002 PARAMETER\002,i4,\002 SET TO\002,e11.3\
,\002 + \002,e10.3,\002 = \002,e12.3)";
    static char fmt_808[] = "(/9x,\002THE \002,a4,\002TIVE MINOS ERROR OF PA\
RAMETER\002,i3,\002 ,\002,a10,\002, IS\002,e12.4)";
    static char fmt_807[] = "(5x,\002THE \002,a4,\002TIVE MINOS ERROR OF PAR\
AMETER\002,i3,\002, \002,a,\002, EXCEEDS ITS LIMIT.\002/)";
    static char fmt_802[] = "(9x,\002THE \002,a,\002TIVE MINOS ERROR\002,i4\
,\002 REQUIRES MORE THAN\002,i5,\002 FUNCTION CALLS.\002/)";
    static char fmt_805[] = "(25x,a,\002TIVE MINOS ERROR NOT CALCULATED FOR \
PARAMETER\002,i4/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;
    static doublereal w[100], dc;
    static integer it;
    static doublereal ut, du1, fac, gcc[100], sig, sav;
    static integer isw2, isw4;
    static char csig[4];
    static integer marc;
    static doublereal delu;
    static integer isig, mpar, ndex;
    static doublereal xdev[100];
    static integer imax, indx, ierr;
    static doublereal aopt, eros, abest;
    static integer iercr;
    static doublereal xunit;
    extern /* Subroutine */ int mnfree_(integer *);
    static doublereal sigsav;
    static integer istrav, nfmxin;
    extern /* Subroutine */ int mninex_(doublereal *), mnfixp_(integer *, 
	    integer *), mnwarn_(char *, char *, char *, ftnlen, ftnlen, 
	    ftnlen);
    static integer nlimit;
    extern /* Subroutine */ int mncros_(S_fp, doublereal *, integer *, U_fp), 
	    mnexin_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___757 = { 0, 0, 0, "(A,I5,A,I5)", 0 };
    static cilist io___761 = { 0, 0, 0, fmt_806, 0 };
    static cilist io___766 = { 0, 0, 0, fmt_801, 0 };
    static cilist io___770 = { 0, 0, 0, fmt_808, 0 };
    static cilist io___771 = { 0, 0, 0, fmt_807, 0 };
    static cilist io___772 = { 0, 0, 0, fmt_802, 0 };
    static cilist io___773 = { 0, 0, 0, fmt_805, 0 };
    static cilist io___774 = { 0, 0, 0, "(5X, 74(1H*))", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Performs a MINOS error analysis on one parameter. */
/* C        The parameter ILAX is varied, and the minimum of the */
/* C        function with respect to the other parameters is followed */
/* C        until it crosses the value FMIN+UP. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/*                                        . . save and prepare start vals */
    isw2 = mn7flg_1.isw[1];
    isw4 = mn7flg_1.isw[3];
    sigsav = mn7min_1.edm;
    istrav = mn7cnv_1.istrat;
    dc = mn7min_1.dcovar;
    mn7log_1.lnewmn = FALSE_;
    mn7min_1.apsi = mn7min_1.epsi * (float).5;
    abest = mn7min_1.amin;
    mpar = mn7npr_1.npar;
    nfmxin = mn7cnv_1.nfcnmx;
    i__1 = mpar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L125: */
	mn7int_1.xt[i__ - 1] = mn7int_1.x[i__ - 1];
    }
    i__1 = mpar * (mpar + 1) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L130: */
	mn7vat_1.vthmat[j - 1] = mn7var_1.vhmat[j - 1];
    }
    i__1 = mpar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gcc[i__ - 1] = mn7err_1.globcc[i__ - 1];
/* L135: */
	w[i__ - 1] = mn7err_1.werr[i__ - 1];
    }
    it = mn7inx_1.niofex[*ilax - 1];
    mn7err_1.erp[it - 1] = (float)0.;
    mn7err_1.ern[it - 1] = (float)0.;
    mninex_(mn7int_1.xt);
    ut = mn7ext_1.u[*ilax - 1];
    if (mn7inx_1.nvarl[*ilax - 1] == 1) {
	mn7ext_1.alim[*ilax - 1] = ut - w[it - 1] * (float)100.;
	mn7ext_1.blim[*ilax - 1] = ut + w[it - 1] * (float)100.;
    }
    ndex = it * (it + 1) / 2;
    xunit = sqrt(mn7min_1.up / mn7vat_1.vthmat[ndex - 1]);
    marc = 0;
    i__1 = mpar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == it) {
	    goto L162;
	}
	++marc;
	imax = max(it,i__);
	indx = imax * (imax - 1) / 2 + min(it,i__);
	xdev[marc - 1] = xunit * mn7vat_1.vthmat[indx - 1];
L162:
	;
    }
/*                           fix the parameter in question */
    mnfixp_(&it, &ierr);
    if (ierr > 0) {
	io___757.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___757);
	do_fio(&c__1, " MINUIT ERROR. CANNOT FIX PARAMETER", (ftnlen)35);
	do_fio(&c__1, (char *)&(*ilax), (ftnlen)sizeof(integer));
	do_fio(&c__1, "    INTERNAL", (ftnlen)12);
	do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L700;
    }
/*                       . . . . . Nota Bene: from here on, NPAR=MPAR-1 */
/*      Remember: MNFIXP squeezes IT out of X, XT, WERR, and VHMAT, */
/*                                                    not W, VTHMAT */
    for (isig = 1; isig <= 2; ++isig) {
	if (isig == 1) {
	    sig = (float)1.;
	    s_copy(csig, "POSI", (ftnlen)4, (ftnlen)4);
	} else {
	    sig = (float)-1.;
	    s_copy(csig, "NEGA", (ftnlen)4, (ftnlen)4);
	}
/*                                        . sig=sign of error being calcd */
	if (mn7flg_1.isw[4] > 1) {
	    io___761.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___761);
	    do_fio(&c__1, csig, (ftnlen)4);
	    do_fio(&c__1, (char *)&(*ilax), (ftnlen)sizeof(integer));
	    do_fio(&c__1, mn7nam_1.cpnam + (*ilax - 1) * 10, (ftnlen)10);
	    e_wsfe();
	}
	if (mn7flg_1.isw[1] <= 0) {
	    mnwarn_("D", "MINOS", "NO COVARIANCE MATRIX.", (ftnlen)1, (ftnlen)
		    5, (ftnlen)21);
	}
	nlimit = mn7cnv_1.nfcn + nfmxin;
/* Computing MAX */
	i__1 = istrav - 1;
	mn7cnv_1.istrat = max(i__1,0);
	du1 = w[it - 1];
	mn7ext_1.u[*ilax - 1] = ut + sig * du1;
/* Computing MIN */
	d__1 = mn7ext_1.u[*ilax - 1], d__2 = mn7ext_1.blim[*ilax - 1];
	mn7ext_1.u[*ilax - 1] = min(d__1,d__2);
/* Computing MAX */
	d__1 = mn7ext_1.u[*ilax - 1], d__2 = mn7ext_1.alim[*ilax - 1];
	mn7ext_1.u[*ilax - 1] = max(d__1,d__2);
	delu = mn7ext_1.u[*ilax - 1] - ut;
/*         stop if already at limit with negligible step size */
	if (abs(delu) / (abs(ut) + (d__1 = mn7ext_1.u[*ilax - 1], abs(d__1))) 
		< mn7cns_1.epsmac) {
	    goto L440;
	}
	fac = delu / w[it - 1];
	i__1 = mn7npr_1.npar;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L185: */
	    mn7int_1.x[i__ - 1] = mn7int_1.xt[i__ - 1] + fac * xdev[i__ - 1];
	}
	if (mn7flg_1.isw[4] > 1) {
	    io___766.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___766);
	    do_fio(&c__1, (char *)&(*ilax), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ut, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&delu, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&mn7ext_1.u[*ilax - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
/*                                        loop to hit AMIN+UP */
	mn7xcr_1.ke1cr = *ilax;
	mn7xcr_1.ke2cr = 0;
	mn7xcr_1.xmidcr = mn7ext_1.u[*ilax - 1];
	mn7xcr_1.xdircr = delu;

	mn7min_1.amin = abest;
	mn7cnv_1.nfcnmx = nlimit - mn7cnv_1.nfcn;
	mncros_((S_fp)fcn, &aopt, &iercr, (U_fp)futil);
	if (abest - mn7min_1.amin > mn7min_1.up * (float).01) {
	    goto L650;
	}
	if (iercr == 1) {
	    goto L440;
	}
	if (iercr == 2) {
	    goto L450;
	}
	if (iercr == 3) {
	    goto L460;
	}
/*                                        . error successfully calculated */
	eros = mn7xcr_1.xmidcr - ut + aopt * mn7xcr_1.xdircr;
	if (mn7flg_1.isw[4] > 1) {
	    io___770.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___770);
	    do_fio(&c__1, csig, (ftnlen)4);
	    do_fio(&c__1, (char *)&(*ilax), (ftnlen)sizeof(integer));
	    do_fio(&c__1, mn7nam_1.cpnam + (*ilax - 1) * 10, (ftnlen)10);
	    do_fio(&c__1, (char *)&eros, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	goto L480;
/*                                        . . . . . . . . failure returns */
L440:
	if (mn7flg_1.isw[4] >= 1) {
	    io___771.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___771);
	    do_fio(&c__1, csig, (ftnlen)4);
	    do_fio(&c__1, (char *)&(*ilax), (ftnlen)sizeof(integer));
	    do_fio(&c__1, mn7nam_1.cpnam + (*ilax - 1) * 10, (ftnlen)10);
	    e_wsfe();
	}
	eros = mn7cns_1.undefi;
	goto L480;
L450:
	if (mn7flg_1.isw[4] >= 1) {
	    io___772.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___772);
	    do_fio(&c__1, csig, (ftnlen)4);
	    do_fio(&c__1, (char *)&(*ilax), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nfmxin, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	eros = (float)0.;
	goto L480;
L460:
	if (mn7flg_1.isw[4] >= 1) {
	    io___773.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___773);
	    do_fio(&c__1, csig, (ftnlen)4);
	    do_fio(&c__1, (char *)&(*ilax), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	eros = (float)0.;

L480:
	if (mn7flg_1.isw[4] > 1) {
	    io___774.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___774);
	    e_wsfe();
	}
	if (sig < 0.) {
	    mn7err_1.ern[it - 1] = eros;
	    if (*ilax2 > 0 && *ilax2 <= mn7npr_1.nu) {
		*val2mi = mn7ext_1.u[*ilax2 - 1];
	    }
	} else {
	    mn7err_1.erp[it - 1] = eros;
	    if (*ilax2 > 0 && *ilax2 <= mn7npr_1.nu) {
		*val2pl = mn7ext_1.u[*ilax2 - 1];
	    }
	}
/* L500: */
    }
/*                                        . . parameter finished. reset v */
/*                       normal termination */
    mn7cnv_1.itaur = 1;
    mnfree_(&c__1);
    i__1 = mpar * (mpar + 1) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L550: */
	mn7var_1.vhmat[j - 1] = mn7vat_1.vthmat[j - 1];
    }
    i__1 = mpar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7err_1.werr[i__ - 1] = w[i__ - 1];
	mn7err_1.globcc[i__ - 1] = gcc[i__ - 1];
/* L595: */
	mn7int_1.x[i__ - 1] = mn7int_1.xt[i__ - 1];
    }
    mninex_(mn7int_1.x);
    mn7min_1.edm = sigsav;
    mn7min_1.amin = abest;
    mn7flg_1.isw[1] = isw2;
    mn7flg_1.isw[3] = isw4;
    mn7min_1.dcovar = dc;
    goto L700;
/*                       new minimum */
L650:
    mn7log_1.lnewmn = TRUE_;
    mn7flg_1.isw[1] = 0;
    mn7min_1.dcovar = (float)1.;
    mn7flg_1.isw[3] = 0;
    sav = mn7ext_1.u[*ilax - 1];
    mn7cnv_1.itaur = 1;
    mnfree_(&c__1);
    mn7ext_1.u[*ilax - 1] = sav;
    mnexin_(mn7int_1.x);
    mn7min_1.edm = mn7cns_1.bigedm;
/*                       in any case */
L700:
    mn7cnv_1.itaur = 0;
    mn7cnv_1.nfcnmx = nfmxin;
    mn7cnv_1.istrat = istrav;
    return 0;
} /* mnmnot_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:50  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnparm_(integer *k, char *cnamj, doublereal *uk, 
	doublereal *wk, doublereal *a, doublereal *b, integer *ierflg, ftnlen 
	cnamj_len)
{
    /* Format strings */
    static char fmt_9[] = "(/\002 MINUIT USER ERROR.  PARAMETER NUMBER IS\
\002,i11/\002,  ALLOWED RANGE IS ONE TO\002,i4/)";
    static char fmt_61[] = "(/\002 PARAMETER DEFINITIONS:\002/\002    NO.   \
NAME         VALUE      STEP SIZE      LIMITS\002)";
    static char fmt_82[] = "(1x,i5,1x,\002'\002,a10,\002'\002,1x,g13.5,\002 \
 constant\002)";
    static char fmt_127[] = "(1x,i5,1x,\002'\002,a10,\002'\002,1x,2g13.5,\
\002     no limits\002)";
    static char fmt_132[] = "(1x,i5,1x,\002'\002,a10,\002'\002,1x,2g13.5,2x,\
2g13.5)";
    static char fmt_135[] = "(/\002 MINUIT USER ERROR.   TOO MANY VARIABLE P\
ARAMETERS.\002/\002 THIS VERSION OF MINUIT DIMENSIONED FOR\002,i4//)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2[3];
    doublereal d__1, d__2;
    char ch__1[34];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_wsfi(icilist *), e_wsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer in, ix;
    static doublereal sav;
    static integer nvl;
    static doublereal sav2;
    static integer ierr, kint;
    static doublereal vplu;
    static char cnamk[10];
    static doublereal small, gsmin, pinti, vminu;
    static char chbufi[4];
    static doublereal danger;
    extern /* Subroutine */ int mnfree_(integer *);
    static integer ktofix;
    extern /* Subroutine */ int mnwarn_(char *, char *, char *, ftnlen, 
	    ftnlen, ftnlen);
    static integer lastin;
    extern /* Subroutine */ int mnrset_(integer *), mnpint_(doublereal *, 
	    integer *, doublereal *);
    static integer kinfix;
    extern /* Subroutine */ int mnfixp_(integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___778 = { 0, 0, 0, fmt_9, 0 };
    static cilist io___781 = { 0, 0, 0, "(A)", 0 };
    static cilist io___782 = { 0, 0, 0, fmt_61, 0 };
    static cilist io___783 = { 0, 0, 0, fmt_82, 0 };
    static cilist io___785 = { 0, 0, 0, fmt_127, 0 };
    static cilist io___786 = { 0, 0, 0, fmt_132, 0 };
    static cilist io___787 = { 0, 0, 0, fmt_135, 0 };
    static cilist io___788 = { 0, 0, 0, "(/A,A/A/)", 0 };
    static icilist io___791 = { 0, chbufi, 0, "(I4)", 4, 1 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called from MNPARS and user-callable */
/* C    Implements one parameter definition, that is: */
/* C          K     (external) parameter number */
/* C          CNAMK parameter name */
/* C          UK    starting value */
/* C          WK    starting step size or uncertainty */
/* C          A, B  lower and upper physical parameter limits */
/* C    and sets up (updates) the parameter lists. */
/* C    Output: IERFLG=0 if no problems */
/* C                  >0 if MNPARM unable to implement definition */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



    s_copy(cnamk, cnamj, (ftnlen)10, cnamj_len);
    kint = mn7npr_1.npar;
    if (*k < 1 || *k > mn7npr_1.maxext) {
/*                     parameter number exceeds allowed maximum value */
	io___778.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___778);
	do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&mn7npr_1.maxext, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L800;
    }
/*                     normal parameter request */
    ktofix = 0;
    if (mn7inx_1.nvarl[*k - 1] < 0) {
	goto L50;
    }
/*         previously defined parameter is being redefined */
/*                                     find if parameter was fixed */
    i__1 = mn7fx1_1.npfix;
    for (ix = 1; ix <= i__1; ++ix) {
	if (mn7fx1_1.ipfix[ix - 1] == *k) {
	    ktofix = *k;
	}
/* L40: */
    }
    if (ktofix > 0) {
	mnwarn_("W", "PARAM DEF", "REDEFINING A FIXED PARAMETER.", (ftnlen)1, 
		(ftnlen)9, (ftnlen)29);
	if (kint >= mn7npr_1.maxint) {
	    io___781.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___781);
	    do_fio(&c__1, " CANNOT RELEASE. MAX NPAR EXCEEDED.", (ftnlen)35);
	    e_wsfe();
	    goto L800;
	}
	i__1 = -(*k);
	mnfree_(&i__1);
    }
/*                       if redefining previously variable parameter */
    if (mn7inx_1.niofex[*k - 1] > 0) {
	kint = mn7npr_1.npar - 1;
    }
L50:

/*                                      . . .print heading */
    if (mn7log_1.lphead && mn7flg_1.isw[4] >= 0) {
	io___782.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___782);
	e_wsfe();
	mn7log_1.lphead = FALSE_;
    }
    if (*wk > 0.) {
	goto L122;
    }
/*                                        . . .constant parameter . . . . */
    if (mn7flg_1.isw[4] >= 0) {
	io___783.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___783);
	do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
	do_fio(&c__1, cnamk, (ftnlen)10);
	do_fio(&c__1, (char *)&(*uk), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    nvl = 0;
    goto L200;
L122:
    if (*a == 0. && *b == 0.) {
/*                                      variable parameter without limits */
	nvl = 1;
	if (mn7flg_1.isw[4] >= 0) {
	    io___785.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___785);
	    do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
	    do_fio(&c__1, cnamk, (ftnlen)10);
	    do_fio(&c__1, (char *)&(*uk), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*wk), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else {
/*                                         variable parameter with limits */
	nvl = 4;
	mn7log_1.lnolim = FALSE_;
	if (mn7flg_1.isw[4] >= 0) {
	    io___786.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___786);
	    do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
	    do_fio(&c__1, cnamk, (ftnlen)10);
	    do_fio(&c__1, (char *)&(*uk), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*wk), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*a), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*b), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/*                             . . request for another variable parameter */
    ++kint;
    if (kint > mn7npr_1.maxint) {
	io___787.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___787);
	do_fio(&c__1, (char *)&mn7npr_1.maxint, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L800;
    }
    if (nvl == 1) {
	goto L200;
    }
    if (*a == *b) {
	io___788.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___788);
	do_fio(&c__1, " USER ERROR IN MINUIT PARAMETER", (ftnlen)31);
	do_fio(&c__1, " DEFINITION", (ftnlen)11);
	do_fio(&c__1, " UPPER AND LOWER LIMITS EQUAL.", (ftnlen)30);
	e_wsfe();
	goto L800;
    }
    if (*b < *a) {
	sav = *b;
	*b = *a;
	*a = sav;
	mnwarn_("W", "PARAM DEF", "PARAMETER LIMITS WERE REVERSED.", (ftnlen)
		1, (ftnlen)9, (ftnlen)31);
	if (mn7log_1.lwarn) {
	    mn7log_1.lphead = TRUE_;
	}
    }
    if (*b - *a > (float)1e7) {
	s_wsfi(&io___791);
	do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 15, a__1[0] = "LIMITS ON PARAM";
	i__2[1] = 4, a__1[1] = chbufi;
	i__2[2] = 15, a__1[2] = " TOO FAR APART.";
	s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)34);
	mnwarn_("W", "PARAM DEF", ch__1, (ftnlen)1, (ftnlen)9, (ftnlen)34);
	if (mn7log_1.lwarn) {
	    mn7log_1.lphead = TRUE_;
	}
    }
    danger = (*b - *uk) * (*uk - *a);
    if (danger < (float)0.) {
	mnwarn_("W", "PARAM DEF", "STARTING VALUE OUTSIDE LIMITS.", (ftnlen)1,
		 (ftnlen)9, (ftnlen)30);
    }
    if (danger == (float)0.) {
	mnwarn_("W", "PARAM DEF", "STARTING VALUE IS AT LIMIT.", (ftnlen)1, (
		ftnlen)9, (ftnlen)27);
    }
L200:
/*                           . . . input OK, set values, arrange lists, */
/*                                    calculate step sizes GSTEP, DIRIN */
    s_copy(mn7tit_1.cfrom, "PARAMETR", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "NEW VALUES", (ftnlen)10, (ftnlen)10);
    mn7npr_1.nu = max(mn7npr_1.nu,*k);
    s_copy(mn7nam_1.cpnam + (*k - 1) * 10, cnamk, (ftnlen)10, (ftnlen)10);
    mn7ext_1.u[*k - 1] = *uk;
    mn7ext_1.alim[*k - 1] = *a;
    mn7ext_1.blim[*k - 1] = *b;
    mn7inx_1.nvarl[*k - 1] = nvl;
/*                             K is external number of new parameter */
/*           LASTIN is the number of var. params with ext. param. no.< K */
    lastin = 0;
    i__1 = *k - 1;
    for (ix = 1; ix <= i__1; ++ix) {
	if (mn7inx_1.niofex[ix - 1] > 0) {
	    ++lastin;
	}
/* L240: */
    }
/*                 KINT is new number of variable params, NPAR is old */
    if (kint == mn7npr_1.npar) {
	goto L280;
    }
    if (kint > mn7npr_1.npar) {
/*                          insert new variable parameter in list */
	i__1 = lastin + 1;
	for (in = mn7npr_1.npar; in >= i__1; --in) {
	    ix = mn7inx_1.nexofi[in - 1];
	    mn7inx_1.niofex[ix - 1] = in + 1;
	    mn7inx_1.nexofi[in] = ix;
	    mn7int_1.x[in] = mn7int_1.x[in - 1];
	    mn7int_1.xt[in] = mn7int_1.xt[in - 1];
	    mn7int_1.dirin[in] = mn7int_1.dirin[in - 1];
	    mn7der_1.g2[in] = mn7der_1.g2[in - 1];
	    mn7der_1.gstep[in] = mn7der_1.gstep[in - 1];
/* L260: */
	}
    } else {
/*                          remove variable parameter from list */
	i__1 = kint;
	for (in = lastin + 1; in <= i__1; ++in) {
	    ix = mn7inx_1.nexofi[in];
	    mn7inx_1.niofex[ix - 1] = in;
	    mn7inx_1.nexofi[in - 1] = ix;
	    mn7int_1.x[in - 1] = mn7int_1.x[in];
	    mn7int_1.xt[in - 1] = mn7int_1.xt[in];
	    mn7int_1.dirin[in - 1] = mn7int_1.dirin[in];
	    mn7der_1.g2[in - 1] = mn7der_1.g2[in];
	    mn7der_1.gstep[in - 1] = mn7der_1.gstep[in];
/* L270: */
	}
    }
L280:
    ix = *k;
    mn7inx_1.niofex[ix - 1] = 0;
    mn7npr_1.npar = kint;
    mnrset_(&c__1);
/*                                       lists are now arranged . . . . */
    if (nvl > 0) {
	in = lastin + 1;
	mn7inx_1.nexofi[in - 1] = ix;
	mn7inx_1.niofex[ix - 1] = in;
	sav = mn7ext_1.u[ix - 1];
	mnpint_(&sav, &ix, &pinti);
	mn7int_1.x[in - 1] = pinti;
	mn7int_1.xt[in - 1] = mn7int_1.x[in - 1];
	mn7err_1.werr[in - 1] = *wk;
	sav2 = sav + *wk;
	mnpint_(&sav2, &ix, &pinti);
	vplu = pinti - mn7int_1.x[in - 1];
	sav2 = sav - *wk;
	mnpint_(&sav2, &ix, &pinti);
	vminu = pinti - mn7int_1.x[in - 1];
	mn7int_1.dirin[in - 1] = (abs(vplu) + abs(vminu)) * (float).5;
/* Computing 2nd power */
	d__1 = mn7int_1.dirin[in - 1];
	mn7der_1.g2[in - 1] = mn7min_1.up * (float)2. / (d__1 * d__1);
	gsmin = mn7cns_1.epsma2 * (float)8. * (d__1 = mn7int_1.x[in - 1], abs(
		d__1));
/* Computing MAX */
	d__1 = gsmin, d__2 = mn7int_1.dirin[in - 1] * (float).1;
	mn7der_1.gstep[in - 1] = max(d__1,d__2);
	if (mn7min_1.amin != mn7cns_1.undefi) {
	    small = sqrt(mn7cns_1.epsma2 * (mn7min_1.amin + mn7min_1.up) / 
		    mn7min_1.up);
/* Computing MAX */
	    d__1 = gsmin, d__2 = small * mn7int_1.dirin[in - 1];
	    mn7der_1.gstep[in - 1] = max(d__1,d__2);
	}
	mn7der_1.grd[in - 1] = mn7der_1.g2[in - 1] * mn7int_1.dirin[in - 1];
/*                   if parameter has limits */
	if (mn7inx_1.nvarl[*k - 1] > 1) {
	    if (mn7der_1.gstep[in - 1] > (float).5) {
		mn7der_1.gstep[in - 1] = (float).5;
	    }
	    mn7der_1.gstep[in - 1] = -mn7der_1.gstep[in - 1];
	}
    }
    if (ktofix > 0) {
	kinfix = mn7inx_1.niofex[ktofix - 1];
	if (kinfix > 0) {
	    mnfixp_(&kinfix, &ierr);
	}
	if (ierr > 0) {
	    goto L800;
	}
    }
    *ierflg = 0;
    return 0;
/*                   error on input, unable to implement request  . . . . */
L800:
    *ierflg = 1;
    return 0;
} /* mnparm_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnpars_(char *crdbuf, integer *icondn, ftnlen crdbuf_len)
{
    /* Format strings */
    static char fmt_158[] = "(bn,f10.0,a10,4f10.0)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    icilist ici__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen), i_indx(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsfi(icilist *), do_fio(integer *, char *, ftnlen), e_rsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static doublereal a, b;
    static integer k;
    static doublereal fk, uk, wk, xk;
    static integer lnc, icy, ierr, kapo1, kapo2;
    static char cnamk[10];
    static integer llist;
    static doublereal plist[30];
    static integer ibegin;
    static char comand[20], celmnt[20];
    static integer lenbuf;
    extern /* Subroutine */ int mncrck_(char *, integer *, char *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, ftnlen, 
	    ftnlen);
    static integer istart;
    extern /* Subroutine */ int mnparm_(integer *, char *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___808 = { 1, celmnt, 0, "(BN,F20.0)", 20, 1 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called from MNREAD and user-callable */
/* C    Implements one parameter definition, that is: */
/* C       parses the string CRDBUF and calls MNPARM */

/* output conditions: */
/*        ICONDN = 0    all OK */
/*        ICONDN = 1    error, attempt to define parameter is ignored */
/*        ICONDN = 2    end of parameter definitions */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */




    lenbuf = i_len(crdbuf, crdbuf_len);
/*                     find out whether fixed or free-field format */
    kapo1 = i_indx(crdbuf, "'", crdbuf_len, (ftnlen)1);
    if (kapo1 == 0) {
	goto L150;
    }
    i__1 = kapo1;
    kapo2 = i_indx(crdbuf + i__1, "'", crdbuf_len - i__1, (ftnlen)1);
    if (kapo2 == 0) {
	goto L150;
    }
/*          new (free-field) format */
    kapo2 += kapo1;
/*                             skip leading blanks if any */
    i__1 = kapo1 - 1;
    for (istart = 1; istart <= i__1; ++istart) {
	if (*(unsigned char *)&crdbuf[istart - 1] != ' ') {
	    goto L120;
	}
/* L115: */
    }
    goto L210;
L120:
/*                               parameter number integer */
    s_copy(celmnt, crdbuf + (istart - 1), (ftnlen)20, kapo1 - 1 - (istart - 1)
	    );
    i__1 = s_rsfi(&io___808);
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = do_fio(&c__1, (char *)&fk, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = e_rsfi();
    if (i__1 != 0) {
	goto L180;
    }
    k = (integer) fk;
    if (k <= 0) {
	goto L210;
    }
/* Writing concatenation */
    i__2[0] = 6, a__1[0] = "PARAM ";
    i__2[1] = 20, a__1[1] = celmnt;
    s_cat(cnamk, a__1, i__2, &c__2, (ftnlen)10);
    if (kapo2 - kapo1 > 1) {
	i__1 = kapo1;
	s_copy(cnamk, crdbuf + i__1, (ftnlen)10, kapo2 - 1 - i__1);
    }
/*  special handling if comma or blanks and a comma follow 'name' */
    i__1 = lenbuf;
    for (icy = kapo2 + 1; icy <= i__1; ++icy) {
	if (*(unsigned char *)&crdbuf[icy - 1] == ',') {
	    goto L139;
	}
	if (*(unsigned char *)&crdbuf[icy - 1] != ' ') {
	    goto L140;
	}
/* L135: */
    }
    uk = (float)0.;
    wk = (float)0.;
    a = (float)0.;
    b = (float)0.;
    goto L170;
L139:
    ++icy;
L140:
    ibegin = icy;
    mncrck_(crdbuf + (ibegin - 1), &c__20, comand, &lnc, &c__30, plist, &
	    llist, &ierr, &mn7iou_1.isyswr, crdbuf_len - (ibegin - 1), (
	    ftnlen)20);
    if (ierr > 0) {
	goto L180;
    }
    uk = plist[0];
    wk = (float)0.;
    if (llist >= 2) {
	wk = plist[1];
    }
    a = (float)0.;
    if (llist >= 3) {
	a = plist[2];
    }
    b = (float)0.;
    if (llist >= 4) {
	b = plist[3];
    }
    goto L170;
/*          old (fixed-field) format */
L150:
    ici__1.icierr = 1;
    ici__1.iciend = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = crdbuf_len;
    ici__1.iciunit = crdbuf;
    ici__1.icifmt = fmt_158;
    i__1 = s_rsfi(&ici__1);
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = do_fio(&c__1, (char *)&xk, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = do_fio(&c__1, cnamk, (ftnlen)10);
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = do_fio(&c__1, (char *)&uk, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = do_fio(&c__1, (char *)&wk, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = do_fio(&c__1, (char *)&b, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L180;
    }
    i__1 = e_rsfi();
    if (i__1 != 0) {
	goto L180;
    }
    k = (integer) xk;
    if (k == 0) {
	goto L210;
    }
/*          parameter format cracked, implement parameter definition */
L170:
    mnparm_(&k, cnamk, &uk, &wk, &a, &b, &ierr, (ftnlen)10);
    *icondn = ierr;
    return 0;
/*          format or other error */
L180:
    *icondn = 1;
    return 0;
/*        end of data */
L210:
    *icondn = 2;
    return 0;
} /* mnpars_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnpfit_(doublereal *parx2p, doublereal *pary2p, integer *
	npar2p, doublereal *coef2p, doublereal *sdev2p)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a, f;
    static integer i__;
    static doublereal s, t, y, s2, x2, x3, x4, y2, cz[3], xm, xy, x2y;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */

/*     to fit a parabola to npar2p points */

/*   npar2p   no. of points */
/*   parx2p(i)   x value of point i */
/*   pary2p(i)   y value of point i */

/*   coef2p(1...3)  coefficients of the fitted parabola */
/*   y=coef2p(1) + coef2p(2)*x + coef2p(3)*x**2 */
/*   sdev2p= variance */
/*   method : chi**2 = min equation solved explicitly */

    /* Parameter adjustments */
    --coef2p;
    --pary2p;
    --parx2p;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
/* L3: */
	cz[i__ - 1] = (float)0.;
    }
    *sdev2p = (float)0.;
    if (*npar2p < 3) {
	goto L10;
    }
    f = (doublereal) (*npar2p);
/* --- center x values for reasons of machine precision */
    xm = (float)0.;
    i__1 = *npar2p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L2: */
	xm += parx2p[i__];
    }
    xm /= f;
    x2 = (float)0.;
    x3 = (float)0.;
    x4 = (float)0.;
    y = (float)0.;
    y2 = (float)0.;
    xy = (float)0.;
    x2y = (float)0.;
    i__1 = *npar2p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = parx2p[i__] - xm;
	t = pary2p[i__];
	s2 = s * s;
	x2 += s2;
	x3 += s * s2;
	x4 += s2 * s2;
	y += t;
	y2 += t * t;
	xy += s * t;
	x2y += s2 * t;
/* L1: */
    }
/* Computing 2nd power */
    d__1 = x2;
/* Computing 2nd power */
    d__2 = x3;
    a = (f * x4 - d__1 * d__1) * x2 - f * (d__2 * d__2);
    if (a == (float)0.) {
	goto L10;
    }
    cz[2] = (x2 * (f * x2y - x2 * y) - f * x3 * xy) / a;
    cz[1] = (xy - x3 * cz[2]) / x2;
    cz[0] = (y - x2 * cz[2]) / f;
    if (*npar2p == 3) {
	goto L6;
    }
    *sdev2p = y2 - (cz[0] * y + cz[1] * xy + cz[2] * x2y);
    if (*sdev2p < (float)0.) {
	*sdev2p = (float)0.;
    }
    *sdev2p /= f - (float)3.;
L6:
    cz[0] += xm * (xm * cz[2] - cz[1]);
    cz[1] -= xm * (float)2. * cz[2];
L10:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L11: */
	coef2p[i__] = cz[i__ - 1];
    }
    return 0;
} /* mnpfit_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnpint_(doublereal *pexti, integer *i__, doublereal *
	pinti)
{
    /* System generated locals */
    address a__1[3];
    integer i__1[3];
    doublereal d__1;
    char ch__1[42];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sin(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double asin(doublereal);

    /* Local variables */
    static doublereal a, yy, yy2;
    static integer igo;
    static doublereal alimi, blimi;
    static char chbuf2[30], chbufi[4];
    extern /* Subroutine */ int mnwarn_(char *, char *, char *, ftnlen, 
	    ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___847 = { 0, chbufi, 0, "(I4)", 4, 1 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Calculates the internal parameter value PINTI corresponding */
/* C        to the external value PEXTI for parameter I. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    *pinti = *pexti;
    igo = mn7inx_1.nvarl[*i__ - 1];
    if (igo == 4) {
/* --                          there are two limits */
	alimi = mn7ext_1.alim[*i__ - 1];
	blimi = mn7ext_1.blim[*i__ - 1];
	yy = (*pexti - alimi) * (float)2. / (blimi - alimi) - (float)1.;
/* Computing 2nd power */
	d__1 = yy;
	yy2 = d__1 * d__1;
	if (yy2 >= (float)1. - mn7cns_1.epsma2) {
	    if (yy < (float)0.) {
		a = mn7cns_1.vlimlo;
		s_copy(chbuf2, " IS AT ITS LOWER ALLOWED LIMIT.", (ftnlen)30, 
			(ftnlen)31);
	    } else {
		a = mn7cns_1.vlimhi;
		s_copy(chbuf2, " IS AT ITS UPPER ALLOWED LIMIT.", (ftnlen)30, 
			(ftnlen)31);
	    }
	    *pinti = a;
	    *pexti = alimi + (blimi - alimi) * (float).5 * (sin(a) + (float)
		    1.);
	    mn7log_1.limset = TRUE_;
	    s_wsfi(&io___847);
	    do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	    e_wsfi();
	    if (yy2 > (float)1.) {
		s_copy(chbuf2, " BROUGHT BACK INSIDE LIMITS.", (ftnlen)30, (
			ftnlen)28);
	    }
/* Writing concatenation */
	    i__1[0] = 8, a__1[0] = "VARIABLE";
	    i__1[1] = 4, a__1[1] = chbufi;
	    i__1[2] = 30, a__1[2] = chbuf2;
	    s_cat(ch__1, a__1, i__1, &c__3, (ftnlen)42);
	    mnwarn_("W", mn7tit_1.cfrom, ch__1, (ftnlen)1, (ftnlen)8, (ftnlen)
		    42);
	} else {
	    *pinti = asin(yy);
	}
    }
    return 0;
} /* mnpint_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnplot_(doublereal *xpt, doublereal *ypt, char *chpt, 
	integer *nxypt, integer *nunit, integer *npagwd, integer *npagln, 
	ftnlen chpt_len)
{
    /* Initialized data */

    static char cdot[1+1] = ".";
    static char cslash[1+1] = "/";
    static char cblank[1+1] = " ";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, ni;
    static doublereal ax, ay, bx, by;
    static integer ks, ix, nx, ny, km1, ibk;
    static doublereal any, dxx, dyy;
    static integer isp1, iten;
    static doublereal xmin, ymin, xmax, ymax, savx, savy, yprt;
    static char cline[100], chsav[1];
    static doublereal bwidx, bwidy, xbest, ybest;
    static integer maxnx, maxny, iquit;
    static char chbest[1];
    static integer linodd;
    static char chmess[30];
    extern /* Subroutine */ int mnbins_(doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, doublereal *);
    static integer nxbest, nybest;
    static logical overpr;
    static doublereal xvalus[12];

    /* Fortran I/O blocks */
    static cilist io___890 = { 0, 0, 0, "(18X,A)", 0 };
    static cilist io___891 = { 0, 0, 0, "(1X,G14.7,A,A)", 0 };
    static cilist io___892 = { 0, 0, 0, "(18X,A)", 0 };
    static cilist io___895 = { 0, 0, 0, "(12X,12G10.4)", 0 };
    static cilist io___897 = { 0, 0, 0, "(25X,A,G13.7,A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        plots points in array xypt onto one page with labelled axes */
/* C        NXYPT is the number of points to be plotted */
/* C        XPT(I) = x-coord. of ith point */
/* C        YPT(I) = y-coord. of ith point */
/* C        CHPT(I) = character to be plotted at this position */
/* C        the input point arrays XPT, YPT, CHPT are destroyed. */
/* C */
    /* Parameter adjustments */
    --chpt;
    --ypt;
    --xpt;

    /* Function Body */
/* Computing MIN */
    i__1 = *npagwd - 20;
    maxnx = min(i__1,100);
    if (maxnx < 10) {
	maxnx = 10;
    }
    maxny = *npagln;
    if (maxny < 10) {
	maxny = 10;
    }
    if (*nxypt <= 1) {
	return 0;
    }
    xbest = xpt[1];
    ybest = ypt[1];
    *(unsigned char *)chbest = *(unsigned char *)&chpt[1];
/*         order the points by decreasing y */
    km1 = *nxypt - 1;
    i__1 = km1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iquit = 0;
	ni = *nxypt - i__;
	i__2 = ni;
	for (j = 1; j <= i__2; ++j) {
	    if (ypt[j] > ypt[j + 1]) {
		goto L140;
	    }
	    savx = xpt[j];
	    xpt[j] = xpt[j + 1];
	    xpt[j + 1] = savx;
	    savy = ypt[j];
	    ypt[j] = ypt[j + 1];
	    ypt[j + 1] = savy;
	    *(unsigned char *)chsav = *(unsigned char *)&chpt[j];
	    *(unsigned char *)&chpt[j] = *(unsigned char *)&chpt[j + 1];
	    *(unsigned char *)&chpt[j + 1] = *(unsigned char *)chsav;
	    iquit = 1;
L140:
	    ;
	}
	if (iquit == 0) {
	    goto L160;
	}
/* L150: */
    }
L160:
/*         find extreme values */
    xmax = xpt[1];
    xmin = xmax;
    i__1 = *nxypt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (xpt[i__] > xmax) {
	    xmax = xpt[i__];
	}
	if (xpt[i__] < xmin) {
	    xmin = xpt[i__];
	}
/* L200: */
    }
    dxx = (xmax - xmin) * (float).001;
    xmax += dxx;
    xmin -= dxx;
    mnbins_(&xmin, &xmax, &maxnx, &xmin, &xmax, &nx, &bwidx);
    ymax = ypt[1];
    ymin = ypt[*nxypt];
    if (ymax == ymin) {
	ymax = ymin + (float)1.;
    }
    dyy = (ymax - ymin) * (float).001;
    ymax += dyy;
    ymin -= dyy;
    mnbins_(&ymin, &ymax, &maxny, &ymin, &ymax, &ny, &bwidy);
    any = (doublereal) ny;
/*         if first point is blank, it is an 'origin' */
    if (*(unsigned char *)chbest == *(unsigned char *)&cblank[0]) {
	goto L50;
    }
    xbest = (xmax + xmin) * (float).5;
    ybest = (ymax + ymin) * (float).5;
L50:
/*         find scale constants */
    ax = (float)1. / bwidx;
    ay = (float)1. / bwidy;
    bx = -ax * xmin + (float)2.;
    by = -ay * ymin - (float)2.;
/*         convert points to grid positions */
    i__1 = *nxypt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xpt[i__] = ax * xpt[i__] + bx;
/* L300: */
	ypt[i__] = any - ay * ypt[i__] - by;
    }
    nxbest = (integer) (ax * xbest + bx);
    nybest = (integer) (any - ay * ybest - by);
/*         print the points */
    ny += 2;
    nx += 2;
    isp1 = 1;
    linodd = 1;
    overpr = FALSE_;
    i__1 = ny;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nx;
	for (ibk = 1; ibk <= i__2; ++ibk) {
/* L310: */
	    *(unsigned char *)&cline[ibk - 1] = *(unsigned char *)&cblank[0];
	}
	*(unsigned char *)cline = *(unsigned char *)&cdot[0];
	*(unsigned char *)&cline[nx - 1] = *(unsigned char *)&cdot[0];
	*(unsigned char *)&cline[nxbest - 1] = *(unsigned char *)&cdot[0];
	if (i__ != 1 && i__ != nybest && i__ != ny) {
	    goto L320;
	}
	i__2 = nx;
	for (j = 1; j <= i__2; ++j) {
/* L315: */
	    *(unsigned char *)&cline[j - 1] = *(unsigned char *)&cdot[0];
	}
L320:
	yprt = ymax - (real) (i__ - 1) * bwidy;
	if (isp1 > *nxypt) {
	    goto L350;
	}
/*         find the points to be plotted on this line */
	i__2 = *nxypt;
	for (k = isp1; k <= i__2; ++k) {
	    ks = (integer) ypt[k];
	    if (ks > i__) {
		goto L345;
	    }
	    ix = (integer) xpt[k];
	    if (*(unsigned char *)&cline[ix - 1] == *(unsigned char *)&cdot[0]
		    ) {
		goto L340;
	    }
	    if (*(unsigned char *)&cline[ix - 1] == *(unsigned char *)&cblank[
		    0]) {
		goto L340;
	    }
	    if (*(unsigned char *)&cline[ix - 1] == *(unsigned char *)&chpt[k]
		    ) {
		goto L341;
	    }
	    overpr = TRUE_;
/*         OVERPR is true if one or more positions contains more than */
/*            one point */
	    *(unsigned char *)&cline[ix - 1] = '&';
	    goto L341;
L340:
	    *(unsigned char *)&cline[ix - 1] = *(unsigned char *)&chpt[k];
L341:
	    ;
	}
	isp1 = *nxypt + 1;
	goto L350;
L345:
	isp1 = k;
L350:
	if (linodd == 1 || i__ == ny) {
	    goto L380;
	}
	linodd = 1;
	io___890.ciunit = *nunit;
	s_wsfe(&io___890);
	do_fio(&c__1, cline, nx);
	e_wsfe();
	goto L400;
L380:
	io___891.ciunit = *nunit;
	s_wsfe(&io___891);
	do_fio(&c__1, (char *)&yprt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " ..", (ftnlen)3);
	do_fio(&c__1, cline, nx);
	e_wsfe();
	linodd = 0;
L400:
	;
    }
/*         print labels on x-axis every ten columns */
    i__1 = nx;
    for (ibk = 1; ibk <= i__1; ++ibk) {
	*(unsigned char *)&cline[ibk - 1] = *(unsigned char *)&cblank[0];
	if (ibk % 10 == 1) {
	    *(unsigned char *)&cline[ibk - 1] = *(unsigned char *)&cslash[0];
	}
/* L410: */
    }
    io___892.ciunit = *nunit;
    s_wsfe(&io___892);
    do_fio(&c__1, cline, nx);
    e_wsfe();

    for (ibk = 1; ibk <= 12; ++ibk) {
/* L430: */
	xvalus[ibk - 1] = xmin + (real) (ibk - 1) * (float)10. * bwidx;
    }
    iten = (nx + 9) / 10;
    io___895.ciunit = *nunit;
    s_wsfe(&io___895);
    i__1 = iten;
    for (ibk = 1; ibk <= i__1; ++ibk) {
	do_fio(&c__1, (char *)&xvalus[ibk - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_copy(chmess, " ", (ftnlen)30, (ftnlen)1);
    if (overpr) {
	s_copy(chmess, "   Overprint character is &", (ftnlen)30, (ftnlen)27);
    }
    io___897.ciunit = *nunit;
    s_wsfe(&io___897);
    do_fio(&c__1, "ONE COLUMN=", (ftnlen)11);
    do_fio(&c__1, (char *)&bwidx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, chmess, (ftnlen)30);
    e_wsfe();
/* L500: */
    return 0;
} /* mnplot_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnpout_(integer *iuext, char *chnam, doublereal *val, 
	doublereal *err, doublereal *xlolim, doublereal *xuplim, integer *
	iuint, ftnlen chnam_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer nvl, iint, iext;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C     User-called */
/* C   Provides the user with information concerning the current status */
/* C          of parameter number IUEXT. Namely, it returns: */
/* C        CHNAM: the name of the parameter */
/* C        VAL: the current (external) value of the parameter */
/* C        ERR: the current estimate of the parameter uncertainty */
/* C        XLOLIM: the lower bound (or zero if no limits) */
/* C        XUPLIM: the upper bound (or zero if no limits) */
/* C        IUINT: the internal parameter number (or zero if not variable, */
/* C           or negative if undefined). */
/* C  Note also:  If IUEXT is negative, then it is -internal parameter */
/* C           number, and IUINT is returned as the EXTERNAL number. */
/* C     Except for IUINT, this is exactly the inverse of MNPARM */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    *xlolim = (float)0.;
    *xuplim = (float)0.;
    *err = (float)0.;
    if (*iuext == 0) {
	goto L100;
    }
    if (*iuext < 0) {
/*                   internal parameter number specified */
	iint = -(*iuext);
	if (iint > mn7npr_1.npar) {
	    goto L100;
	}
	iext = mn7inx_1.nexofi[iint - 1];
	*iuint = iext;
    } else {
/*                    external parameter number specified */
	iext = *iuext;
	if (iext == 0) {
	    goto L100;
	}
	if (iext > mn7npr_1.nu) {
	    goto L100;
	}
	iint = mn7inx_1.niofex[iext - 1];
	*iuint = iint;
    }
/*                     in both cases */
    nvl = mn7inx_1.nvarl[iext - 1];
    if (nvl < 0) {
	goto L100;
    }
    s_copy(chnam, mn7nam_1.cpnam + (iext - 1) * 10, chnam_len, (ftnlen)10);
    *val = mn7ext_1.u[iext - 1];
    if (iint > 0) {
	*err = mn7err_1.werr[iint - 1];
    }
    if (nvl == 4) {
	*xlolim = mn7ext_1.alim[iext - 1];
	*xuplim = mn7ext_1.blim[iext - 1];
    }
    return 0;
/*                parameter is undefined */
L100:
    *iuint = -1;
    s_copy(chnam, "undefined", chnam_len, (ftnlen)9);
    *val = (float)0.;
    return 0;
} /* mnpout_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnprin_(integer *inkode, doublereal *fval)
{
    /* Initialized data */

    static char cblank[11+1] = "           ";

    /* Format strings */
    static char fmt_905[] = "(/\002 FCN=\002,a,\002 FROM \002,a8,\002  STATU\
S=\002,a10,i6,\002 CALLS\002,i9,\002 TOTAL\002)";
    static char fmt_907[] = "(21x,\002EDM=\002,a,\002    STRATEGY=\002,i2,6x\
,a)";
    static char fmt_908[] = "(21x,\002EDM=\002,a,\002  STRATEGY=\002,i1,\002\
  ERROR MATRIX\002,\002 UNCERTAINTY=\002,f5.1,\002%\002)";
    static char fmt_910[] = "(/\002  EXT PARAMETER \002,13x,6a14)";
    static char fmt_911[] = "(\002  NO.   NAME    \002,\002    VALUE    \002\
,6a14)";
    static char fmt_952[] = "(i4,1x,a11,2g14.5,2a)";
    static char fmt_1004[] = "(\002 \002,32x,\002WARNING -   - ABOVE PARAMET\
ER IS AT LIMIT.\002)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfi(icilist *), e_wsfi()
	    ;
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double cos(doublereal);

    /* Local variables */
    static integer i__, k, l, m;
    static doublereal x1, x2, x3, dc;
    static integer ic, nc, kk;
    static char cx2[14], cx3[14];
    static integer lbl, nadd, ncol;
    static char chedm[10];
    static integer ikode;
    static doublereal dcmax;
    static char cnambf[11], colhdl[14*6], cheval[15], colhdu[14*6];
    static integer ntrail;

    /* Fortran I/O blocks */
    static cilist io___902 = { 0, 0, 0, "(A)", 0 };
    static cilist io___907 = { 0, 0, 0, "(/A,A)", 0 };
    static icilist io___909 = { 0, cheval, 0, "(G15.7)", 15, 1 };
    static icilist io___911 = { 0, chedm, 0, "(E10.2)", 10, 1 };
    static cilist io___913 = { 0, 0, 0, fmt_905, 0 };
    static cilist io___915 = { 0, 0, 0, fmt_907, 0 };
    static cilist io___918 = { 0, 0, 0, fmt_908, 0 };
    static cilist io___925 = { 0, 0, 0, fmt_910, 0 };
    static cilist io___927 = { 0, 0, 0, fmt_911, 0 };
    static cilist io___933 = { 0, 0, 0, fmt_952, 0 };
    static icilist io___936 = { 0, cx2, 0, "(G14.5)", 14, 1 };
    static icilist io___937 = { 0, cx3, 0, "(G14.5)", 14, 1 };
    static cilist io___938 = { 0, 0, 0, fmt_952, 0 };
    static cilist io___939 = { 0, 0, 0, fmt_1004, 0 };
    static cilist io___940 = { 0, 0, 0, "(I4,1X,A11,G14.5,A,2G14.5)", 0 };
    static cilist io___941 = { 0, 0, 0, "(I4,1X,A11,G14.5,A)", 0 };
    static cilist io___942 = { 0, 0, 0, "(31X,A,G10.3)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Prints the values of the parameters at the time of the call. */
/* C        also prints other relevant information such as function value, */
/* C        estimated distance to minimum, parameter errors, step sizes. */
/* C */
/*         According to the value of IKODE, the printout is: */
/*    IKODE=INKODE= 0    only info about function value */
/*                  1    parameter values, errors, limits */
/*                  2    values, errors, step sizes, internal values */
/*                  3    values, errors, step sizes, first derivs. */
/*                  4    values, parabolic errors, MINOS errors */
/*    when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to ISW(2) */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */




    if (mn7npr_1.nu == 0) {
	io___902.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___902);
	do_fio(&c__1, " THERE ARE CURRENTLY NO PARAMETERS DEFINED", (ftnlen)
		42);
	e_wsfe();
	goto L700;
    }
/*                  get value of IKODE based in INKODE, ISW(2) */
    ikode = *inkode;
    if (*inkode == 5) {
	ikode = mn7flg_1.isw[1] + 1;
	if (ikode > 3) {
	    ikode = 3;
	}
    }
/*                  set 'default' column headings */
    for (k = 1; k <= 6; ++k) {
	s_copy(colhdu + (k - 1) * 14, "UNDEFINED", (ftnlen)14, (ftnlen)9);
/* L5: */
	s_copy(colhdl + (k - 1) * 14, "COLUMN HEAD", (ftnlen)14, (ftnlen)11);
    }
/*              print title if Minos errors, and title exists. */
    if (ikode == 4 && s_cmp(mn7tit_1.ctitl, mn7tit_1.cundef, (ftnlen)50, (
	    ftnlen)10) != 0) {
	io___907.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___907);
	do_fio(&c__1, " MINUIT TASK: ", (ftnlen)14);
	do_fio(&c__1, mn7tit_1.ctitl, (ftnlen)50);
	e_wsfe();
    }
/*              report function value and status */
    if (*fval == mn7cns_1.undefi) {
	s_copy(cheval, " unknown       ", (ftnlen)15, (ftnlen)15);
    } else {
	s_wsfi(&io___909);
	do_fio(&c__1, (char *)&(*fval), (ftnlen)sizeof(doublereal));
	e_wsfi();
    }
    if (mn7min_1.edm == mn7cns_1.bigedm) {
	s_copy(chedm, " unknown  ", (ftnlen)10, (ftnlen)10);
    } else {
	s_wsfi(&io___911);
	do_fio(&c__1, (char *)&mn7min_1.edm, (ftnlen)sizeof(doublereal));
	e_wsfi();
    }
    nc = mn7cnv_1.nfcn - mn7cnv_1.nfcnfr;
    io___913.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___913);
    do_fio(&c__1, cheval, (ftnlen)15);
    do_fio(&c__1, mn7tit_1.cfrom, (ftnlen)8);
    do_fio(&c__1, mn7tit_1.cstatu, (ftnlen)10);
    do_fio(&c__1, (char *)&nc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&mn7cnv_1.nfcn, (ftnlen)sizeof(integer));
    e_wsfe();
    m = mn7flg_1.isw[1];
    if (m == 0 || m == 2 || mn7min_1.dcovar == 0.) {
	io___915.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___915);
	do_fio(&c__1, chedm, (ftnlen)10);
	do_fio(&c__1, (char *)&mn7cnv_1.istrat, (ftnlen)sizeof(integer));
	do_fio(&c__1, mn7tit_1.covmes + m * 22, (ftnlen)22);
	e_wsfe();
    } else {
	dcmax = (float)1.;
	dc = min(mn7min_1.dcovar,dcmax) * (float)100.;
	io___918.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___918);
	do_fio(&c__1, chedm, (ftnlen)10);
	do_fio(&c__1, (char *)&mn7cnv_1.istrat, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&dc, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

    if (ikode == 0) {
	goto L700;
    }
/*               find longest name (for Rene!) */
    ntrail = 10;
    i__1 = mn7npr_1.nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (mn7inx_1.nvarl[i__ - 1] < 0) {
	    goto L20;
	}
	for (ic = 10; ic >= 1; --ic) {
	    if (*(unsigned char *)&mn7nam_1.cpnam[(i__ - 1) * 10 + (ic - 1)] 
		    != ' ') {
		goto L16;
	    }
/* L15: */
	}
	ic = 1;
L16:
	lbl = 10 - ic;
	if (lbl < ntrail) {
	    ntrail = lbl;
	}
L20:
	;
    }
    nadd = ntrail / 2 + 1;
    if (ikode == 1) {
	s_copy(colhdu, "              ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl, "      ERROR   ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdu + 14, "      PHYSICAL", (ftnlen)14, (ftnlen)14);
	s_copy(colhdu + 28, " LIMITS       ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl + 14, "    NEGATIVE  ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl + 28, "    POSITIVE  ", (ftnlen)14, (ftnlen)14);
    }
    if (ikode == 2) {
	s_copy(colhdu, "              ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl, "      ERROR   ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdu + 14, "    INTERNAL  ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl + 14, "    STEP SIZE ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdu + 28, "    INTERNAL  ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl + 28, "      VALUE   ", (ftnlen)14, (ftnlen)14);
    }
    if (ikode == 3) {
	s_copy(colhdu, "              ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl, "      ERROR   ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdu + 14, "       STEP   ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl + 14, "       SIZE   ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdu + 28, "      FIRST   ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl + 28, "   DERIVATIVE ", (ftnlen)14, (ftnlen)14);
    }
    if (ikode == 4) {
	s_copy(colhdu, "    PARABOLIC ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl, "      ERROR   ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdu + 14, "        MINOS ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdu + 28, "ERRORS        ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl + 14, "   NEGATIVE   ", (ftnlen)14, (ftnlen)14);
	s_copy(colhdl + 28, "   POSITIVE   ", (ftnlen)14, (ftnlen)14);
    }

    if (ikode != 4) {
	if (mn7flg_1.isw[1] < 3) {
	    s_copy(colhdu, "  APPROXIMATE ", (ftnlen)14, (ftnlen)14);
	}
	if (mn7flg_1.isw[1] < 1) {
	    s_copy(colhdu, " CURRENT GUESS", (ftnlen)14, (ftnlen)14);
	}
    }
    ncol = 3;
    io___925.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___925);
    i__1 = ncol;
    for (kk = 1; kk <= i__1; ++kk) {
	do_fio(&c__1, colhdu + (kk - 1) * 14, (ftnlen)14);
    }
    e_wsfe();
    io___927.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___927);
    i__1 = ncol;
    for (kk = 1; kk <= i__1; ++kk) {
	do_fio(&c__1, colhdl + (kk - 1) * 14, (ftnlen)14);
    }
    e_wsfe();

/*                                        . . . loop over parameters . . */
    i__1 = mn7npr_1.nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (mn7inx_1.nvarl[i__ - 1] < 0) {
	    goto L200;
	}
	l = mn7inx_1.niofex[i__ - 1];
/* Writing concatenation */
	i__2[0] = nadd, a__1[0] = cblank;
	i__2[1] = 10, a__1[1] = mn7nam_1.cpnam + (i__ - 1) * 10;
	s_cat(cnambf, a__1, i__2, &c__2, (ftnlen)11);
	if (l == 0) {
	    goto L55;
	}
/*              variable parameter. */
	x1 = mn7err_1.werr[l - 1];
	s_copy(cx2, "PLEASE GET X..", (ftnlen)14, (ftnlen)14);
	s_copy(cx3, "PLEASE GET X..", (ftnlen)14, (ftnlen)14);
	if (ikode == 1) {
	    if (mn7inx_1.nvarl[i__ - 1] <= 1) {
		io___933.ciunit = mn7iou_1.isyswr;
		s_wsfe(&io___933);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, cnambf, (ftnlen)11);
		do_fio(&c__1, (char *)&mn7ext_1.u[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&x1, (ftnlen)sizeof(doublereal));
		e_wsfe();
		goto L200;
	    } else {
		x2 = mn7ext_1.alim[i__ - 1];
		x3 = mn7ext_1.blim[i__ - 1];
	    }
	}
	if (ikode == 2) {
	    x2 = mn7int_1.dirin[l - 1];
	    x3 = mn7int_1.x[l - 1];
	}
	if (ikode == 3) {
	    x2 = mn7int_1.dirin[l - 1];
	    x3 = mn7der_1.grd[l - 1];
	    if (mn7inx_1.nvarl[i__ - 1] > 1 && (d__1 = cos(mn7int_1.x[l - 1]),
		     abs(d__1)) < (float).001) {
		s_copy(cx3, "** at limit **", (ftnlen)14, (ftnlen)14);
	    }
	}
	if (ikode == 4) {
	    x2 = mn7err_1.ern[l - 1];
	    if (x2 == 0.) {
		s_copy(cx2, " ", (ftnlen)14, (ftnlen)1);
	    }
	    if (x2 == mn7cns_1.undefi) {
		s_copy(cx2, "   at limit   ", (ftnlen)14, (ftnlen)14);
	    }
	    x3 = mn7err_1.erp[l - 1];
	    if (x3 == 0.) {
		s_copy(cx3, " ", (ftnlen)14, (ftnlen)1);
	    }
	    if (x3 == mn7cns_1.undefi) {
		s_copy(cx3, "   at limit   ", (ftnlen)14, (ftnlen)14);
	    }
	}
	if (s_cmp(cx2, "PLEASE GET X..", (ftnlen)14, (ftnlen)14) == 0) {
	    s_wsfi(&io___936);
	    do_fio(&c__1, (char *)&x2, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	}
	if (s_cmp(cx3, "PLEASE GET X..", (ftnlen)14, (ftnlen)14) == 0) {
	    s_wsfi(&io___937);
	    do_fio(&c__1, (char *)&x3, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	}
	io___938.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___938);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, cnambf, (ftnlen)11);
	do_fio(&c__1, (char *)&mn7ext_1.u[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, (char *)&x1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, cx2, (ftnlen)14);
	do_fio(&c__1, cx3, (ftnlen)14);
	e_wsfe();
/*               check if parameter is at limit */
	if (mn7inx_1.nvarl[i__ - 1] <= 1 || ikode == 3) {
	    goto L200;
	}
	if ((d__1 = cos(mn7int_1.x[l - 1]), abs(d__1)) < (float).001) {
	    io___939.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___939);
	    e_wsfe();
	}
	goto L200;

/*                                print constant or fixed parameter. */
L55:
	s_copy(colhdu, "   constant   ", (ftnlen)14, (ftnlen)14);
	if (mn7inx_1.nvarl[i__ - 1] > 0) {
	    s_copy(colhdu, "     fixed    ", (ftnlen)14, (ftnlen)14);
	}
	if (mn7inx_1.nvarl[i__ - 1] == 4 && ikode == 1) {
	    io___940.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___940);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, cnambf, (ftnlen)11);
	    do_fio(&c__1, (char *)&mn7ext_1.u[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, colhdu, (ftnlen)14);
	    do_fio(&c__1, (char *)&mn7ext_1.alim[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&mn7ext_1.blim[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else {
	    io___941.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___941);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, cnambf, (ftnlen)11);
	    do_fio(&c__1, (char *)&mn7ext_1.u[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, colhdu, (ftnlen)14);
	    e_wsfe();
	}
L200:
	;
    }

    if (mn7min_1.up != mn7cns_1.updflt) {
	io___942.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___942);
	do_fio(&c__1, "ERR DEF=", (ftnlen)8);
	do_fio(&c__1, (char *)&mn7min_1.up, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
L700:
    return 0;
} /* mnprin_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:50  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnpsdf_()
{
    /* Format strings */
    static char fmt_550[] = "(\002 EIGENVALUES OF SECOND-DERIVATIVE MATRIX\
:\002)";
    static char fmt_551[] = "(7x,6e12.4)";

    /* System generated locals */
    address a__1[3], a__2[2];
    integer i__1, i__2[3], i__3[2], i__4;
    doublereal d__1;
    char ch__1[44], ch__2[46], ch__3[57];

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi();
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;
    static doublereal s[100], dg;
    static integer ip;
    static doublereal padd;
    static integer ndex;
    static doublereal pmin, pmax, dgmin;
    extern /* Subroutine */ int mneig_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *);
    static integer ndexd;
    static char chbuff[12];
    static doublereal epspdf;
    static integer ifault;
    static doublereal epsmin;
    extern /* Subroutine */ int mnwarn_(char *, char *, char *, ftnlen, 
	    ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___949 = { 0, chbuff, 0, "(I3)", 3, 1 };
    static icilist io___951 = { 0, chbuff, 0, "(E12.2)", 12, 1 };
    static cilist io___959 = { 0, 0, 0, fmt_550, 0 };
    static cilist io___960 = { 0, 0, 0, fmt_551, 0 };
    static icilist io___962 = { 0, chbuff, 0, "(G12.5)", 12, 1 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        calculates the eigenvalues of v to see if positive-def. */
/* C        if not, adds constant along diagonal to make positive. */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    epsmin = (float)1e-6;
    epspdf = max(epsmin,mn7cns_1.epsma2);
    dgmin = mn7var_1.vhmat[0];
/*                        Check if negative or zero on diagonal */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ndex = i__ * (i__ + 1) / 2;
	if (mn7var_1.vhmat[ndex - 1] <= 0.) {
	    s_wsfi(&io___949);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__2[0] = 25, a__1[0] = "Negative diagonal element";
	    i__2[1] = 3, a__1[1] = chbuff;
	    i__2[2] = 16, a__1[2] = " in Error Matrix";
	    s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)44);
	    mnwarn_("W", mn7tit_1.cfrom, ch__1, (ftnlen)1, (ftnlen)8, (ftnlen)
		    44);
	}
	if (mn7var_1.vhmat[ndex - 1] < dgmin) {
	    dgmin = mn7var_1.vhmat[ndex - 1];
	}
/* L200: */
    }
    if (dgmin <= 0.) {
	dg = epspdf + 1. - dgmin;
	s_wsfi(&io___951);
	do_fio(&c__1, (char *)&dg, (ftnlen)sizeof(doublereal));
	e_wsfi();
/* Writing concatenation */
	i__3[0] = 12, a__2[0] = chbuff;
	i__3[1] = 34, a__2[1] = " added to diagonal of error matrix";
	s_cat(ch__2, a__2, i__3, &c__2, (ftnlen)46);
	mnwarn_("W", mn7tit_1.cfrom, ch__2, (ftnlen)1, (ftnlen)8, (ftnlen)46);
    } else {
	dg = 0.;
    }
/*                    Store VHMAT in P, make sure diagonal pos. */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ndex = i__ * (i__ - 1) / 2;
	ndexd = ndex + i__;
	mn7var_1.vhmat[ndexd - 1] += dg;
	if (mn7var_1.vhmat[ndexd - 1] <= 0.) {
	    mn7var_1.vhmat[ndexd - 1] = (float)1.;
	}
	s[i__ - 1] = (float)1. / sqrt(mn7var_1.vhmat[ndexd - 1]);
	i__4 = i__;
	for (j = 1; j <= i__4; ++j) {
	    ++ndex;
/* L213: */
	    mn7sim_1.p[i__ + j * 100 - 101] = mn7var_1.vhmat[ndex - 1] * s[
		    i__ - 1] * s[j - 1];
	}
    }
/*      call eigen (p,p,maxint,npar,pstar,-npar) */
    mneig_(mn7sim_1.p, &mn7npr_1.maxint, &mn7npr_1.npar, &mn7npr_1.maxint, 
	    mn7sim_1.pstar, &epspdf, &ifault);
    pmin = mn7sim_1.pstar[0];
    pmax = mn7sim_1.pstar[0];
    i__4 = mn7npr_1.npar;
    for (ip = 2; ip <= i__4; ++ip) {
	if (mn7sim_1.pstar[ip - 1] < pmin) {
	    pmin = mn7sim_1.pstar[ip - 1];
	}
	if (mn7sim_1.pstar[ip - 1] > pmax) {
	    pmax = mn7sim_1.pstar[ip - 1];
	}
/* L215: */
    }
/* Computing MAX */
    d__1 = abs(pmax);
    pmax = max(d__1,1.);
    if (pmin <= 0. && mn7log_1.lwarn || mn7flg_1.isw[4] >= 2) {
	io___959.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___959);
	e_wsfe();
	io___960.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___960);
	i__4 = mn7npr_1.npar;
	for (ip = 1; ip <= i__4; ++ip) {
	    do_fio(&c__1, (char *)&mn7sim_1.pstar[ip - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
    }
    if (pmin > epspdf * pmax) {
	goto L217;
    }
    if (mn7flg_1.isw[1] == 3) {
	mn7flg_1.isw[1] = 2;
    }
    padd = pmax * (float).001 - pmin;
    i__4 = mn7npr_1.npar;
    for (ip = 1; ip <= i__4; ++ip) {
	ndex = ip * (ip + 1) / 2;
/* L216: */
	mn7var_1.vhmat[ndex - 1] *= padd + (float)1.;
    }
    s_copy(mn7tit_1.cstatu, "NOT POSDEF", (ftnlen)10, (ftnlen)10);
    s_wsfi(&io___962);
    do_fio(&c__1, (char *)&padd, (ftnlen)sizeof(doublereal));
    e_wsfi();
/* Writing concatenation */
    i__2[0] = 32, a__1[0] = "MATRIX FORCED POS-DEF BY ADDING ";
    i__2[1] = 12, a__1[1] = chbuff;
    i__2[2] = 13, a__1[2] = " TO DIAGONAL.";
    s_cat(ch__3, a__1, i__2, &c__3, (ftnlen)57);
    mnwarn_("W", mn7tit_1.cfrom, ch__3, (ftnlen)1, (ftnlen)8, (ftnlen)57);
L217:

    return 0;
} /* mnpsdf_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnrazz_(doublereal *ynew, doublereal *pnew, doublereal *
	y, integer *jh, integer *jl)
{
    /* Format strings */
    static char fmt_1000[] = "(\002   FUNCTION VALUE DOES NOT SEEM TO DEPEND\
 ON ANY OF THE\002,i3,\002 VARIABLE PARAMETERS.\002/10x,\002VERIFY THAT STEP\
 SIZES ARE\002,\002 BIG ENOUGH AND CHECK FCN LOGIC.\002/1x,79(\002*\002)/1x,\
79(\002*\002)/)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer i__, j;
    static doublereal pbig, plit;
    static integer nparp1;
    extern /* Subroutine */ int mninex_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___968 = { 0, 0, 0, fmt_1000, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called only by MNSIMP (and MNIMPR) to add a new point */
/* C        and remove an old one from the current simplex, and get the */
/* C        estimated distance to minimum. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    /* Parameter adjustments */
    --y;
    --pnew;

    /* Function Body */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	mn7sim_1.p[i__ + *jh * 100 - 101] = pnew[i__];
    }
    y[*jh] = *ynew;
    if (*ynew < mn7min_1.amin) {
	i__1 = mn7npr_1.npar;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	    mn7int_1.x[i__ - 1] = pnew[i__];
	}
	mninex_(mn7int_1.x);
	mn7min_1.amin = *ynew;
	s_copy(mn7tit_1.cstatu, "PROGRESS  ", (ftnlen)10, (ftnlen)10);
	*jl = *jh;
    }
    *jh = 1;
    nparp1 = mn7npr_1.npar + 1;
/* L20: */
    i__1 = nparp1;
    for (j = 2; j <= i__1; ++j) {
	if (y[j] > y[*jh]) {
	    *jh = j;
	}
/* L25: */
    }
    mn7min_1.edm = y[*jh] - y[*jl];
    if (mn7min_1.edm <= 0.) {
	goto L45;
    }
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pbig = mn7sim_1.p[i__ - 1];
	plit = pbig;
	i__2 = nparp1;
	for (j = 2; j <= i__2; ++j) {
	    if (mn7sim_1.p[i__ + j * 100 - 101] > pbig) {
		pbig = mn7sim_1.p[i__ + j * 100 - 101];
	    }
	    if (mn7sim_1.p[i__ + j * 100 - 101] < plit) {
		plit = mn7sim_1.p[i__ + j * 100 - 101];
	    }
/* L30: */
	}
	mn7int_1.dirin[i__ - 1] = pbig - plit;
/* L35: */
    }
L40:
    return 0;
L45:
    io___968.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___968);
    do_fio(&c__1, (char *)&mn7npr_1.npar, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L40;
} /* mnrazz_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnread_(S_fp fcn, integer *iflgin, integer *iflgut, U_fp 
	futil)
{
    /* Initialized data */

    static char cpromt[40*3+1] = " ENTER MINUIT TITLE, or \"SET INPUT n\" : \
 ENTER MINUIT PARAMETER DEFINITION:      ENTER MINUIT COMMAND:              \
    ";
    static char clower[26+1] = "abcdefghijklmnopqrstuvwxyz";
    static char cupper[26+1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    /* Format strings */
    static char fmt_21[] = "(\002 **********\002/\002 **\002,i5,\002 **\002,\
a/\002 **********\002)";
    static char fmt_405[] = "(\002 \002,10(\002*\002)/\002 **\002,i5,\002 *\
*\002,a)";
    static char fmt_420[] = "(bn,7e11.4,3x)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsfe(cilist *), e_rsfe(), i_indx(char *, char *, ftnlen, ftnlen)
	    ;

    /* Local variables */
    static integer i__, ic;
    static logical leof;
    static integer ierr, npar2;
    static char crdbuf[80];
    static integer iflgdo, icondp;
    static char cupbuf[10];
    static integer incomp;
    extern /* Subroutine */ int mnstin_(char *, integer *, ftnlen), mnseti_(
	    char *, ftnlen), mnpars_(char *, integer *, ftnlen), mncomd_(S_fp,
	     char *, integer *, U_fp, ftnlen);
    static integer icondn;
    extern /* Subroutine */ int mnmatu_(integer *), mnprin_(integer *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___975 = { 0, 0, 0, "(A)", 0 };
    static cilist io___977 = { 1, 0, 1, "(A)", 0 };
    static cilist io___981 = { 0, 0, 0, "(A,I3)", 0 };
    static cilist io___982 = { 0, 0, 0, fmt_21, 0 };
    static cilist io___983 = { 0, 0, 0, "(A,I3)", 0 };
    static cilist io___985 = { 0, 0, 0, "(A,A/)", 0 };
    static cilist io___986 = { 0, 0, 0, "(1X,A50)", 0 };
    static cilist io___987 = { 0, 0, 0, "(1X,78(1H*))", 0 };
    static cilist io___989 = { 0, 0, 0, "(A)", 0 };
    static cilist io___990 = { 0, 0, 0, "(A)", 0 };
    static cilist io___991 = { 0, 0, 0, "(4X,75(1H*))", 0 };
    static cilist io___993 = { 0, 0, 0, fmt_405, 0 };
    static cilist io___994 = { 0, 0, 0, "(1H ,10(1H*))", 0 };
    static cilist io___996 = { 1, 0, 1, fmt_420, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called from MINUIT.  Reads all user input to MINUIT. */
/* C     This routine is highly unstructured and defies normal logic. */
/* C */
/* C     IFLGIN indicates the function originally requested: */
/* C           = 1: read one-line title */
/* C             2: read parameter definitions */
/* C             3: read MINUIT commands */
/* C */
/* C     IFLGUT= 1: reading terminated normally */
/* C             2: end-of-data on input */
/* C             3: unrecoverable read error */
/* C             4: unable to process parameter requests */
/* C             5: more than 100 incomprehensible commands */
/* C internally, */
/* C     IFLGDO indicates the subfunction to be performed on the next */
/* C         input record: 1: read a one-line title */
/* C                       2: read a parameter definition */
/* C                       3: read a command */
/* C                       4: read in covariance matrix */
/* C     for example, when IFLGIN=3, but IFLGDO=1, then it should read */
/* C       a title, but this was requested by a command, not by MINUIT. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */




    *iflgut = 1;
    iflgdo = *iflgin;
    leof = FALSE_;
    incomp = 0;
/*                                           . . . . read next record */
L10:
    if (mn7flg_1.isw[5] == 1) {
	io___975.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___975);
	do_fio(&c__1, cpromt + (iflgdo - 1) * 40, (ftnlen)40);
	e_wsfe();
	if (iflgdo == 2) {
	    mn7log_1.lphead = FALSE_;
	}
    }
    s_copy(crdbuf, "   ", (ftnlen)80, (ftnlen)3);
    io___977.ciunit = mn7iou_1.isysrd;
    i__1 = s_rsfe(&io___977);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_fio(&c__1, crdbuf, (ftnlen)80);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = e_rsfe();
L100001:
    if (i__1 < 0) {
	goto L45;
    }
    if (i__1 > 0) {
	goto L500;
    }

/*                 CUPBUF is the first few characters in upper case */
    s_copy(cupbuf, crdbuf, (ftnlen)10, (ftnlen)10);
    for (i__ = 1; i__ <= 10; ++i__) {
	if (*(unsigned char *)&crdbuf[i__ - 1] == '\'') {
	    goto L13;
	}
	for (ic = 1; ic <= 26; ++ic) {
	    if (*(unsigned char *)&crdbuf[i__ - 1] == *(unsigned char *)&
		    clower[ic - 1]) {
		*(unsigned char *)&cupbuf[i__ - 1] = *(unsigned char *)&
			cupper[ic - 1];
	    }
/* L11: */
	}
/* L12: */
    }
L13:
/*                                           . .   preemptive commands */
    leof = FALSE_;
    if (i_indx(cupbuf, "*EOF", (ftnlen)10, (ftnlen)4) == 1) {
	io___981.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___981);
	do_fio(&c__1, " *EOF ENCOUNTERED ON UNIT NO.", (ftnlen)29);
	do_fio(&c__1, (char *)&mn7iou_1.isysrd, (ftnlen)sizeof(integer));
	e_wsfe();
	mn7log_1.lphead = TRUE_;
	goto L50;
    }
    if (i_indx(cupbuf, "SET INP", (ftnlen)10, (ftnlen)7) == 1) {
	++mn7flg_1.icomnd;
	io___982.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___982);
	do_fio(&c__1, (char *)&mn7flg_1.icomnd, (ftnlen)sizeof(integer));
	do_fio(&c__1, crdbuf, (ftnlen)50);
	e_wsfe();
	mn7log_1.lphead = TRUE_;
	goto L50;
    }
    goto L80;
/*                                    . . hardware EOF on current ISYSRD */
L45:
    s_copy(crdbuf, "*EOF ", (ftnlen)80, (ftnlen)5);
    io___983.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___983);
    do_fio(&c__1, " END OF DATA ON UNIT NO.", (ftnlen)24);
    do_fio(&c__1, (char *)&mn7iou_1.isysrd, (ftnlen)sizeof(integer));
    e_wsfe();
/*                                     or SET INPUT command */
L50:
    mnstin_(crdbuf, &ierr, (ftnlen)80);
    if (ierr == 0) {
	goto L10;
    }
    if (ierr == 2) {
	if (! leof) {
	    io___985.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___985);
	    do_fio(&c__1, " TWO CONSECUTIVE EOFs ON ", (ftnlen)25);
	    do_fio(&c__1, "PRIMARY INPUT FILE WILL TERMINATE EXECUTION.", (
		    ftnlen)44);
	    e_wsfe();
	    leof = TRUE_;
	    goto L10;
	}
    }
    *iflgut = ierr;
    goto L900;
L80:
    if (iflgdo > 1) {
	goto L100;
    }
/*                            read title        . . . . .   IFLGDO = 1 */
/*              if title is 'SET TITLE', skip and read again */
    if (i_indx(cupbuf, "SET TIT", (ftnlen)10, (ftnlen)7) == 1) {
	goto L10;
    }
    mnseti_(crdbuf, (ftnlen)50);
    io___986.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___986);
    do_fio(&c__1, mn7tit_1.ctitl, (ftnlen)50);
    e_wsfe();
    io___987.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___987);
    e_wsfe();
    mn7log_1.lphead = TRUE_;
    if (*iflgin == iflgdo) {
	goto L900;
    }
    iflgdo = *iflgin;
    goto L10;
/*                            data record is not a title. */
L100:
    if (iflgdo > 2) {
	goto L300;
    }
/*                          expect parameter definitions.   IFLGDO = 2 */
/*              if parameter def is 'PARAMETER', skip and read again */
    if (i_indx(cupbuf, "PAR", (ftnlen)10, (ftnlen)3) == 1) {
	goto L10;
    }
/*              if line starts with SET TITLE, read a title first */
    if (i_indx(cupbuf, "SET TIT", (ftnlen)10, (ftnlen)7) == 1) {
	iflgdo = 1;
	goto L10;
    }
/*                      we really have parameter definitions now */
    mnpars_(crdbuf, &icondp, (ftnlen)80);
    if (icondp == 0) {
	goto L10;
    }
/*          format error */
    if (icondp == 1) {
	if (mn7flg_1.isw[5] == 1) {
	    io___989.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___989);
	    do_fio(&c__1, " FORMAT ERROR.  IGNORED.  ENTER AGAIN.", (ftnlen)
		    38);
	    e_wsfe();
	    goto L10;
	} else {
	    io___990.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___990);
	    do_fio(&c__1, " ERROR IN PARAMETER DEFINITION", (ftnlen)30);
	    e_wsfe();
	    *iflgut = 4;
	    goto L900;
	}
    }
/*                     ICONDP = 2            . . . end parameter requests */
    if (mn7flg_1.isw[4] >= 0 && mn7flg_1.isw[5] < 1) {
	io___991.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___991);
	e_wsfe();
    }
    mn7log_1.lphead = TRUE_;
    if (*iflgin == iflgdo) {
	goto L900;
    }
    iflgdo = *iflgin;
    goto L10;
/*                                              . . . . .   IFLGDO = 3 */
/*                                           read commands */
L300:
    mncomd_((S_fp)fcn, crdbuf, &icondn, (U_fp)futil, (ftnlen)80);
/* C     ICONDN = 0: command executed normally */
/* C              1: command is blank, ignored */
/* C              2: command line unreadable, ignored */
/* C              3: unknown command, ignored */
/* C              4: abnormal termination (e.g., MIGRAD not converged) */
/* C              5: command is a request to read PARAMETER definitions */
/* C              6: 'SET INPUT' command */
/* C              7: 'SET TITLE' command */
/* C              8: 'SET COVAR' command */
/* C              9: reserved */
/* C             10: END command */
/* C             11: EXIT or STOP command */
/* C             12: RETURN command */
    if (icondn == 2 || icondn == 3) {
	++incomp;
	if (incomp > 100) {
	    *iflgut = 5;
	    goto L900;
	}
    }
/*                         parameter */
    if (icondn == 5) {
	iflgdo = 2;
    }
/*                         SET INPUT */
    if (icondn == 6) {
	goto L50;
    }
/*                         SET TITLE */
    if (icondn == 7) {
	iflgdo = 1;
    }
/*                                        . . . . . . . . . . set covar */
    if (icondn == 8) {
	++mn7flg_1.icomnd;
	io___993.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___993);
	do_fio(&c__1, (char *)&mn7flg_1.icomnd, (ftnlen)sizeof(integer));
	do_fio(&c__1, crdbuf, (ftnlen)50);
	e_wsfe();
	io___994.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___994);
	e_wsfe();
	npar2 = mn7npr_1.npar * (mn7npr_1.npar + 1) / 2;
	io___996.ciunit = mn7iou_1.isysrd;
	i__1 = s_rsfe(&io___996);
	if (i__1 != 0) {
	    goto L100002;
	}
	i__2 = npar2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = do_fio(&c__1, (char *)&mn7var_1.vhmat[i__ - 1], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L100002;
	    }
	}
	i__1 = e_rsfe();
L100002:
	if (i__1 < 0) {
	    goto L45;
	}
	if (i__1 > 0) {
	    goto L500;
	}
	mn7flg_1.isw[1] = 3;
	mn7min_1.dcovar = (float)0.;
	if (mn7flg_1.isw[4] >= 0) {
	    mnmatu_(&c__1);
	}
	if (mn7flg_1.isw[4] >= 1) {
	    mnprin_(&c__2, &mn7min_1.amin);
	}
	goto L10;
    }
    if (icondn < 10) {
	goto L10;
    }
    goto L900;
/*                                              . . . . error conditions */
L500:
    *iflgut = 3;
L900:
    return 0;
} /* mnread_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnrn15_(doublereal *val, integer *inseed)
{
    /* Initialized data */

    static integer iseed = 12345;

    static integer k;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/*         This is a super-portable random number generator. */
/*         It should not overflow on any 32-bit machine. */
/*         The cycle is only ~10**9, so use with care! */
/*         Note especially that VAL must not be undefined on input. */
/*                    Set Default Starting Seed */
    if (*val == 3.) {
	goto L100;
    }

    *inseed = iseed;
    k = iseed / 53668;
    iseed = (iseed - k * 53668) * 40014 - k * 12211;
    if (iseed < 0) {
	iseed += 2147483563;
    }
    *val = (real) iseed * (float)4.656613e-10;
    return 0;
/*               "entry" to set seed, flag is VAL=3. */
L100:
    iseed = *inseed;
    return 0;
} /* mnrn15_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnrset_(integer *iopt)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, iext;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called from MNCLER and whenever problem changes, for example */
/* C        after SET LIMITS, SET PARAM, CALL FCN 6 */
/* C    If IOPT=1, */
/* C        Resets function value and errors to UNDEFINED */
/* C    If IOPT=0, sets only MINOS errors to undefined */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    s_copy(mn7tit_1.cstatu, "RESET     ", (ftnlen)10, (ftnlen)10);
    if (*iopt >= 1) {
	mn7min_1.amin = mn7cns_1.undefi;
	mn7min_1.fval3 = abs(mn7min_1.amin) * (float)2. + (float)1.;
	mn7min_1.edm = mn7cns_1.bigedm;
	mn7flg_1.isw[3] = 0;
	mn7flg_1.isw[1] = 0;
	mn7min_1.dcovar = (float)1.;
	mn7flg_1.isw[0] = 0;
    }
    mn7log_1.lnolim = TRUE_;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iext = mn7inx_1.nexofi[i__ - 1];
	if (mn7inx_1.nvarl[iext - 1] >= 4) {
	    mn7log_1.lnolim = FALSE_;
	}
	mn7err_1.erp[i__ - 1] = 0.;
	mn7err_1.ern[i__ - 1] = 0.;
	mn7err_1.globcc[i__ - 1] = 0.;
/* L10: */
    }
    if (mn7flg_1.isw[1] >= 1) {
	mn7flg_1.isw[1] = 1;
	mn7min_1.dcovar = max(mn7min_1.dcovar,.5);
    }
    return 0;
} /* mnrset_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnsave_()
{
    /* Format strings */
    static char fmt_32[] = "(\002 CURRENT VALUES WILL BE SAVED ON UNIT\002,i\
3,\002: \002,a/)";
    static char fmt_35[] = "(\002 UNIT\002,i3,\002 IS NOT OPENED.\002)";
    static char fmt_37[] = "(\002 SHOULD UNIT\002,i3,\002 BE REWOUND BEFORE \
WRITING TO IT?\002)";
    static char fmt_1001[] = "(1x,i5,\002'\002,a10,\002'\002,4e13.5)";
    static char fmt_1003[] = "(\002SET COVARIANCE\002,i6)";
    static char fmt_1004[] = "(bn,7e11.4,3x)";
    static char fmt_501[] = "(1x,i5,\002 RECORDS WRITTEN TO UNIT\002,i4\
,\002:\002,a)";
    static char fmt_502[] = "(\002 INCLUDING\002,i5,\002 RECORDS FOR THE COV\
ARIANCE MATRIX.\002/)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_rsfe(cilist *), e_rsfe(), f_open(olist *), f_rew(alist *);

    /* Local variables */
    static integer i__, iint, npar2;
    static logical lname, lopen;
    static char cgname[64], cfname[64], canswr[1];
    static integer nlines, ncovar;

    /* Fortran I/O blocks */
    static cilist io___1004 = { 0, 0, 0, fmt_32, 0 };
    static cilist io___1005 = { 0, 0, 0, fmt_35, 0 };
    static cilist io___1006 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1007 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1009 = { 0, 0, 0, fmt_37, 0 };
    static cilist io___1010 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1012 = { 1, 0, 0, "(10HSET TITLE )", 0 };
    static cilist io___1013 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1014 = { 0, 0, 0, "(10HPARAMETERS)", 0 };
    static cilist io___1018 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___1019 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___1020 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1021 = { 1, 0, 0, fmt_1003, 0 };
    static cilist io___1023 = { 0, 0, 0, fmt_1004, 0 };
    static cilist io___1025 = { 0, 0, 0, fmt_501, 0 };
    static cilist io___1026 = { 0, 0, 0, fmt_502, 0 };
    static cilist io___1027 = { 0, 0, 0, "(A,I4)", 0 };
    static cilist io___1028 = { 0, 0, 0, "(A,I4,A)", 0 };
    static cilist io___1029 = { 0, 0, 0, "(A,I4)", 0 };
    static cilist io___1030 = { 0, 0, 0, "(A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       Writes current parameter values and step sizes onto file ISYSSA */
/* C          in format which can be reread by Minuit for restarting. */
/* C       The covariance matrix is also output if it exists. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



    ioin__1.inerr = 0;
    ioin__1.inunit = mn7iou_1.isyssa;
    ioin__1.infile = 0;
    ioin__1.inex = 0;
    ioin__1.inopen = &lopen;
    ioin__1.innum = 0;
    ioin__1.innamed = &lname;
    ioin__1.innamlen = 64;
    ioin__1.inname = cgname;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (lopen) {
	if (! lname) {
	    s_copy(cgname, "UNNAMED FILE", (ftnlen)64, (ftnlen)12);
	}
	io___1004.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1004);
	do_fio(&c__1, (char *)&mn7iou_1.isyssa, (ftnlen)sizeof(integer));
	do_fio(&c__1, cgname, (ftnlen)64);
	e_wsfe();
    } else {
/*                new file, open it */
	io___1005.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1005);
	do_fio(&c__1, (char *)&mn7iou_1.isyssa, (ftnlen)sizeof(integer));
	e_wsfe();
	if (mn7flg_1.isw[5] == 1) {
	    io___1006.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1006);
	    do_fio(&c__1, " PLEASE GIVE FILE NAME:", (ftnlen)23);
	    e_wsfe();
	    io___1007.ciunit = mn7iou_1.isysrd;
	    s_rsfe(&io___1007);
	    do_fio(&c__1, cfname, (ftnlen)64);
	    e_rsfe();
	    o__1.oerr = 1;
	    o__1.ounit = mn7iou_1.isyssa;
	    o__1.ofnmlen = 64;
	    o__1.ofnm = cfname;
	    o__1.orl = 0;
	    o__1.osta = "NEW";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L600;
	    }
	    s_copy(cgname, cfname, (ftnlen)64, (ftnlen)64);
	} else {
	    goto L650;
	}
    }
/*                               file is now correctly opened */
    if (mn7flg_1.isw[5] == 1) {
	io___1009.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1009);
	do_fio(&c__1, (char *)&mn7iou_1.isyssa, (ftnlen)sizeof(integer));
	e_wsfe();
	io___1010.ciunit = mn7iou_1.isysrd;
	s_rsfe(&io___1010);
	do_fio(&c__1, canswr, (ftnlen)1);
	e_rsfe();
	if (*(unsigned char *)canswr == 'Y' || *(unsigned char *)canswr == 
		'y') {
	    al__1.aerr = 0;
	    al__1.aunit = mn7iou_1.isyssa;
	    f_rew(&al__1);
	}
    }
/*                               and rewound if requested */
    io___1012.ciunit = mn7iou_1.isyssa;
    i__1 = s_wsfe(&io___1012);
    if (i__1 != 0) {
	goto L700;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L700;
    }
    io___1013.ciunit = mn7iou_1.isyssa;
    s_wsfe(&io___1013);
    do_fio(&c__1, mn7tit_1.ctitl, (ftnlen)50);
    e_wsfe();
    io___1014.ciunit = mn7iou_1.isyssa;
    s_wsfe(&io___1014);
    e_wsfe();
    nlines = 3;
/*                                write out parameter values */
    i__1 = mn7npr_1.nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (mn7inx_1.nvarl[i__ - 1] < 0) {
	    goto L200;
	}
	++nlines;
	iint = mn7inx_1.niofex[i__ - 1];
	if (mn7inx_1.nvarl[i__ - 1] > 1) {
	    goto L100;
	}
/*         parameter without limits */
	io___1018.ciunit = mn7iou_1.isyssa;
	s_wsfe(&io___1018);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, mn7nam_1.cpnam + (i__ - 1) * 10, (ftnlen)10);
	do_fio(&c__1, (char *)&mn7ext_1.u[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, (char *)&mn7err_1.werr[iint - 1], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
	goto L200;
/*         parameter with limits */
L100:
	io___1019.ciunit = mn7iou_1.isyssa;
	s_wsfe(&io___1019);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, mn7nam_1.cpnam + (i__ - 1) * 10, (ftnlen)10);
	do_fio(&c__1, (char *)&mn7ext_1.u[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, (char *)&mn7err_1.werr[iint - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&mn7ext_1.alim[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&mn7ext_1.blim[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
L200:
	;
    }
    io___1020.ciunit = mn7iou_1.isyssa;
    s_wsfe(&io___1020);
    do_fio(&c__1, " ", (ftnlen)1);
    e_wsfe();
    ++nlines;
/*                                  write out covariance matrix, if any */
    if (mn7flg_1.isw[1] < 1) {
	goto L750;
    }
    io___1021.ciunit = mn7iou_1.isyssa;
    i__1 = s_wsfe(&io___1021);
    if (i__1 != 0) {
	goto L700;
    }
    i__1 = do_fio(&c__1, (char *)&mn7npr_1.npar, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L700;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L700;
    }
    npar2 = mn7npr_1.npar * (mn7npr_1.npar + 1) / 2;
    io___1023.ciunit = mn7iou_1.isyssa;
    s_wsfe(&io___1023);
    i__1 = npar2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&mn7var_1.vhmat[i__ - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_wsfe();
    ncovar = npar2 / 7 + 1;
    if (npar2 % 7 > 0) {
	++ncovar;
    }
    nlines += ncovar;
    io___1025.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1025);
    do_fio(&c__1, (char *)&nlines, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&mn7iou_1.isyssa, (ftnlen)sizeof(integer));
    do_fio(&c__1, cgname, (ftnlen)45);
    e_wsfe();
    if (ncovar > 0) {
	io___1026.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1026);
	do_fio(&c__1, (char *)&ncovar, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L900;
/*                                           some error conditions */
L600:
    io___1027.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1027);
    do_fio(&c__1, " I/O ERROR: UNABLE TO OPEN UNIT", (ftnlen)31);
    do_fio(&c__1, (char *)&mn7iou_1.isyssa, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L900;
L650:
    io___1028.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1028);
    do_fio(&c__1, " UNIT", (ftnlen)5);
    do_fio(&c__1, (char *)&mn7iou_1.isyssa, (ftnlen)sizeof(integer));
    do_fio(&c__1, " IS NOT OPENED.", (ftnlen)15);
    e_wsfe();
    goto L900;
L700:
    io___1029.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1029);
    do_fio(&c__1, " ERROR: UNABLE TO WRITE TO UNIT", (ftnlen)31);
    do_fio(&c__1, (char *)&mn7iou_1.isyssa, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L900;
L750:
    io___1030.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1030);
    do_fio(&c__1, " THERE IS NO COVARIANCE MATRIX TO SAVE.", (ftnlen)39);
    e_wsfe();

L900:
    return 0;
} /* mnsave_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:51  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnscan_(S_fp fcn, U_fp futil)
{
    /* Format strings */
    static char fmt_1001[] = "(i1,\002SCAN OF PARAMETER NO.\002,i3,\002, \
 \002,a10)";
    static char fmt_1000[] = "(\002 REQUESTED RANGE OUTSIDE LIMITS FOR PARAM\
ETER \002,i3/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static doublereal xh, xl;
    static integer ipar, iint;
    static doublereal step;
    static integer icall, ncall;
    static doublereal uhigh;
    static integer nbins;
    static doublereal xhreq, ubest, xlreq;
    static integer nparx;
    static doublereal fnext;
    static integer nunit;
    static doublereal unext;
    static integer nxypt, nccall;
    extern /* Subroutine */ int mnamin_(S_fp, U_fp);
    static integer iparwd;
    extern /* Subroutine */ int mnbins_(doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, doublereal *), mnexin_(
	    doublereal *), mnplot_(doublereal *, doublereal *, char *, 
	    integer *, integer *, integer *, integer *, ftnlen), mnprin_(
	    integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___1049 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___1051 = { 0, 0, 0, fmt_1000, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Scans the values of FCN as a function of one parameter */
/* C        and plots the resulting values as a curve using MNPLOT. */
/* C        It may be called to scan one parameter or all parameters. */
/* C        retains the best function and parameter values found. */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    xlreq = min(mn7arg_1.word7[2],mn7arg_1.word7[3]);
    xhreq = max(mn7arg_1.word7[2],mn7arg_1.word7[3]);
    ncall = (integer) (mn7arg_1.word7[1] + (float).01);
    if (ncall <= 1) {
	ncall = 41;
    }
    if (ncall > 101) {
	ncall = 101;
    }
    nccall = ncall;
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    iparwd = (integer) (mn7arg_1.word7[0] + (float).1);
    ipar = max(iparwd,0);
    iint = mn7inx_1.niofex[ipar - 1];
    s_copy(mn7tit_1.cstatu, "NO CHANGE", (ftnlen)10, (ftnlen)9);
    if (iparwd > 0) {
	goto L200;
    }

/*         equivalent to a loop over parameters requested */
L100:
    ++ipar;
    if (ipar > mn7npr_1.nu) {
	goto L900;
    }
    iint = mn7inx_1.niofex[ipar - 1];
    if (iint <= 0) {
	goto L100;
    }
/*         set up range for parameter IPAR */
L200:
    ubest = mn7ext_1.u[ipar - 1];
    mn7rpt_1.xpt[0] = ubest;
    mn7rpt_1.ypt[0] = mn7min_1.amin;
    *(unsigned char *)&mn7cpt_1.chpt[0] = ' ';
    mn7rpt_1.xpt[1] = ubest;
    mn7rpt_1.ypt[1] = mn7min_1.amin;
    *(unsigned char *)&mn7cpt_1.chpt[1] = 'X';
    nxypt = 2;
    if (mn7inx_1.nvarl[ipar - 1] > 1) {
	goto L300;
    }
/*         no limits on parameter */
    if (xlreq == xhreq) {
	goto L250;
    }
    unext = xlreq;
    step = (xhreq - xlreq) / (real) (ncall - 1);
    goto L500;
L250:
    xl = ubest - mn7err_1.werr[iint - 1];
    xh = ubest + mn7err_1.werr[iint - 1];
    mnbins_(&xl, &xh, &ncall, &unext, &uhigh, &nbins, &step);
    nccall = nbins + 1;
    goto L500;
/*         limits on parameter */
L300:
    if (xlreq == xhreq) {
	goto L350;
    }
/* Computing MAX */
    d__1 = xlreq, d__2 = mn7ext_1.alim[ipar - 1];
    xl = max(d__1,d__2);
/* Computing MIN */
    d__1 = xhreq, d__2 = mn7ext_1.blim[ipar - 1];
    xh = min(d__1,d__2);
    if (xl >= xh) {
	goto L700;
    }
    unext = xl;
    step = (xh - xl) / (real) (ncall - 1);
    goto L500;
L350:
    unext = mn7ext_1.alim[ipar - 1];
    step = (mn7ext_1.blim[ipar - 1] - mn7ext_1.alim[ipar - 1]) / (real) (
	    ncall - 1);
/*         main scanning loop over parameter IPAR */
L500:
    i__1 = nccall;
    for (icall = 1; icall <= i__1; ++icall) {
	mn7ext_1.u[ipar - 1] = unext;
	nparx = mn7npr_1.npar;
	(*fcn)(&nparx, mn7der_1.gin, &fnext, mn7ext_1.u, &c__4, (U_fp)futil);
	++mn7cnv_1.nfcn;
	++nxypt;
	mn7rpt_1.xpt[nxypt - 1] = unext;
	mn7rpt_1.ypt[nxypt - 1] = fnext;
	*(unsigned char *)&mn7cpt_1.chpt[nxypt - 1] = '*';
	if (fnext < mn7min_1.amin) {
	    mn7min_1.amin = fnext;
	    ubest = unext;
	    s_copy(mn7tit_1.cstatu, "IMPROVED  ", (ftnlen)10, (ftnlen)10);
	}
/* L530: */
	unext += step;
/* L600: */
    }
/*         finished with scan of parameter IPAR */
    mn7ext_1.u[ipar - 1] = ubest;
    mnexin_(mn7int_1.x);
    if (mn7flg_1.isw[4] >= 1) {
	io___1049.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1049);
	do_fio(&c__1, (char *)&mn7iou_1.newpag, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipar, (ftnlen)sizeof(integer));
	do_fio(&c__1, mn7nam_1.cpnam + (ipar - 1) * 10, (ftnlen)10);
	e_wsfe();
	nunit = mn7iou_1.isyswr;
	mnplot_(mn7rpt_1.xpt, mn7rpt_1.ypt, mn7cpt_1.chpt, &nxypt, &nunit, &
		mn7iou_1.npagwd, &mn7iou_1.npagln, (ftnlen)1);
    }
    goto L800;
L700:
    io___1051.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1051);
    do_fio(&c__1, (char *)&ipar, (ftnlen)sizeof(integer));
    e_wsfe();
L800:
    if (iparwd <= 0) {
	goto L100;
    }
/*         finished with all parameters */
L900:
    if (mn7flg_1.isw[4] >= 0) {
	mnprin_(&c__5, &mn7min_1.amin);
    }
    return 0;
} /* mnscan_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnseek_(S_fp fcn, U_fp futil)
{
    /* Format strings */
    static char fmt_3[] = "(\002 MNSEEK: MONTE CARLO MINIMIZATION USING METR\
OPOLIS\002,\002 ALGORITHM\002/\002 TO STOP AFTER\002,i6,\002 SUCCESSIVE FAIL\
URES, OR\002,i7,\002 STEPS\002/\002 MAXIMUM STEP SIZE IS\002,f9.3,\002 ERROR\
 BARS.\002)";
    static char fmt_601[] = "(\002 MNSEEK:\002,i5,\002 SUCCESSIVE UNSUCCESSF\
UL TRIALS.\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double log(doublereal);

    /* Local variables */
    static integer j, ib;
    static doublereal bar, dxdi;
    static integer ipar;
    static doublereal xmid[100];
    static integer iext;
    static doublereal rnum, ftry, rnum1, rnum2;
    static integer ifail;
    static doublereal alpha;
    static integer iseed;
    static doublereal flast, xbest[100];
    static integer nparx, istep;
    extern /* Subroutine */ int mnrn15_(doublereal *, integer *);
    static integer mxfail;
    extern /* Subroutine */ int mnamin_(S_fp, U_fp);
    static integer mxstep;
    extern /* Subroutine */ int mnprin_(integer *, doublereal *), mndxdi_(
	    doublereal *, integer *, doublereal *), mninex_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___1055 = { 0, 0, 0, fmt_3, 0 };
    static cilist io___1073 = { 0, 0, 0, fmt_601, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C   Performs a rough (but global) minimization by monte carlo search. */
/* C        Each time a new minimum is found, the search area is shifted */
/* C        to be centered at the best value.  Random points are chosen */
/* C        uniformly over a hypercube determined by current step sizes. */
/* C   The Metropolis algorithm accepts a worse point with probability */
/* C      exp(-d/UP), where d is the degradation.  Improved points */
/* C      are of course always accepted.  Actual steps are random */
/* C      multiples of the nominal steps (DIRIN). */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    mxfail = (integer) mn7arg_1.word7[0];
    if (mxfail <= 0) {
	mxfail = mn7npr_1.npar * 20 + 100;
    }
    mxstep = mxfail * 10;
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    alpha = mn7arg_1.word7[1];
    if (alpha <= 0.) {
	alpha = (float)3.;
    }
    if (mn7flg_1.isw[4] >= 1) {
	io___1055.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1055);
	do_fio(&c__1, (char *)&mxfail, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&mxstep, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&alpha, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    s_copy(mn7tit_1.cstatu, "INITIAL  ", (ftnlen)10, (ftnlen)9);
    if (mn7flg_1.isw[4] >= 2) {
	mnprin_(&c__2, &mn7min_1.amin);
    }
    s_copy(mn7tit_1.cstatu, "UNCHANGED ", (ftnlen)10, (ftnlen)10);
    ifail = 0;
    rnum = 0.;
    rnum1 = 0.;
    rnum2 = 0.;
    nparx = mn7npr_1.npar;
    flast = mn7min_1.amin;
/*              set up step sizes, starting values */
    i__1 = mn7npr_1.npar;
    for (ipar = 1; ipar <= i__1; ++ipar) {
	iext = mn7inx_1.nexofi[ipar - 1];
	mn7int_1.dirin[ipar - 1] = alpha * (float)2. * mn7err_1.werr[ipar - 1]
		;
	if (mn7inx_1.nvarl[iext - 1] > 1) {
/*              parameter with limits */
	    mndxdi_(&mn7int_1.x[ipar - 1], &ipar, &dxdi);
	    if (dxdi == 0.) {
		dxdi = (float)1.;
	    }
	    mn7int_1.dirin[ipar - 1] = alpha * (float)2. * mn7err_1.werr[ipar 
		    - 1] / dxdi;
	    if ((d__1 = mn7int_1.dirin[ipar - 1], abs(d__1)) > 
		    6.2831859999999997) {
		mn7int_1.dirin[ipar - 1] = 6.2831859999999997;
	    }
	}
	xmid[ipar - 1] = mn7int_1.x[ipar - 1];
/* L10: */
	xbest[ipar - 1] = mn7int_1.x[ipar - 1];
    }
/*                              search loop */
    i__1 = mxstep;
    for (istep = 1; istep <= i__1; ++istep) {
	if (ifail >= mxfail) {
	    goto L600;
	}
	i__2 = mn7npr_1.npar;
	for (ipar = 1; ipar <= i__2; ++ipar) {
	    mnrn15_(&rnum1, &iseed);
	    mnrn15_(&rnum2, &iseed);
/* L100: */
	    mn7int_1.x[ipar - 1] = xmid[ipar - 1] + (rnum1 + rnum2 - (float)
		    1.) * (float).5 * mn7int_1.dirin[ipar - 1];
	}
	mninex_(mn7int_1.x);
	(*fcn)(&nparx, mn7der_1.gin, &ftry, mn7ext_1.u, &c__4, (U_fp)futil);
	++mn7cnv_1.nfcn;
	if (ftry < flast) {
	    if (ftry < mn7min_1.amin) {
		s_copy(mn7tit_1.cstatu, "IMPROVEMNT", (ftnlen)10, (ftnlen)10);
		mn7min_1.amin = ftry;
		i__2 = mn7npr_1.npar;
		for (ib = 1; ib <= i__2; ++ib) {
/* L200: */
		    xbest[ib - 1] = mn7int_1.x[ib - 1];
		}
		ifail = 0;
		if (mn7flg_1.isw[4] >= 2) {
		    mnprin_(&c__2, &mn7min_1.amin);
		}
	    }
	    goto L300;
	} else {
	    ++ifail;
/*                   Metropolis algorithm */
	    bar = (mn7min_1.amin - ftry) / mn7min_1.up;
	    mnrn15_(&rnum, &iseed);
	    if (bar < log(rnum)) {
		goto L500;
	    }
	}
/*                    Accept new point, move there */
L300:
	i__2 = mn7npr_1.npar;
	for (j = 1; j <= i__2; ++j) {
	    xmid[j - 1] = mn7int_1.x[j - 1];
/* L350: */
	}
	flast = ftry;
L500:
	;
    }
/*                               end search loop */
L600:
    if (mn7flg_1.isw[4] > 1) {
	io___1073.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1073);
	do_fio(&c__1, (char *)&ifail, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    i__1 = mn7npr_1.npar;
    for (ib = 1; ib <= i__1; ++ib) {
/* L700: */
	mn7int_1.x[ib - 1] = xbest[ib - 1];
    }
    mninex_(mn7int_1.x);
    if (mn7flg_1.isw[4] >= 1) {
	mnprin_(&c__2, &mn7min_1.amin);
    }
    if (mn7flg_1.isw[4] == 0) {
	mnprin_(&c__0, &mn7min_1.amin);
    }
    return 0;
} /* mnseek_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:52  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni */
/* Minuit */


/* Subroutine */ int mnset_(S_fp fcn, U_fp futil)
{
    /* Initialized data */

    static char cname[10*30+1] = "FCN value PARametersLIMits    COVarianceCO\
RrelatioPRInt levlNOGradientGRAdient  ERRor def INPut fileWIDth pageLINes pa\
geNOWarningsWARnings  RANdom genTITle     STRategy  EIGenvaluePAGe throwMINo\
s errsEPSmachineOUTputfileBATch     INTeractivVERsion   reserve   NODebug   \
DEBug     SHOw      SET       ";
    static integer nname = 25;
    static integer nntot = 30;
    static char cprlev[34*5+1] = "-1: NO OUTPUT EXCEPT FROM \"SHOW\"   0: RE\
DUCED OUTPUT                 1: NORMAL OUTPUT                  2: EXTRA OUTP\
UT FOR PROBLEM CASES 3: MAXIMUM OUTPUT                ";
    static char cstrat[44*3+1] = " 0: MINIMIZE THE NUMBER OF CALLS TO FUNCTI\
ON 1: TRY TO BALANCE SPEED AGAINST RELIABILITY 2: MAKE SURE MINIMUM TRUE, ER\
RORS CORRECT  ";
    static char cdbopt[40*7+1] = "REPORT ALL EXCEPTIONAL CONDITIONS       MN\
LINE: LINE SEARCH MINIMIZATION        MNDERI: FIRST DERIVATIVE CALCULATIONS \
  MNHESS: SECOND DERIVATIVE CALCULATIONS  MNMIGR: COVARIANCE MATRIX UPDATES \
      MNHES1: FIRST DERIVATIVE UNCERTAINTIES  MNCONT: MNCONTOUR PLOT (MNCROS\
 SEARCH)  ";

    /* Format strings */
    static char fmt_151[] = "(\002 MINUIT RANDOM NUMBER SEED SET TO \002,i10)"
	    ;
    static char fmt_289[] = "(\002 UNKNOWN DEBUG OPTION\002,i6,\002 REQUESTE\
D. IGNORED\002)";
    static char fmt_1061[] = "(/\002 CURRENT PRINTOUT LEVEL IS \002,a)";
    static char fmt_1081[] = "(\002 NOGRAD IS SET.  DERIVATIVES NOT COMPUTED\
 IN FCN.\002)";
    static char fmt_1082[] = "(\002   GRAD IS SET.  USER COMPUTES DERIVATIVE\
S IN FCN.\002)";
    static char fmt_1091[] = "(\002 ERRORS CORRESPOND TO FUNCTION CHANGE O\
F\002,g13.5)";
    static char fmt_1002[] = "(\002 INPUT NOW BEING READ IN \002,a,\002 FROM\
 UNIT NO.\002,i3/\002 FILENAME: \002,a)";
    static char fmt_1111[] = "(10x,\002PAGE WIDTH IS SET TO\002,i4,\002 COLU\
MNS\002)";
    static char fmt_1121[] = "(10x,\002PAGE LENGTH IS SET TO\002,i4,\002 LIN\
ES\002)";
    static char fmt_1141[] = "(\002 MINUIT WARNING MESSAGES ARE \002,a)";
    static char fmt_1151[] = "(\002 MINUIT RNDM SEED IS CURRENTLY=\002,i10/)";
    static char fmt_1175[] = "(/\002 NOW USING STRATEGY \002,a/)";
    static char fmt_1286[] = "(10x,\002DEBUG OPTION\002,i3,\002 IS \002,a3\
,\002 :\002,a)";
    static char fmt_1901[] = "(\002 THE COMMAND:\002,a10,\002 IS UNKNOWN.\
\002/)";
    static char fmt_2101[] = "(\002 THE FORMAT OF THE \002,a4,\002 COMMAND I\
S:\002//1x,a4,\002 xxx    [numerical arguments if any]\002//\002 WHERE xxx M\
AY BE ONE OF THE FOLLOWING:\002/(7x,6a12))";

    /* System generated locals */
    integer i__1;
    inlist ioin__1;

    /* Builtin functions */
    integer i_indx(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double sqrt(doublereal);
    integer f_inqu(inlist *);

    /* Local variables */
    static integer i__, id, ii, kk;
    static doublereal val;
    static integer isw2, iset;
    static char copt[3];
    static integer iprm;
    static char cmode[16], ckind[4];
    static integer iseed, jseed, kname;
    static logical lname;
    static char cwarn[10];
    extern /* Subroutine */ int mnrn15_(doublereal *, integer *);
    static integer iunit;
    static char cfname[64];
    static integer ikseed;
    extern /* Subroutine */ int mngrad_(S_fp, U_fp);
    static integer idbopt;
    extern /* Subroutine */ int mnexin_(doublereal *), mnrset_(integer *), 
	    mnlims_(), mnwerr_(), mnwarn_(char *, char *, char *, ftnlen, 
	    ftnlen, ftnlen), mnamin_(S_fp, U_fp), mnprin_(integer *, 
	    doublereal *), mnmatu_(integer *);
    static integer igrain, iswsav;
    extern /* Subroutine */ int mnpsdf_();

    /* Fortran I/O blocks */
    static cilist io___1085 = { 0, 0, 0, "(A/)", 0 };
    static cilist io___1088 = { 0, 0, 0, fmt_151, 0 };
    static cilist io___1093 = { 0, 0, 0, fmt_289, 0 };
    static cilist io___1094 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1095 = { 0, 0, 0, "(27X,A)", 0 };
    static cilist io___1096 = { 0, 0, 0, fmt_1061, 0 };
    static cilist io___1097 = { 0, 0, 0, fmt_1081, 0 };
    static cilist io___1098 = { 0, 0, 0, fmt_1082, 0 };
    static cilist io___1099 = { 0, 0, 0, fmt_1091, 0 };
    static cilist io___1103 = { 0, 0, 0, fmt_1002, 0 };
    static cilist io___1104 = { 0, 0, 0, fmt_1111, 0 };
    static cilist io___1105 = { 0, 0, 0, fmt_1121, 0 };
    static cilist io___1107 = { 0, 0, 0, fmt_1141, 0 };
    static cilist io___1110 = { 0, 0, 0, fmt_1151, 0 };
    static cilist io___1112 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___1113 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1114 = { 0, 0, 0, "(20X,A)", 0 };
    static cilist io___1115 = { 0, 0, 0, fmt_1175, 0 };
    static cilist io___1117 = { 0, 0, 0, "(1X,A)", 0 };
    static cilist io___1118 = { 0, 0, 0, "(A,I3)", 0 };
    static cilist io___1119 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1121 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1122 = { 0, 0, 0, "(A,E12.3)", 0 };
    static cilist io___1123 = { 0, 0, 0, "(A,I4)", 0 };
    static cilist io___1124 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___1126 = { 0, 0, 0, fmt_1286, 0 };
    static cilist io___1127 = { 0, 0, 0, fmt_1901, 0 };
    static cilist io___1128 = { 0, 0, 0, fmt_2101, 0 };
    static cilist io___1130 = { 0, 0, 0, "(' ABOVE COMMAND IS ILLEGAL.   IGN\
ORED')", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Called from MNEXCM */
/* C        Interprets the commands that start with SET and SHOW */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



/*        file characteristics for SET INPUT */
/*       'SET ' or 'SHOW',  'ON ' or 'OFF', 'SUPPRESSED' or 'REPORTED  ' */
/*        explanation of print level numbers -1:3  and strategies 0:2 */
/*        identification of debug options */
/*        things that can be set or shown */
/*        options not intended for normal users */





    i__1 = nntot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i_indx(mn7tit_1.cword + 3, cname + (i__ - 1) * 10, (ftnlen)7, (
		ftnlen)3) > 0) {
	    goto L5;
	}
/* L2: */
    }
    i__ = 0;
L5:
    kname = i__;

/*           Command could be SET xxx, SHOW xxx,  HELP SET or HELP SHOW */
    if (i_indx(mn7tit_1.cword, "HEL", (ftnlen)4, (ftnlen)3) > 0) {
	goto L2000;
    }
    if (i_indx(mn7tit_1.cword, "SHO", (ftnlen)4, (ftnlen)3) > 0) {
	goto L1000;
    }
    if (i_indx(mn7tit_1.cword, "SET", (ftnlen)4, (ftnlen)3) == 0) {
	goto L1900;
    }
/*                           --- */
    s_copy(ckind, "SET ", (ftnlen)4, (ftnlen)4);
/*                                        . . . . . . . . . . set unknown */
    if (kname <= 0) {
	goto L1900;
    }
/*                                        . . . . . . . . . . set known */
    switch (kname) {
	case 1:  goto L3000;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	case 5:  goto L3000;
	case 6:  goto L60;
	case 7:  goto L70;
	case 8:  goto L80;
	case 9:  goto L90;
	case 10:  goto L100;
	case 11:  goto L110;
	case 12:  goto L120;
	case 13:  goto L130;
	case 14:  goto L140;
	case 15:  goto L150;
	case 16:  goto L160;
	case 17:  goto L170;
	case 18:  goto L3000;
	case 19:  goto L190;
	case 20:  goto L3000;
	case 21:  goto L210;
	case 22:  goto L220;
	case 23:  goto L230;
	case 24:  goto L240;
	case 25:  goto L3000;
	case 26:  goto L1900;
	case 27:  goto L270;
	case 28:  goto L280;
	case 29:  goto L290;
	case 30:  goto L300;
    }

/*                                        . . . . . . . . . . set param */
L20:
    iprm = (integer) mn7arg_1.word7[0];
    if (iprm > mn7npr_1.nu) {
	goto L25;
    }
    if (iprm <= 0) {
	goto L25;
    }
    if (mn7inx_1.nvarl[iprm - 1] < 0) {
	goto L25;
    }
    mn7ext_1.u[iprm - 1] = mn7arg_1.word7[1];
    mnexin_(mn7int_1.x);
    isw2 = mn7flg_1.isw[1];
    mnrset_(&c__1);
/*        Keep approximate covariance matrix, even if new param value */
    mn7flg_1.isw[1] = min(isw2,1);
    s_copy(mn7tit_1.cfrom, "SET PARM", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "NEW VALUES", (ftnlen)10, (ftnlen)10);
    goto L4000;
L25:
    io___1085.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1085);
    do_fio(&c__1, " UNDEFINED PARAMETER NUMBER.  IGNORED.", (ftnlen)38);
    e_wsfe();
    goto L4000;
/*                                        . . . . . . . . . . set limits */
L30:
    mnlims_();
    goto L4000;
/*                                        . . . . . . . . . . set covar */
L40:
/*   this command must be handled by MNREAD, and is not Fortran-callable */
    goto L3000;
/*                                        . . . . . . . . . . set print */
L60:
    mn7flg_1.isw[4] = (integer) mn7arg_1.word7[0];
    goto L4000;
/*                                        . . . . . . . . . . set nograd */
L70:
    mn7flg_1.isw[2] = 0;
    goto L4000;
/*                                        . . . . . . . . . . set grad */
L80:
    mngrad_((S_fp)fcn, (U_fp)futil);
    goto L4000;
/*                                        . . . . . . . . . . set errdef */
L90:
    if (mn7arg_1.word7[0] == mn7min_1.up) {
	goto L4000;
    }
    if (mn7arg_1.word7[0] <= 0.) {
	if (mn7min_1.up == mn7cns_1.updflt) {
	    goto L4000;
	}
	mn7min_1.up = mn7cns_1.updflt;
    } else {
	mn7min_1.up = mn7arg_1.word7[0];
    }
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7err_1.ern[i__ - 1] = (float)0.;
/* L95: */
	mn7err_1.erp[i__ - 1] = (float)0.;
    }
    mnwerr_();
    goto L4000;
/*                                        . . . . . . . . . . set input */
/* This command must be handled by MNREAD. If it gets this far, */
/*         it is illegal. */
L100:
    goto L3000;
/*                                        . . . . . . . . . . set width */
L110:
    mn7iou_1.npagwd = (integer) mn7arg_1.word7[0];
    mn7iou_1.npagwd = max(mn7iou_1.npagwd,50);
    goto L4000;
/*                                        . . . . . . . . . . set lines */
L120:
    mn7iou_1.npagln = (integer) mn7arg_1.word7[0];
    goto L4000;
/*                                        . . . . . . . . . . set nowarn */
L130:
    mn7log_1.lwarn = FALSE_;
    goto L4000;
/*                                        . . . . . . . . . . set warn */
L140:
    mn7log_1.lwarn = TRUE_;
    mnwarn_("W", "SHO", "SHO", (ftnlen)1, (ftnlen)3, (ftnlen)3);
    goto L4000;
/*                                        . . . . . . . . . . set random */
L150:
    jseed = (integer) mn7arg_1.word7[0];
    val = (float)3.;
    mnrn15_(&val, &jseed);
    if (mn7flg_1.isw[4] > 0) {
	io___1088.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1088);
	do_fio(&c__1, (char *)&jseed, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L4000;
/*                                        . . . . . . . . . . set title */
L160:
/*   this command must be handled by MNREAD, and is not Fortran-callable */
    goto L3000;
/*                                        . . . . . . . . . set strategy */
L170:
    mn7cnv_1.istrat = (integer) mn7arg_1.word7[0];
    mn7cnv_1.istrat = max(mn7cnv_1.istrat,0);
    mn7cnv_1.istrat = min(mn7cnv_1.istrat,2);
    if (mn7flg_1.isw[4] > 0) {
	goto L1172;
    }
    goto L4000;
/*                                       . . . . . . . . . set page throw */
L190:
    mn7iou_1.newpag = (integer) mn7arg_1.word7[0];
    goto L1190;
/*                                        . . . . . . . . . . set epsmac */
L210:
    if (mn7arg_1.word7[0] > 0. && mn7arg_1.word7[0] < (float).1) {
	mn7cns_1.epsmac = mn7arg_1.word7[0];
    }
    mn7cns_1.epsma2 = sqrt(mn7cns_1.epsmac);
    goto L1210;
/*                                        . . . . . . . . . . set outputfile */
L220:
    iunit = (integer) mn7arg_1.word7[0];
    mn7iou_1.isyswr = iunit;
    mn7io2_1.istkwr[0] = iunit;
    if (mn7flg_1.isw[4] >= 0) {
	goto L1220;
    }
    goto L4000;
/*                                        . . . . . . . . . . set batch */
L230:
    mn7flg_1.isw[5] = 0;
    if (mn7flg_1.isw[4] >= 0) {
	goto L1100;
    }
    goto L4000;
/*                                        . . . . . . . . . . set interactive */
L240:
    mn7flg_1.isw[5] = 1;
    if (mn7flg_1.isw[4] >= 0) {
	goto L1100;
    }
    goto L4000;
/*                                        . . . . . . . . . . set nodebug */
L270:
    iset = 0;
    goto L281;
/*                                        . . . . . . . . . . set debug */
L280:
    iset = 1;
L281:
    idbopt = (integer) mn7arg_1.word7[0];
    if (idbopt > 6) {
	goto L288;
    }
    if (idbopt >= 0) {
	mn7flg_1.idbg[idbopt] = iset;
	if (iset == 1) {
	    mn7flg_1.idbg[0] = 1;
	}
    } else {
/*             SET DEBUG -1  sets all debug options */
	for (id = 0; id <= 6; ++id) {
/* L285: */
	    mn7flg_1.idbg[id] = iset;
	}
    }
    mn7log_1.lrepor = mn7flg_1.idbg[0] >= 1;
    mnwarn_("D", "SHO", "SHO", (ftnlen)1, (ftnlen)3, (ftnlen)3);
    goto L4000;
L288:
    io___1093.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1093);
    do_fio(&c__1, (char *)&idbopt, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L4000;
/*                                        . . . . . . . . . . set show */
L290:
/*                                        . . . . . . . . . . set set */
L300:
    goto L3000;
/*                ----------------------------------------------------- */
L1000:
/*               at this point, CWORD must be 'SHOW' */
    s_copy(ckind, "SHOW", (ftnlen)4, (ftnlen)4);
    if (kname <= 0) {
	goto L1900;
    }
    switch (kname) {
	case 1:  goto L1010;
	case 2:  goto L1020;
	case 3:  goto L1030;
	case 4:  goto L1040;
	case 5:  goto L1050;
	case 6:  goto L1060;
	case 7:  goto L1070;
	case 8:  goto L1070;
	case 9:  goto L1090;
	case 10:  goto L1100;
	case 11:  goto L1110;
	case 12:  goto L1120;
	case 13:  goto L1130;
	case 14:  goto L1130;
	case 15:  goto L1150;
	case 16:  goto L1160;
	case 17:  goto L1170;
	case 18:  goto L1180;
	case 19:  goto L1190;
	case 20:  goto L1200;
	case 21:  goto L1210;
	case 22:  goto L1220;
	case 23:  goto L1100;
	case 24:  goto L1100;
	case 25:  goto L1250;
	case 26:  goto L1900;
	case 27:  goto L1270;
	case 28:  goto L1270;
	case 29:  goto L1290;
	case 30:  goto L1300;
    }

/*                                        . . . . . . . . . . show fcn */
L1010:
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    mnprin_(&c__0, &mn7min_1.amin);
    goto L4000;
/*                                        . . . . . . . . . . show param */
L1020:
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    mnprin_(&c__5, &mn7min_1.amin);
    goto L4000;
/*                                        . . . . . . . . . . show limits */
L1030:
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    mnprin_(&c__1, &mn7min_1.amin);
    goto L4000;
/*                                        . . . . . . . . . . show covar */
L1040:
    mnmatu_(&c__1);
    goto L4000;
/*                                        . . . . . . . . . . show corre */
L1050:
    mnmatu_(&c__0);
    goto L4000;
/*                                        . . . . . . . . . . show print */
L1060:
    if (mn7flg_1.isw[4] < -1) {
	mn7flg_1.isw[4] = -1;
    }
    if (mn7flg_1.isw[4] > 3) {
	mn7flg_1.isw[4] = 3;
    }
    io___1094.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1094);
    do_fio(&c__1, " ALLOWED PRINT LEVELS ARE:", (ftnlen)26);
    e_wsfe();
    io___1095.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1095);
    do_fio(&c__5, cprlev, (ftnlen)34);
    e_wsfe();
    io___1096.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1096);
    do_fio(&c__1, cprlev + (mn7flg_1.isw[4] + 1) * 34, (ftnlen)34);
    e_wsfe();
    goto L4000;
/*                                        . . . . . . . show nograd, grad */
L1070:
    if (mn7flg_1.isw[2] <= 0) {
	io___1097.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1097);
	e_wsfe();
    } else {
	io___1098.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1098);
	e_wsfe();
    }
    goto L4000;
/*                                       . . . . . . . . . . show errdef */
L1090:
    io___1099.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1099);
    do_fio(&c__1, (char *)&mn7min_1.up, (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L4000;
/*                                       . . . . . . . . . . show input, */
/*                                                batch, or interactive */
L1100:
    ioin__1.inerr = 0;
    ioin__1.inunit = mn7iou_1.isysrd;
    ioin__1.infile = 0;
    ioin__1.inex = 0;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = &lname;
    ioin__1.innamlen = 64;
    ioin__1.inname = cfname;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    s_copy(cmode, "BATCH MODE      ", (ftnlen)16, (ftnlen)16);
    if (mn7flg_1.isw[5] == 1) {
	s_copy(cmode, "INTERACTIVE MODE", (ftnlen)16, (ftnlen)16);
    }
    if (! lname) {
	s_copy(cfname, "unknown", (ftnlen)64, (ftnlen)7);
    }
    io___1103.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1103);
    do_fio(&c__1, cmode, (ftnlen)16);
    do_fio(&c__1, (char *)&mn7iou_1.isysrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, cfname, (ftnlen)64);
    e_wsfe();
    goto L4000;
/*                                       . . . . . . . . . . show width */
L1110:
    io___1104.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1104);
    do_fio(&c__1, (char *)&mn7iou_1.npagwd, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L4000;
/*                                       . . . . . . . . . . show lines */
L1120:
    io___1105.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1105);
    do_fio(&c__1, (char *)&mn7iou_1.npagln, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L4000;
/*                                       . . . . . . .show nowarn, warn */
L1130:
    s_copy(cwarn, "SUPPRESSED", (ftnlen)10, (ftnlen)10);
    if (mn7log_1.lwarn) {
	s_copy(cwarn, "REPORTED  ", (ftnlen)10, (ftnlen)10);
    }
    io___1107.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1107);
    do_fio(&c__1, cwarn, (ftnlen)10);
    e_wsfe();
    if (! mn7log_1.lwarn) {
	mnwarn_("W", "SHO", "SHO", (ftnlen)1, (ftnlen)3, (ftnlen)3);
    }
    goto L4000;
/*                                      . . . . . . . . . . show random */
L1150:
    val = (float)0.;
    mnrn15_(&val, &igrain);
    ikseed = igrain;
    io___1110.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1110);
    do_fio(&c__1, (char *)&ikseed, (ftnlen)sizeof(integer));
    e_wsfe();
    val = (float)3.;
    iseed = ikseed;
    mnrn15_(&val, &iseed);
    goto L4000;
/*                                        . . . . . . . . . show title */
L1160:
    io___1112.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1112);
    do_fio(&c__1, " TITLE OF CURRENT TASK IS:", (ftnlen)26);
    do_fio(&c__1, mn7tit_1.ctitl, (ftnlen)50);
    e_wsfe();
    goto L4000;
/*                                        . . . . . . . show strategy */
L1170:
    io___1113.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1113);
    do_fio(&c__1, " ALLOWED STRATEGIES ARE:", (ftnlen)24);
    e_wsfe();
    io___1114.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1114);
    do_fio(&c__3, cstrat, (ftnlen)44);
    e_wsfe();
L1172:
    io___1115.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1115);
    do_fio(&c__1, cstrat + mn7cnv_1.istrat * 44, (ftnlen)44);
    e_wsfe();
    goto L4000;
/*                                          . . . . . show eigenvalues */
L1180:
    iswsav = mn7flg_1.isw[4];
    mn7flg_1.isw[4] = 3;
    if (mn7flg_1.isw[1] < 1) {
	io___1117.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1117);
	do_fio(&c__1, mn7tit_1.covmes, (ftnlen)22);
	e_wsfe();
    } else {
	mnpsdf_();
    }
    mn7flg_1.isw[4] = iswsav;
    goto L4000;
/*                                            . . . . . show page throw */
L1190:
    io___1118.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1118);
    do_fio(&c__1, " PAGE THROW CARRIAGE CONTROL =", (ftnlen)30);
    do_fio(&c__1, (char *)&mn7iou_1.newpag, (ftnlen)sizeof(integer));
    e_wsfe();
    if (mn7iou_1.newpag == 0) {
	io___1119.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1119);
	do_fio(&c__1, " NO PAGE THROWS IN MINUIT OUTPUT", (ftnlen)32);
	e_wsfe();
    }
    goto L4000;
/*                                        . . . . . . show minos errors */
L1200:
    i__1 = mn7npr_1.npar;
    for (ii = 1; ii <= i__1; ++ii) {
	if (mn7err_1.erp[ii - 1] > 0. || mn7err_1.ern[ii - 1] < 0.) {
	    goto L1204;
	}
/* L1202: */
    }
    io___1121.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1121);
    do_fio(&c__1, "       THERE ARE NO MINOS ERRORS CURRENTLY VALID.", (
	    ftnlen)49);
    e_wsfe();
    goto L4000;
L1204:
    mnprin_(&c__4, &mn7min_1.amin);
    goto L4000;
/*                                        . . . . . . . . . show epsmac */
L1210:
    io___1122.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1122);
    do_fio(&c__1, " FLOATING-POINT NUMBERS ASSUMED ACCURATE TO", (ftnlen)43);
    do_fio(&c__1, (char *)&mn7cns_1.epsmac, (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L4000;
/*                                        . . . . . . show outputfiles */
L1220:
    io___1123.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1123);
    do_fio(&c__1, "  MINUIT PRIMARY OUTPUT TO UNIT", (ftnlen)31);
    do_fio(&c__1, (char *)&mn7iou_1.isyswr, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L4000;
/*                                        . . . . . . show version */
L1250:
    io___1124.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1124);
    do_fio(&c__1, " THIS IS MINUIT VERSION:", (ftnlen)24);
    do_fio(&c__1, mn7tit_1.cvrsn, (ftnlen)6);
    e_wsfe();
    goto L4000;
/*                                        . . . . . . show nodebug, debug */
L1270:
    for (id = 0; id <= 6; ++id) {
	s_copy(copt, "OFF", (ftnlen)3, (ftnlen)3);
	if (mn7flg_1.idbg[id] >= 1) {
	    s_copy(copt, "ON ", (ftnlen)3, (ftnlen)3);
	}
/* L1285: */
	io___1126.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1126);
	do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	do_fio(&c__1, copt, (ftnlen)3);
	do_fio(&c__1, cdbopt + id * 40, (ftnlen)40);
	e_wsfe();
    }
    if (! mn7log_1.lrepor) {
	mnwarn_("D", "SHO", "SHO", (ftnlen)1, (ftnlen)3, (ftnlen)3);
    }
    goto L4000;
/*                                        . . . . . . . . . . show show */
L1290:
    s_copy(ckind, "SHOW", (ftnlen)4, (ftnlen)4);
    goto L2100;
/*                                        . . . . . . . . . . show set */
L1300:
    s_copy(ckind, "SET ", (ftnlen)4, (ftnlen)4);
    goto L2100;
/*                ----------------------------------------------------- */
/*                              UNKNOWN COMMAND */
L1900:
    io___1127.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1127);
    do_fio(&c__1, mn7tit_1.cword, (ftnlen)20);
    e_wsfe();
    goto L2100;
/*                ----------------------------------------------------- */
/*                    HELP SHOW,  HELP SET,  SHOW SET, or SHOW SHOW */
L2000:
    s_copy(ckind, "SET ", (ftnlen)4, (ftnlen)4);
    if (i_indx(mn7tit_1.cword + 3, "SHO", (ftnlen)7, (ftnlen)3) > 0) {
	s_copy(ckind, "SHOW", (ftnlen)4, (ftnlen)4);
    }
L2100:
    io___1128.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1128);
    do_fio(&c__1, ckind, (ftnlen)4);
    do_fio(&c__1, ckind, (ftnlen)4);
    i__1 = nname;
    for (kk = 1; kk <= i__1; ++kk) {
	do_fio(&c__1, cname + (kk - 1) * 10, (ftnlen)10);
    }
    e_wsfe();
    goto L4000;
/*                ----------------------------------------------------- */
/*                               ILLEGAL COMMAND */
L3000:
    io___1130.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1130);
    e_wsfe();
L4000:
    return 0;
} /* mnset_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnseti_(char *tit, ftnlen tit_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       Called by user to set or change title of current task. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    s_copy(mn7tit_1.ctitl, tit, (ftnlen)50, tit_len);
    return 0;
} /* mnseti_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:54  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnsimp_(S_fp fcn, U_fp futil)
{
    /* Initialized data */

    static doublereal alpha = 1.;
    static doublereal beta = .5;
    static doublereal gamma = 2.;
    static doublereal rhomin = 4.;
    static doublereal rhomax = 8.;

    /* Format strings */
    static char fmt_100[] = "(\002 START SIMPLEX MINIMIZATION.    CONVERGENC\
E WHEN EDM .LT.\002,e10.2)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static doublereal f;
    static integer i__, j, k;
    static doublereal y[101], y1, y2;
    static integer kg, jh, nf;
    static doublereal pb;
    static integer jl;
    static doublereal wg;
    static integer ns;
    static doublereal rho, sig2, rho1, rho2, dmin__, dxdi;
    static integer npfn;
    static doublereal yrho, ynpp1, aming;
    static integer jhold, ncycl;
    static doublereal ypbar;
    static integer nparx;
    static doublereal bestx, ystar, ystst;
    static integer nparp1;
    static doublereal absmin;
    extern /* Subroutine */ int mnamin_(S_fp, U_fp), mndxdi_(doublereal *, 
	    integer *, doublereal *), mninex_(doublereal *), mnrazz_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    mnprin_(integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___1142 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___1169 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1170 = { 0, 0, 0, "(A)", 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Performs a minimization using the simplex method of Nelder */
/* C        and Mead (ref. -- Comp. J. 7,308 (1965)). */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    if (mn7npr_1.npar <= 0) {
	return 0;
    }
    if (mn7min_1.amin == mn7cns_1.undefi) {
	mnamin_((S_fp)fcn, (U_fp)futil);
    }
    s_copy(mn7tit_1.cfrom, "SIMPLEX ", (ftnlen)8, (ftnlen)8);
    mn7cnv_1.nfcnfr = mn7cnv_1.nfcn;
    s_copy(mn7tit_1.cstatu, "UNCHANGED ", (ftnlen)10, (ftnlen)10);
    npfn = mn7cnv_1.nfcn;
    nparp1 = mn7npr_1.npar + 1;
    nparx = mn7npr_1.npar;
    rho1 = alpha + (float)1.;
    rho2 = rho1 + alpha * gamma;
    wg = (float)1. / (real) mn7npr_1.npar;
    if (mn7flg_1.isw[4] >= 0) {
	io___1142.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1142);
	do_fio(&c__1, (char *)&mn7min_1.epsi, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mn7int_1.dirin[i__ - 1] = mn7err_1.werr[i__ - 1];
	mndxdi_(&mn7int_1.x[i__ - 1], &i__, &dxdi);
	if (dxdi != 0.) {
	    mn7int_1.dirin[i__ - 1] = mn7err_1.werr[i__ - 1] / dxdi;
	}
	dmin__ = mn7cns_1.epsma2 * (d__1 = mn7int_1.x[i__ - 1], abs(d__1));
	if (mn7int_1.dirin[i__ - 1] < dmin__) {
	    mn7int_1.dirin[i__ - 1] = dmin__;
	}
/* L2: */
    }
/* **       choose the initial simplex using single-parameter searches */
L1:
    ynpp1 = mn7min_1.amin;
    jl = nparp1;
    y[nparp1 - 1] = mn7min_1.amin;
    absmin = mn7min_1.amin;
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aming = mn7min_1.amin;
	mn7sim_1.pbar[i__ - 1] = mn7int_1.x[i__ - 1];
	bestx = mn7int_1.x[i__ - 1];
	kg = 0;
	ns = 0;
	nf = 0;
L4:
	mn7int_1.x[i__ - 1] = bestx + mn7int_1.dirin[i__ - 1];
	mninex_(mn7int_1.x);
	(*fcn)(&nparx, mn7der_1.gin, &f, mn7ext_1.u, &c__4, (U_fp)futil);
	++mn7cnv_1.nfcn;
	if (f < aming) {
	    goto L6;
	}
/*         failure */
	if (kg == 1) {
	    goto L8;
	}
	kg = -1;
	++nf;
	mn7int_1.dirin[i__ - 1] *= (float)-.4;
	if (nf < 3) {
	    goto L4;
	}
/*         stop after three failures */
	bestx = mn7int_1.x[i__ - 1];
	mn7int_1.dirin[i__ - 1] *= (float)3.;
	aming = f;
	goto L8;

/*         success */
L6:
	bestx = mn7int_1.x[i__ - 1];
	mn7int_1.dirin[i__ - 1] *= (float)3.;
	aming = f;
	s_copy(mn7tit_1.cstatu, "PROGRESS  ", (ftnlen)10, (ftnlen)10);
	kg = 1;
	++ns;
	if (ns < 6) {
	    goto L4;
	}

/*         3 failures or 6 successes or */
/*         local minimum found in ith direction */
L8:
	y[i__ - 1] = aming;
	if (aming < absmin) {
	    jl = i__;
	}
	if (aming < absmin) {
	    absmin = aming;
	}
	mn7int_1.x[i__ - 1] = bestx;
	i__2 = mn7npr_1.npar;
	for (k = 1; k <= i__2; ++k) {
/* L9: */
	    mn7sim_1.p[k + i__ * 100 - 101] = mn7int_1.x[k - 1];
	}
/* L10: */
    }
    jh = nparp1;
    mn7min_1.amin = y[jl - 1];
    mnrazz_(&ynpp1, mn7sim_1.pbar, y, &jh, &jl);
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	mn7int_1.x[i__ - 1] = mn7sim_1.p[i__ + jl * 100 - 101];
    }
    mninex_(mn7int_1.x);
    if (mn7flg_1.isw[4] >= 1) {
	mnprin_(&c__5, &mn7min_1.amin);
    }
    mn7min_1.edm = mn7cns_1.bigedm;
    sig2 = mn7min_1.edm;
    ncycl = 0;
/*                                        . . . . .  start main loop */
L50:
    if (sig2 < mn7min_1.epsi && mn7min_1.edm < mn7min_1.epsi) {
	goto L76;
    }
    sig2 = mn7min_1.edm;
    if (mn7cnv_1.nfcn - npfn > mn7cnv_1.nfcnmx) {
	goto L78;
    }
/*         calculate new point * by reflection */
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pb = (float)0.;
	i__2 = nparp1;
	for (j = 1; j <= i__2; ++j) {
/* L59: */
	    pb += wg * mn7sim_1.p[i__ + j * 100 - 101];
	}
	mn7sim_1.pbar[i__ - 1] = pb - wg * mn7sim_1.p[i__ + jh * 100 - 101];
/* L60: */
	mn7sim_1.pstar[i__ - 1] = (alpha + (float)1.) * mn7sim_1.pbar[i__ - 1]
		 - alpha * mn7sim_1.p[i__ + jh * 100 - 101];
    }
    mninex_(mn7sim_1.pstar);
    (*fcn)(&nparx, mn7der_1.gin, &ystar, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    if (ystar >= mn7min_1.amin) {
	goto L70;
    }
/*         point * better than jl, calculate new point ** */
    s_copy(mn7tit_1.cstatu, "PROGRESS  ", (ftnlen)10, (ftnlen)10);
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L61: */
	mn7sim_1.pstst[i__ - 1] = gamma * mn7sim_1.pstar[i__ - 1] + ((float)
		1. - gamma) * mn7sim_1.pbar[i__ - 1];
    }
    mninex_(mn7sim_1.pstst);
    (*fcn)(&nparx, mn7der_1.gin, &ystst, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
/*         try a parabola through ph, pstar, pstst.  min = prho */
    y1 = (ystar - y[jh - 1]) * rho2;
    y2 = (ystst - y[jh - 1]) * rho1;
    rho = (rho2 * y1 - rho1 * y2) * (float).5 / (y1 - y2);
    if (rho < rhomin) {
	goto L66;
    }
    if (rho > rhomax) {
	rho = rhomax;
    }
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L64: */
	mn7sim_1.prho[i__ - 1] = rho * mn7sim_1.pbar[i__ - 1] + ((float)1. - 
		rho) * mn7sim_1.p[i__ + jh * 100 - 101];
    }
    mninex_(mn7sim_1.prho);
    (*fcn)(&nparx, mn7der_1.gin, &yrho, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    if (yrho < mn7min_1.amin) {
	s_copy(mn7tit_1.cstatu, "PROGRESS  ", (ftnlen)10, (ftnlen)10);
    }
    if (yrho < y[jl - 1] && yrho < ystst) {
	goto L65;
    }
    if (ystst < y[jl - 1]) {
	goto L67;
    }
    if (yrho > y[jl - 1]) {
	goto L66;
    }
/*         accept minimum point of parabola, PRHO */
L65:
    mnrazz_(&yrho, mn7sim_1.prho, y, &jh, &jl);
    goto L68;
L66:
    if (ystst < y[jl - 1]) {
	goto L67;
    }
    mnrazz_(&ystar, mn7sim_1.pstar, y, &jh, &jl);
    goto L68;
L67:
    mnrazz_(&ystst, mn7sim_1.pstst, y, &jh, &jl);
L68:
    ++ncycl;
    if (mn7flg_1.isw[4] < 2) {
	goto L50;
    }
    if (mn7flg_1.isw[4] >= 3 || ncycl % 10 == 0) {
	mnprin_(&c__5, &mn7min_1.amin);
    }
    goto L50;
/*         point * is not as good as jl */
L70:
    if (ystar >= y[jh - 1]) {
	goto L73;
    }
    jhold = jh;
    mnrazz_(&ystar, mn7sim_1.pstar, y, &jh, &jl);
    if (jhold != jh) {
	goto L50;
    }
/*         calculate new point ** */
L73:
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L74: */
	mn7sim_1.pstst[i__ - 1] = beta * mn7sim_1.p[i__ + jh * 100 - 101] + ((
		float)1. - beta) * mn7sim_1.pbar[i__ - 1];
    }
    mninex_(mn7sim_1.pstst);
    (*fcn)(&nparx, mn7der_1.gin, &ystst, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    if (ystst > y[jh - 1]) {
	goto L1;
    }
/*     point ** is better than jh */
    if (ystst < mn7min_1.amin) {
	s_copy(mn7tit_1.cstatu, "PROGRESS  ", (ftnlen)10, (ftnlen)10);
    }
    if (ystst < mn7min_1.amin) {
	goto L67;
    }
    mnrazz_(&ystst, mn7sim_1.pstst, y, &jh, &jl);
    goto L50;
/*                                        . . . . . .  end main loop */
L76:
    if (mn7flg_1.isw[4] >= 0) {
	io___1169.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1169);
	do_fio(&c__1, " SIMPLEX MINIMIZATION HAS CONVERGED.", (ftnlen)36);
	e_wsfe();
    }
    mn7flg_1.isw[3] = 1;
    goto L80;
L78:
    if (mn7flg_1.isw[4] >= 0) {
	io___1170.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1170);
	do_fio(&c__1, " SIMPLEX TERMINATES WITHOUT CONVERGENCE.", (ftnlen)40);
	e_wsfe();
    }
    s_copy(mn7tit_1.cstatu, "CALL LIMIT", (ftnlen)10, (ftnlen)10);
    mn7flg_1.isw[3] = -1;
    mn7flg_1.isw[0] = 1;
L80:
    i__1 = mn7npr_1.npar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pb = (float)0.;
	i__2 = nparp1;
	for (j = 1; j <= i__2; ++j) {
/* L81: */
	    pb += wg * mn7sim_1.p[i__ + j * 100 - 101];
	}
/* L82: */
	mn7sim_1.pbar[i__ - 1] = pb - wg * mn7sim_1.p[i__ + jh * 100 - 101];
    }
    mninex_(mn7sim_1.pbar);
    (*fcn)(&nparx, mn7der_1.gin, &ypbar, mn7ext_1.u, &c__4, (U_fp)futil);
    ++mn7cnv_1.nfcn;
    if (ypbar < mn7min_1.amin) {
	mnrazz_(&ypbar, mn7sim_1.pbar, y, &jh, &jl);
    }
    mninex_(mn7int_1.x);
    if (mn7cnv_1.nfcnmx + npfn - mn7cnv_1.nfcn < mn7npr_1.npar * 3) {
	goto L90;
    }
    if (mn7min_1.edm > mn7min_1.epsi * (float)2.) {
	goto L1;
    }
L90:
    if (mn7flg_1.isw[4] >= 0) {
	mnprin_(&c__5, &mn7min_1.amin);
    }
    return 0;
} /* mnsimp_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni */
/* Minuit */


/* Subroutine */ int mnstat_(doublereal *fmin, doublereal *fedm, doublereal *
	errdef, integer *npari, integer *nparx, integer *istat)
{

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C       User-called */
/* C       Provides the user with information concerning the current status */
/* C          of the current minimization. Namely, it returns: */
/* C        FMIN: the best function value found so far */
/* C        FEDM: the estimated vertical distance remaining to minimum */
/* C        ERRDEF: the value of UP defining parameter uncertainties */
/* C        NPARI: the number of currently variable parameters */
/* C        NPARX: the highest (external) parameter number defined by user */
/* C        ISTAT: a status integer indicating how good is the covariance */
/* C           matrix:  0= not calculated at all */
/* C                    1= approximation only, not accurate */
/* C                    2= full matrix, but forced positive-definite */
/* C                    3= full accurate covariance matrix */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    *fmin = mn7min_1.amin;
    *fedm = mn7min_1.edm;
    *errdef = mn7min_1.up;
    *npari = mn7npr_1.npar;
    *nparx = mn7npr_1.nu;
    *istat = mn7flg_1.isw[1];
    if (mn7min_1.edm == mn7cns_1.bigedm) {
	*fedm = mn7min_1.up;
    }
    if (mn7min_1.amin == mn7cns_1.undefi) {
	*fmin = (float)0.;
	*fedm = mn7min_1.up;
	*istat = 0;
    }
    return 0;
} /* mnstat_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */


/* Subroutine */ int mnstin_(char *crdbuf, integer *ierr, ftnlen crdbuf_len)
{
    /* Format strings */
    static char fmt_132[] = "(\002 UNIT\002,i3,\002 ALREADY OPENED WITH NA\
ME:\002,a/\002                 NEW NAME IGNORED:\002,a)";
    static char fmt_135[] = "(\002 UNIT\002,i3,\002 IS NOT OPENED.\002)";
    static char fmt_137[] = "(\002 SHOULD UNIT\002,i3,\002 BE REWOUND?\002)";
    static char fmt_290[] = "(\002 INPUT WILL NOW BE READ IN \002,a,\002 FRO\
M UNIT NO.\002,i3/\002 FILENAME: \002,a)";
    static char fmt_601[] = "(\002 SYSTEM IS UNABLE TO OPEN FILE:\002,a)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer i_indx(char *, char *, ftnlen, ftnlen), i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_rsfi(icilist *), e_rsfi(), f_inqu(inlist *), s_rsfe(cilist *), 
	    e_rsfe(), f_open(olist *), f_rew(alist *);

    /* Local variables */
    static integer ic, ic1, ic2, lend, icol;
    static char cmode[16];
    static logical lname, lopen;
    static char cunit__[10];
    static doublereal funit;
    static integer iunit;
    static char cfname[64], cgname[64];
    static logical noname;
    static char canswr[1];
    static logical lrewin;
    extern logical mnunpt_(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1180 = { 0, 0, 0, "(A,A)", 0 };
    static icilist io___1181 = { 1, cunit__, 0, "(BN,F10.0)", 10, 1 };
    static cilist io___1185 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___1189 = { 0, 0, 0, fmt_132, 0 };
    static cilist io___1190 = { 0, 0, 0, fmt_135, 0 };
    static cilist io___1191 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1192 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1193 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1194 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1195 = { 0, 0, 0, fmt_137, 0 };
    static cilist io___1196 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1198 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___1199 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1201 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___1202 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1203 = { 0, 0, 0, "(A,A)", 0 };
    static cilist io___1204 = { 0, 0, 0, fmt_601, 0 };



/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C Called from MNREAD. */
/* C Implements the SET INPUT command to change input units. */
/* C If command is: 'SET INPUT'   'SET INPUT 0'   or  '*EOF', */
/* C                 or 'SET INPUT , ,  ', */
/* C                reverts to previous input unit number,if any. */
/* C */
/* C      If it is: 'SET INPUT n'  or  'SET INPUT n filename', */
/* C                changes to new input file, added to stack */
/* C */
/* C      IERR = 0: reading terminated normally */
/* C             2: end-of-data on primary input file */
/* C             3: unrecoverable read error */
/* C             4: unable to process request */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    noname = TRUE_;
    *ierr = 0;
    if (i_indx(crdbuf, "*EOF", crdbuf_len, (ftnlen)4) == 1) {
	goto L190;
    }
    if (i_indx(crdbuf, "*eof", crdbuf_len, (ftnlen)4) == 1) {
	goto L190;
    }
    lend = i_len(crdbuf, crdbuf_len);
/*                               look for end of SET INPUT command */
    i__1 = lend;
    for (ic = 8; ic <= i__1; ++ic) {
	if (*(unsigned char *)&crdbuf[ic - 1] == ' ') {
	    goto L25;
	}
	if (*(unsigned char *)&crdbuf[ic - 1] == ',') {
	    goto L53;
	}
/* L20: */
    }
    goto L200;
L25:
/*         look for end of separator between command and first argument */
    icol = ic + 1;
    i__1 = lend;
    for (ic = icol; ic <= i__1; ++ic) {
	if (*(unsigned char *)&crdbuf[ic - 1] == ' ') {
	    goto L50;
	}
	if (*(unsigned char *)&crdbuf[ic - 1] == ',') {
	    goto L53;
	}
	goto L55;
L50:
	;
    }
    goto L200;
L53:
    ++ic;
L55:
    ic1 = ic;
/*                      see if "REWIND" was requested in command */
    lrewin = FALSE_;
    if (i_indx(crdbuf, "REW", ic1, (ftnlen)3) > 5) {
	lrewin = TRUE_;
    }
    if (i_indx(crdbuf, "rew", ic1, (ftnlen)3) > 5) {
	lrewin = TRUE_;
    }
/*                      first argument begins in or after col IC1 */
    i__1 = lend;
    for (ic = ic1; ic <= i__1; ++ic) {
	if (*(unsigned char *)&crdbuf[ic - 1] == ' ') {
	    goto L75;
	}
	if (*(unsigned char *)&crdbuf[ic - 1] == ',') {
	    goto L200;
	}
	goto L80;
L75:
	;
    }
    goto L200;
L80:
    ic1 = ic;
/*                        first argument really begins in col IC1 */
    i__1 = lend;
    for (ic = ic1 + 1; ic <= i__1; ++ic) {
	if (*(unsigned char *)&crdbuf[ic - 1] == ' ') {
	    goto L108;
	}
	if (*(unsigned char *)&crdbuf[ic - 1] == ',') {
	    goto L108;
	}
/* L100: */
    }
    ic = lend + 1;
L108:
    ic2 = ic - 1;
/*                            end of first argument is in col IC2 */
/* L110: */
    s_copy(cunit__, crdbuf + (ic1 - 1), (ftnlen)10, ic2 - (ic1 - 1));
    io___1180.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1180);
    do_fio(&c__1, " UNIT NO. :", (ftnlen)11);
    do_fio(&c__1, cunit__, (ftnlen)10);
    e_wsfe();
    i__1 = s_rsfi(&io___1181);
    if (i__1 != 0) {
	goto L500;
    }
    i__1 = do_fio(&c__1, (char *)&funit, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L500;
    }
    i__1 = e_rsfi();
    if (i__1 != 0) {
	goto L500;
    }
    iunit = (integer) funit;
    if (iunit == 0) {
	goto L200;
    }
/*                             skip blanks and commas, find file name */
    i__1 = lend;
    for (ic = ic2 + 1; ic <= i__1; ++ic) {
	if (*(unsigned char *)&crdbuf[ic - 1] == ' ') {
	    goto L120;
	}
	if (*(unsigned char *)&crdbuf[ic - 1] == ',') {
	    goto L120;
	}
	goto L130;
L120:
	;
    }
    goto L131;
L130:
    s_copy(cfname, crdbuf + (ic - 1), (ftnlen)64, lend - (ic - 1));
    noname = FALSE_;
    io___1185.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1185);
    do_fio(&c__1, " FILE NAME IS:", (ftnlen)14);
    do_fio(&c__1, cfname, (ftnlen)64);
    e_wsfe();
/*              ask if file exists, if not ask for name and open it */
L131:
    ioin__1.inerr = 0;
    ioin__1.inunit = iunit;
    ioin__1.infile = 0;
    ioin__1.inex = 0;
    ioin__1.inopen = &lopen;
    ioin__1.innum = 0;
    ioin__1.innamed = &lname;
    ioin__1.innamlen = 64;
    ioin__1.inname = cgname;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (lopen) {
	if (noname) {
	    goto L136;
	} else {
	    if (! lname) {
		s_copy(cgname, "unknown", (ftnlen)64, (ftnlen)7);
	    }
	    io___1189.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1189);
	    do_fio(&c__1, (char *)&iunit, (ftnlen)sizeof(integer));
	    do_fio(&c__1, cgname, (ftnlen)64);
	    do_fio(&c__1, cfname, (ftnlen)64);
	    e_wsfe();
	}
    } else {
/*                new file, open it */
	io___1190.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1190);
	do_fio(&c__1, (char *)&iunit, (ftnlen)sizeof(integer));
	e_wsfe();
	if (noname) {
	    io___1191.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1191);
	    do_fio(&c__1, " NO FILE NAME GIVEN IN COMMAND.", (ftnlen)31);
	    e_wsfe();
	    if (mn7flg_1.isw[5] < 1) {
		goto L800;
	    }
	    io___1192.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1192);
	    do_fio(&c__1, " PLEASE GIVE FILE NAME:", (ftnlen)23);
	    e_wsfe();
	    io___1193.ciunit = mn7iou_1.isysrd;
	    s_rsfe(&io___1193);
	    do_fio(&c__1, cfname, (ftnlen)64);
	    e_rsfe();
	}
	o__1.oerr = 1;
	o__1.ounit = iunit;
	o__1.ofnmlen = 64;
	o__1.ofnm = cfname;
	o__1.orl = 0;
	o__1.osta = "OLD";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L600;
	}
	io___1194.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1194);
	do_fio(&c__1, " FILE OPENED SUCCESSFULLY.", (ftnlen)26);
	e_wsfe();
    }
/*                                     . .   file is correctly opened */
L136:
    if (lrewin) {
	goto L150;
    }
    if (mn7flg_1.isw[5] < 1) {
	goto L300;
    }
    io___1195.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1195);
    do_fio(&c__1, (char *)&iunit, (ftnlen)sizeof(integer));
    e_wsfe();
    io___1196.ciunit = mn7iou_1.isysrd;
    s_rsfe(&io___1196);
    do_fio(&c__1, canswr, (ftnlen)1);
    e_rsfe();
    if (*(unsigned char *)canswr != 'Y' && *(unsigned char *)canswr != 'y') {
	goto L300;
    }
L150:
    al__1.aerr = 0;
    al__1.aunit = iunit;
    f_rew(&al__1);
    goto L300;
/*                      *EOF */
L190:
    if (mn7io2_1.nstkrd == 0) {
	*ierr = 2;
	goto L900;
    }
/*                      revert to previous input file */
L200:
    if (mn7io2_1.nstkrd == 0) {
	io___1198.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1198);
	do_fio(&c__1, " COMMAND IGNORED:", (ftnlen)17);
	do_fio(&c__1, crdbuf, crdbuf_len);
	e_wsfe();
	io___1199.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1199);
	do_fio(&c__1, " ALREADY READING FROM PRIMARY INPUT", (ftnlen)35);
	e_wsfe();
    } else {
	mn7iou_1.isysrd = mn7io2_1.istkrd[mn7io2_1.nstkrd - 1];
	--mn7io2_1.nstkrd;
	if (mn7io2_1.nstkrd == 0) {
	    mn7flg_1.isw[5] = abs(mn7flg_1.isw[5]);
	}
	if (mn7flg_1.isw[4] >= 0) {
	    ioin__1.inerr = 0;
	    ioin__1.inunit = mn7iou_1.isysrd;
	    ioin__1.infile = 0;
	    ioin__1.inex = 0;
	    ioin__1.inopen = 0;
	    ioin__1.innum = 0;
	    ioin__1.innamed = &lname;
	    ioin__1.innamlen = 64;
	    ioin__1.inname = cfname;
	    ioin__1.inacc = 0;
	    ioin__1.inseq = 0;
	    ioin__1.indir = 0;
	    ioin__1.infmt = 0;
	    ioin__1.inform = 0;
	    ioin__1.inunf = 0;
	    ioin__1.inrecl = 0;
	    ioin__1.innrec = 0;
	    ioin__1.inblank = 0;
	    f_inqu(&ioin__1);
	    s_copy(cmode, "BATCH MODE      ", (ftnlen)16, (ftnlen)16);
	    if (mn7flg_1.isw[5] == 1) {
		s_copy(cmode, "INTERACTIVE MODE", (ftnlen)16, (ftnlen)16);
	    }
	    if (! lname) {
		s_copy(cfname, "unknown", (ftnlen)64, (ftnlen)7);
	    }
	    if (mnunpt_(cfname, (ftnlen)64)) {
		s_copy(cfname, "unprintable", (ftnlen)64, (ftnlen)11);
	    }
	    io___1201.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1201);
	    do_fio(&c__1, cmode, (ftnlen)16);
	    do_fio(&c__1, (char *)&mn7iou_1.isysrd, (ftnlen)sizeof(integer));
	    do_fio(&c__1, cfname, (ftnlen)64);
	    e_wsfe();
	}
    }
    goto L900;
/*                      switch to new input file, add to stack */
L300:
    if (mn7io2_1.nstkrd >= 10) {
	io___1202.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1202);
	do_fio(&c__1, " INPUT FILE STACK SIZE EXCEEDED.", (ftnlen)32);
	e_wsfe();
	goto L800;
    }
    ++mn7io2_1.nstkrd;
    mn7io2_1.istkrd[mn7io2_1.nstkrd - 1] = mn7iou_1.isysrd;
    mn7iou_1.isysrd = iunit;
/*                   ISW(6) = 0 for batch, =1 for interactive, and */
/*                      =-1 for originally interactive temporarily batch */
    if (mn7flg_1.isw[5] == 1) {
	mn7flg_1.isw[5] = -1;
    }
    goto L900;
/*                      format error */
L500:
    io___1203.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1203);
    do_fio(&c__1, " CANNOT READ FOLLOWING AS INTEGER:", (ftnlen)34);
    do_fio(&c__1, cunit__, (ftnlen)10);
    e_wsfe();
    goto L800;
L600:
    io___1204.ciunit = mn7iou_1.isyswr;
    s_wsfe(&io___1204);
    do_fio(&c__1, cfname, (ftnlen)64);
    e_wsfe();
/*                      serious error */
L800:
    *ierr = 3;
L900:
    return 0;
} /* mnstin_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */


/* Subroutine */ int mntiny_(doublereal *epsp1, doublereal *epsbak)
{

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        Compares its argument with the value 1.0, and returns */
/* C        the value .TRUE. if they are equal.  To find EPSMAC */
/* C        safely by foiling the Fortran optimizer */
/* C */
    *epsbak = *epsp1 - 1.;
    return 0;
} /* mntiny_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */


logical mnunpt_(char *cfname, ftnlen cfname_len)
{
    /* System generated locals */
    integer i__1;
    logical ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, l, ic;
    static char cpt[80];

/*           is .TRUE. if CFNAME contains unprintable characters. */
    s_copy(cpt, " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456\
7890./;:[]$%*_!@#&+()", (ftnlen)80, (ftnlen)80);
    ret_val = FALSE_;
    l = i_len(cfname, cfname_len);
    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (ic = 1; ic <= 80; ++ic) {
	    if (*(unsigned char *)&cfname[i__ - 1] == *(unsigned char *)&cpt[
		    ic - 1]) {
		goto L100;
	    }
/* L50: */
	}
	ret_val = TRUE_;
	goto L150;
L100:
	;
    }
L150:
    return ret_val;
} /* mnunpt_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */


/* Subroutine */ int mnvers_(char *cv, ftnlen cv_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C         Returns the Minuit version in CV, char*6 */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    s_copy(cv, mn7tit_1.cvrsn, cv_len, (ftnlen)6);
    return 0;
} /* mnvers_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.2  1996/03/15 18:02:54  james */
/*     Modified Files: */
/* mnderi.F eliminate possible division by zero */
/* mnexcm.F suppress print on STOP when print flag=-1 */
/*          set FVAL3 to flag if FCN already called with IFLAG=3 */
/* mninit.F set version 96.03 */
/* mnlims.F remove arguments, not needed */
/* mnmigr.F VLEN -> LENV in debug print statement */
/* mnparm.F move call to MNRSET to after NPAR redefined, to zero all */
/* mnpsdf.F eliminate possible division by zero */
/* mnscan.F suppress printout when print flag =-1 */
/* mnset.F  remove arguments in call to MNLIMS */
/* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum */
/* mnvert.F eliminate possible division by zero */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */


/* Subroutine */ int mnvert_(doublereal *a, integer *l, integer *m, integer *
	n, integer *ifail)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal q[100], s[100], si, pp[100];
    static integer km1, kp1;


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        inverts a symmetric matrix.   matrix is first scaled to */
/* C        have all ones on the diagonal (equivalent to change of units) */
/* C        but no pivoting is done since matrix is positive-definite. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


    /* Parameter adjustments */
    a_dim1 = *l;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    *ifail = 0;
    if (*n < 1) {
	goto L100;
    }
    if (*n > mn7npr_1.maxint) {
	goto L100;
    }
/*                   scale matrix by sqrt of diag elements */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	si = a[i__ + i__ * a_dim1];
	if (si <= 0.) {
	    goto L100;
	} else {
	    goto L8;
	}
L8:
	s[i__ - 1] = (float)1. / sqrt(si);
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L20: */
	    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] * s[i__ - 1] * s[j - 1];
	}
    }
/*                                        . . . start main loop . . . . */
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	k = i__;
/*                   preparation for elimination step1 */
	if (a[k + k * a_dim1] == 0.) {
	    goto L100;
	}
	q[k - 1] = (float)1. / a[k + k * a_dim1];
	pp[k - 1] = (float)1.;
	a[k + k * a_dim1] = (float)0.;
	kp1 = k + 1;
	km1 = k - 1;
	if (km1 < 0) {
	    goto L100;
	} else if (km1 == 0) {
	    goto L50;
	} else {
	    goto L40;
	}
L40:
	i__1 = km1;
	for (j = 1; j <= i__1; ++j) {
	    pp[j - 1] = a[j + k * a_dim1];
	    q[j - 1] = a[j + k * a_dim1] * q[k - 1];
/* L49: */
	    a[j + k * a_dim1] = (float)0.;
	}
L50:
	if ((i__1 = k - *n) < 0) {
	    goto L51;
	} else if (i__1 == 0) {
	    goto L60;
	} else {
	    goto L100;
	}
L51:
	i__1 = *n;
	for (j = kp1; j <= i__1; ++j) {
	    pp[j - 1] = a[k + j * a_dim1];
	    q[j - 1] = -a[k + j * a_dim1] * q[k - 1];
/* L59: */
	    a[k + j * a_dim1] = (float)0.;
	}
/*                   elimination proper */
L60:
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__3 = *n;
	    for (k = j; k <= i__3; ++k) {
/* L65: */
		a[j + k * a_dim1] += pp[j - 1] * q[k - 1];
	    }
	}
    }
/*                   elements of left diagonal and unscaling */
    i__3 = *n;
    for (j = 1; j <= i__3; ++j) {
	i__1 = j;
	for (k = 1; k <= i__1; ++k) {
	    a[k + j * a_dim1] = a[k + j * a_dim1] * s[k - 1] * s[j - 1];
/* L70: */
	    a[j + k * a_dim1] = a[k + j * a_dim1];
	}
    }
    return 0;
/*                   failure return */
L100:
    *ifail = 1;
    return 0;
} /* mnvert_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */


/* Subroutine */ int mnwarn_(char *copt, char *corg, char *cmes, ftnlen 
	copt_len, ftnlen corg_len, ftnlen cmes_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, ic, nm;
    static char ctyp[7];
    static integer ityp;
    static char englsh[20];

    /* Fortran I/O blocks */
    static cilist io___1219 = { 0, 0, 0, "(A,A/A,A)", 0 };
    static cilist io___1220 = { 0, 0, 0, "(A,A/A,A)", 0 };
    static cilist io___1224 = { 0, 0, 0, "(/1X,I5,A,A,A,A/)", 0 };
    static cilist io___1226 = { 0, 0, 0, "(A,I2,A)", 0 };
    static cilist io___1227 = { 0, 0, 0, "(A)", 0 };
    static cilist io___1229 = { 0, 0, 0, "(1X,I6,1X,A,1X,A)", 0 };
    static cilist io___1230 = { 0, 0, 0, "(1H )", 0 };


/*     If COPT='W', CMES is a WARning message from CORG. */
/*     If COPT='D', CMES is a DEBug message from CORG. */
/*         If SET WARnings is in effect (the default), this routine */
/*             prints the warning message CMES coming from CORG. */
/*         If SET NOWarnings is in effect, the warning message is */
/*             stored in a circular buffer of length MAXMES. */
/*         If called with CORG=CMES='SHO', it prints the messages in */
/*             the circular buffer, FIFO, and empties the buffer. */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */



    if (s_cmp(corg, "SHO", (ftnlen)3, (ftnlen)3) == 0 && s_cmp(cmes, "SHO", (
	    ftnlen)3, (ftnlen)3) == 0) {
	goto L200;
    }
/*             Either print warning or put in buffer */
    if (*(unsigned char *)copt == 'W') {
	ityp = 1;
	if (mn7log_1.lwarn) {
	    io___1219.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1219);
	    do_fio(&c__1, " MINUIT WARNING IN ", (ftnlen)19);
	    do_fio(&c__1, corg, corg_len);
	    do_fio(&c__1, " ============== ", (ftnlen)16);
	    do_fio(&c__1, cmes, cmes_len);
	    e_wsfe();
	    return 0;
	}
    } else {
	ityp = 2;
	if (mn7log_1.lrepor) {
	    io___1220.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1220);
	    do_fio(&c__1, " MINUIT DEBUG FOR  ", (ftnlen)19);
	    do_fio(&c__1, corg, corg_len);
	    do_fio(&c__1, " ============== ", (ftnlen)16);
	    do_fio(&c__1, cmes, cmes_len);
	    e_wsfe();
	    return 0;
	}
    }
/*                 if appropriate flag is off, fill circular buffer */
    if (mn7cnv_1.nwrmes[ityp - 1] == 0) {
	mn7wri_1.icirc[ityp - 1] = 0;
    }
    ++mn7cnv_1.nwrmes[ityp - 1];
    ++mn7wri_1.icirc[ityp - 1];
    if (mn7wri_1.icirc[ityp - 1] > 10) {
	mn7wri_1.icirc[ityp - 1] = 1;
    }
    ic = mn7wri_1.icirc[ityp - 1];
    s_copy(mn7wrc_1.origin + (ic + ityp * 10 - 11) * 10, corg, (ftnlen)10, 
	    corg_len);
    s_copy(mn7wrc_1.warmes + (ic + ityp * 10 - 11) * 60, cmes, (ftnlen)60, 
	    cmes_len);
    mn7wri_1.nfcwar[ic + ityp * 10 - 11] = mn7cnv_1.nfcn;
    return 0;

/*             'SHO WARnings', ask if any suppressed mess in buffer */
L200:
    if (*(unsigned char *)copt == 'W') {
	ityp = 1;
	s_copy(ctyp, "WARNING", (ftnlen)7, (ftnlen)7);
    } else {
	ityp = 2;
	s_copy(ctyp, "*DEBUG*", (ftnlen)7, (ftnlen)7);
    }
    if (mn7cnv_1.nwrmes[ityp - 1] > 0) {
	s_copy(englsh, " WAS SUPPRESSED.  ", (ftnlen)20, (ftnlen)18);
	if (mn7cnv_1.nwrmes[ityp - 1] > 1) {
	    s_copy(englsh, "S WERE SUPPRESSED.", (ftnlen)20, (ftnlen)18);
	}
	io___1224.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1224);
	do_fio(&c__1, (char *)&mn7cnv_1.nwrmes[ityp - 1], (ftnlen)sizeof(
		integer));
	do_fio(&c__1, " MINUIT ", (ftnlen)8);
	do_fio(&c__1, ctyp, (ftnlen)7);
	do_fio(&c__1, " MESSAGE", (ftnlen)8);
	do_fio(&c__1, englsh, (ftnlen)20);
	e_wsfe();
	nm = mn7cnv_1.nwrmes[ityp - 1];
	ic = 0;
	if (nm > 10) {
	    io___1226.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1226);
	    do_fio(&c__1, " ONLY THE MOST RECENT ", (ftnlen)22);
	    do_fio(&c__1, (char *)&c__10, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " WILL BE LISTED BELOW.", (ftnlen)22);
	    e_wsfe();
	    nm = 10;
	    ic = mn7wri_1.icirc[ityp - 1];
	}
	io___1227.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1227);
	do_fio(&c__1, "  CALLS  ORIGIN         MESSAGE", (ftnlen)31);
	e_wsfe();
	i__1 = nm;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++ic;
	    if (ic > 10) {
		ic = 1;
	    }
	    io___1229.ciunit = mn7iou_1.isyswr;
	    s_wsfe(&io___1229);
	    do_fio(&c__1, (char *)&mn7wri_1.nfcwar[ic + ityp * 10 - 11], (
		    ftnlen)sizeof(integer));
	    do_fio(&c__1, mn7wrc_1.origin + (ic + ityp * 10 - 11) * 10, (
		    ftnlen)10);
	    do_fio(&c__1, mn7wrc_1.warmes + (ic + ityp * 10 - 11) * 60, (
		    ftnlen)60);
	    e_wsfe();
/* L300: */
	}
	mn7cnv_1.nwrmes[ityp - 1] = 0;
	io___1230.ciunit = mn7iou_1.isyswr;
	s_wsfe(&io___1230);
	e_wsfe();
    }
    return 0;
} /* mnwarn_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */


/* Subroutine */ int mnwerr_()
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k, l, k1;
    static doublereal ba, al, dx, du1, du2;
    static integer iin, ndex, ierr, ndiag;
    static doublereal denom;
    extern /* Subroutine */ int mnvert_(doublereal *, integer *, integer *, 
	    integer *, integer *);


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C          Calculates the WERR, external parameter errors, */
/* C      and the global correlation coefficients, to be called */
/* C      whenever a new covariance matrix is available. */
/* C */

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.2  2001/01/02 08:35:54  andras */
/* *** empty log message *** */

/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506cm.inc */


/*                         calculate external error if v exists */
    if (mn7flg_1.isw[1] >= 1) {
	i__1 = mn7npr_1.npar;
	for (l = 1; l <= i__1; ++l) {
	    ndex = l * (l + 1) / 2;
	    dx = sqrt((d__1 = mn7var_1.vhmat[ndex - 1] * mn7min_1.up, abs(
		    d__1)));
	    i__ = mn7inx_1.nexofi[l - 1];
	    if (mn7inx_1.nvarl[i__ - 1] > 1) {
		al = mn7ext_1.alim[i__ - 1];
		ba = mn7ext_1.blim[i__ - 1] - al;
		du1 = al + (sin(mn7int_1.x[l - 1] + dx) + (float)1.) * (float)
			.5 * ba - mn7ext_1.u[i__ - 1];
		du2 = al + (sin(mn7int_1.x[l - 1] - dx) + (float)1.) * (float)
			.5 * ba - mn7ext_1.u[i__ - 1];
		if (dx > (float)1.) {
		    du1 = ba;
		}
		dx = (abs(du1) + abs(du2)) * (float).5;
	    }
	    mn7err_1.werr[l - 1] = dx;
/* L100: */
	}
    }
/*                          global correlation coefficients */
    if (mn7flg_1.isw[1] >= 1) {
	i__1 = mn7npr_1.npar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    mn7err_1.globcc[i__ - 1] = (float)0.;
	    k1 = i__ * (i__ - 1) / 2;
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		k = k1 + j;
		mn7sim_1.p[i__ + j * 100 - 101] = mn7var_1.vhmat[k - 1];
/* L130: */
		mn7sim_1.p[j + i__ * 100 - 101] = mn7sim_1.p[i__ + j * 100 - 
			101];
	    }
	}
	mnvert_(mn7sim_1.p, &mn7npr_1.maxint, &mn7npr_1.maxint, &
		mn7npr_1.npar, &ierr);
	if (ierr == 0) {
	    i__2 = mn7npr_1.npar;
	    for (iin = 1; iin <= i__2; ++iin) {
		ndiag = iin * (iin + 1) / 2;
		denom = mn7sim_1.p[iin + iin * 100 - 101] * mn7var_1.vhmat[
			ndiag - 1];
		if (denom <= 1. && denom >= 0.) {
		    mn7err_1.globcc[iin - 1] = (float)0.;
		} else {
		    mn7err_1.globcc[iin - 1] = sqrt((float)1. - (float)1. / 
			    denom);
		}
/* L150: */
	    }
	}
    }
    return 0;
} /* mnwerr_ */


/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:20  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni */
/* Minuit */


/* Subroutine */ int stand_()
{

/* $Id: minuit_routines.c,v 1.1 2003/05/20 21:07:38 pln Exp $ */

/* $Log: minuit_routines.c,v $
/* Revision 1.1  2003/05/20 21:07:38  pln
/* Minuit optimizer from CERN
/* */
/* Revision 1.1.1.1  2000/06/08 11:19:21  andras */
/* import of MINUIT from CERNlib 2000 */

/* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni */
/* Minuit */




/* d506dp.inc */

/* ************ DOUBLE PRECISION VERSION ************* */
/* C        optional user-supplied subroutine is called whenever the */
/* C        command "standard" appears. */
/* C */
    return 0;
} /* stand_ */

#ifdef __cplusplus
	}
#endif
