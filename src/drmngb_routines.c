/*  -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c/f2c.h"

/* Table of constant values */

static integer c__3 = 3;
static doublereal c_b34 = 0.;
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_false = FALSE_;
static doublereal c_b52 = 1.;
static integer c__2 = 2;
static integer c__6 = 6;
static integer c__5 = 5;
static integer c__4 = 4;
static doublereal c_b403 = -1.;
static doublereal c_b433 = .33333333333333331;

/* Subroutine */ int da7sst_(integer *iv, integer *liv, integer *lv, 
	doublereal *v)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, nfc;
    static doublereal gts, emax, xmax, rfac1, emaxs;
    static logical goodx;


/*  ***  ASSESS CANDIDATE STEP (***SOL VERSION 2.3)  *** */


/*  ***  PURPOSE  *** */

/*        THIS SUBROUTINE IS CALLED BY AN UNCONSTRAINED MINIMIZATION */
/*     ROUTINE TO ASSESS THE NEXT CANDIDATE STEP.  IT MAY RECOMMEND ONE */
/*     OF SEVERAL COURSES OF ACTION, SUCH AS ACCEPTING THE STEP, RECOM- */
/*     PUTING IT USING THE SAME OR A NEW QUADRATIC MODEL, OR HALTING DUE */
/*     TO CONVERGENCE OR FALSE CONVERGENCE.  SEE THE RETURN CODE LISTING */
/*     BELOW. */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/*  IV (I/O) INTEGER PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION */
/*             BELOW OF IV VALUES REFERENCED. */
/* LIV (IN)  LENGTH OF IV ARRAY. */
/*  LV (IN)  LENGTH OF V ARRAY. */
/*   V (I/O) REAL PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION */
/*             BELOW OF V VALUES REFERENCED. */

/*  ***  IV VALUES REFERENCED  *** */

/*    IV(IRC) (I/O) ON INPUT FOR THE FIRST STEP TRIED IN A NEW ITERATION, */
/*             IV(IRC) SHOULD BE SET TO 3 OR 4 (THE VALUE TO WHICH IT IS */
/*             SET WHEN STEP IS DEFINITELY TO BE ACCEPTED).  ON INPUT */
/*             AFTER STEP HAS BEEN RECOMPUTED, IV(IRC) SHOULD BE */
/*             UNCHANGED SINCE THE PREVIOUS RETURN OF DA7SST. */
/*                ON OUTPUT, IV(IRC) IS A RETURN CODE HAVING ONE OF THE */
/*             FOLLOWING VALUES... */
/*                  1 = SWITCH MODELS OR TRY SMALLER STEP. */
/*                  2 = SWITCH MODELS OR ACCEPT STEP. */
/*                  3 = ACCEPT STEP AND DETERMINE V(RADFAC) BY GRADIENT */
/*                       TESTS. */
/*                  4 = ACCEPT STEP, V(RADFAC) HAS BEEN DETERMINED. */
/*                  5 = RECOMPUTE STEP (USING THE SAME MODEL). */
/*                  6 = RECOMPUTE STEP WITH RADIUS = V(LMAXS) BUT DO NOT */
/*                       EVALUATE THE OBJECTIVE FUNCTION. */
/*                  7 = X-CONVERGENCE (SEE V(XCTOL)). */
/*                  8 = RELATIVE FUNCTION CONVERGENCE (SEE V(RFCTOL)). */
/*                  9 = BOTH X- AND RELATIVE FUNCTION CONVERGENCE. */
/*                 10 = ABSOLUTE FUNCTION CONVERGENCE (SEE V(AFCTOL)). */
/*                 11 = SINGULAR CONVERGENCE (SEE V(LMAXS)). */
/*                 12 = FALSE CONVERGENCE (SEE V(XFTOL)). */
/*                 13 = IV(IRC) WAS OUT OF RANGE ON INPUT. */
/*             RETURN CODE I HAS PRECEDENCE OVER I+1 FOR I = 9, 10, 11. */
/* IV(MLSTGD) (I/O) SAVED VALUE OF IV(MODEL). */
/*  IV(MODEL) (I/O) ON INPUT, IV(MODEL) SHOULD BE AN INTEGER IDENTIFYING */
/*             THE CURRENT QUADRATIC MODEL OF THE OBJECTIVE FUNCTION. */
/*             IF A PREVIOUS STEP YIELDED A BETTER FUNCTION REDUCTION, */
/*             THEN IV(MODEL) WILL BE SET TO IV(MLSTGD) ON OUTPUT. */
/* IV(NFCALL) (IN)  INVOCATION COUNT FOR THE OBJECTIVE FUNCTION. */
/* IV(NFGCAL) (I/O) VALUE OF IV(NFCALL) AT STEP THAT GAVE THE BIGGEST */
/*             FUNCTION REDUCTION THIS ITERATION.  IV(NFGCAL) REMAINS */
/*             UNCHANGED UNTIL A FUNCTION REDUCTION IS OBTAINED. */
/* IV(RADINC) (I/O) THE NUMBER OF RADIUS INCREASES (OR MINUS THE NUMBER */
/*             OF DECREASES) SO FAR THIS ITERATION. */
/* IV(RESTOR) (OUT) SET TO 1 IF V(F) HAS BEEN RESTORED AND X SHOULD BE */
/*             RESTORED TO ITS INITIAL VALUE, TO 2 IF X SHOULD BE SAVED, */
/*             TO 3 IF X SHOULD BE RESTORED FROM THE SAVED VALUE, AND TO */
/*             0 OTHERWISE. */
/*  IV(STAGE) (I/O) COUNT OF THE NUMBER OF MODELS TRIED SO FAR IN THE */
/*             CURRENT ITERATION. */
/* IV(STGLIM) (IN)  MAXIMUM NUMBER OF MODELS TO CONSIDER. */
/* IV(SWITCH) (OUT) SET TO 0 UNLESS A NEW MODEL IS BEING TRIED AND IT */
/*             GIVES A SMALLER FUNCTION VALUE THAN THE PREVIOUS MODEL, */
/*             IN WHICH CASE DA7SST SETS IV(SWITCH) = 1. */
/* IV(TOOBIG) (I/O)  IS NONZERO ON INPUT IF STEP WAS TOO BIG (E.G., IF */
/*             IT WOULD CAUSE OVERFLOW).  IT IS SET TO 0 ON RETURN. */
/*   IV(XIRC) (I/O) VALUE THAT IV(IRC) WOULD HAVE IN THE ABSENCE OF */
/*             CONVERGENCE, FALSE CONVERGENCE, AND OVERSIZED STEPS. */

/*  ***  V VALUES REFERENCED  *** */

/* V(AFCTOL) (IN)  ABSOLUTE FUNCTION CONVERGENCE TOLERANCE.  IF THE */
/*             ABSOLUTE VALUE OF THE CURRENT FUNCTION VALUE V(F) IS LESS */
/*             THAN V(AFCTOL) AND DA7SST DOES NOT RETURN WITH */
/*             IV(IRC) = 11, THEN DA7SST RETURNS WITH IV(IRC) = 10. */
/* V(DECFAC) (IN)  FACTOR BY WHICH TO DECREASE RADIUS WHEN IV(TOOBIG) IS */
/*             NONZERO. */
/* V(DSTNRM) (IN)  THE 2-NORM OF D*STEP. */
/* V(DSTSAV) (I/O) VALUE OF V(DSTNRM) ON SAVED STEP. */
/*   V(DST0) (IN)  THE 2-NORM OF D TIMES THE NEWTON STEP (WHEN DEFINED, */
/*             I.E., FOR V(NREDUC) .GE. 0). */
/*      V(F) (I/O) ON BOTH INPUT AND OUTPUT, V(F) IS THE OBJECTIVE FUNC- */
/*             TION VALUE AT X.  IF X IS RESTORED TO A PREVIOUS VALUE, */
/*             THEN V(F) IS RESTORED TO THE CORRESPONDING VALUE. */
/*   V(FDIF) (OUT) THE FUNCTION REDUCTION V(F0) - V(F) (FOR THE OUTPUT */
/*             VALUE OF V(F) IF AN EARLIER STEP GAVE A BIGGER FUNCTION */
/*             DECREASE, AND FOR THE INPUT VALUE OF V(F) OTHERWISE). */
/* V(FLSTGD) (I/O) SAVED VALUE OF V(F). */
/*     V(F0) (IN)  OBJECTIVE FUNCTION VALUE AT START OF ITERATION. */
/* V(GTSLST) (I/O) VALUE OF V(GTSTEP) ON SAVED STEP. */
/* V(GTSTEP) (IN)  INNER PRODUCT BETWEEN STEP AND GRADIENT. */
/* V(INCFAC) (IN)  MINIMUM FACTOR BY WHICH TO INCREASE RADIUS. */
/*  V(LMAXS) (IN)  MAXIMUM REASONABLE STEP SIZE (AND INITIAL STEP BOUND). */
/*             IF THE ACTUAL FUNCTION DECREASE IS NO MORE THAN TWICE */
/*             WHAT WAS PREDICTED, IF A RETURN WITH IV(IRC) = 7, 8, OR 9 */
/*             DOES NOT OCCUR, IF V(DSTNRM) .GT. V(LMAXS) OR THE CURRENT */
/*             STEP IS A NEWTON STEP, AND IF */
/*             V(PREDUC) .LE. V(SCTOL) * ABS(V(F0)), THEN DA7SST RETURNS */
/*             WITH IV(IRC) = 11.  IF SO DOING APPEARS WORTHWHILE, THEN */
/*            DA7SST REPEATS THIS TEST (DISALLOWING A FULL NEWTON STEP) */
/*             WITH V(PREDUC) COMPUTED FOR A STEP OF LENGTH V(LMAXS) */
/*             (BY A RETURN WITH IV(IRC) = 6). */
/* V(NREDUC) (I/O)  FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR */
/*             NEWTON STEP.  IF DA7SST IS CALLED WITH IV(IRC) = 6, I.E., */
/*             IF V(PREDUC) HAS BEEN COMPUTED WITH RADIUS = V(LMAXS) FOR */
/*             USE IN THE SINGULAR CONVERGENCE TEST, THEN V(NREDUC) IS */
/*             SET TO -V(PREDUC) BEFORE THE LATTER IS RESTORED. */
/* V(PLSTGD) (I/O) VALUE OF V(PREDUC) ON SAVED STEP. */
/* V(PREDUC) (I/O) FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR */
/*             CURRENT STEP. */
/* V(RADFAC) (OUT) FACTOR TO BE USED IN DETERMINING THE NEW RADIUS, */
/*             WHICH SHOULD BE V(RADFAC)*DST, WHERE  DST  IS EITHER THE */
/*             OUTPUT VALUE OF V(DSTNRM) OR THE 2-NORM OF */
/*             DIAG(NEWD)*STEP  FOR THE OUTPUT VALUE OF STEP AND THE */
/*             UPDATED VERSION, NEWD, OF THE SCALE VECTOR D.  FOR */
/*             IV(IRC) = 3, V(RADFAC) = 1.0 IS RETURNED. */
/* V(RDFCMN) (IN)  MINIMUM VALUE FOR V(RADFAC) IN TERMS OF THE INPUT */
/*             VALUE OF V(DSTNRM) -- SUGGESTED VALUE = 0.1. */
/* V(RDFCMX) (IN)  MAXIMUM VALUE FOR V(RADFAC) -- SUGGESTED VALUE = 4.0. */
/*  V(RELDX) (IN) SCALED RELATIVE CHANGE IN X CAUSED BY STEP, COMPUTED */
/*             (E.G.) BY FUNCTION  DRLDST  AS */
/*                 MAX (D(I)*ABS(X(I)-X0(I)), 1 .LE. I .LE. P) / */
/*                    MAX (D(I)*(ABS(X(I))+ABS(X0(I))), 1 .LE. I .LE. P). */
/* V(RFCTOL) (IN)  RELATIVE FUNCTION CONVERGENCE TOLERANCE.  IF THE */
/*             ACTUAL FUNCTION REDUCTION IS AT MOST TWICE WHAT WAS PRE- */
/*             DICTED AND  V(NREDUC) .LE. V(RFCTOL)*ABS(V(F0)),  THEN */
/*            DA7SST RETURNS WITH IV(IRC) = 8 OR 9. */
/*  V(SCTOL) (IN)  SINGULAR CONVERGENCE TOLERANCE -- SEE V(LMAXS). */
/* V(STPPAR) (IN)  MARQUARDT PARAMETER -- 0 MEANS FULL NEWTON STEP. */
/* V(TUNER1) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION */
/*             REDUCTION WAS MUCH LESS THAN EXPECTED.  SUGGESTED */
/*             VALUE = 0.1. */
/* V(TUNER2) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION */
/*             REDUCTION WAS LARGE ENOUGH TO ACCEPT STEP.  SUGGESTED */
/*             VALUE = 10**-4. */
/* V(TUNER3) (IN)  TUNING CONSTANT USED TO DECIDE IF THE RADIUS */
/*             SHOULD BE INCREASED.  SUGGESTED VALUE = 0.75. */
/*  V(XCTOL) (IN)  X-CONVERGENCE CRITERION.  IF STEP IS A NEWTON STEP */
/*             (V(STPPAR) = 0) HAVING V(RELDX) .LE. V(XCTOL) AND GIVING */
/*             AT MOST TWICE THE PREDICTED FUNCTION DECREASE, THEN */
/*            DA7SST RETURNS IV(IRC) = 7 OR 9. */
/*  V(XFTOL) (IN)  FALSE CONVERGENCE TOLERANCE.  IF STEP GAVE NO OR ONLY */
/*             A SMALL FUNCTION DECREASE AND V(RELDX) .LE. V(XFTOL), */
/*             THEN DA7SST RETURNS WITH IV(IRC) = 12. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  APPLICATION AND USAGE RESTRICTIONS  *** */

/*        THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR */
/*     LEAST-SQUARES) PACKAGE.  IT MAY BE USED IN ANY UNCONSTRAINED */
/*     MINIMIZATION SOLVER THAT USES DOGLEG, GOLDFELD-QUANDT-TROTTER, */
/*     OR LEVENBERG-MARQUARDT STEPS. */

/*  ***  ALGORITHM NOTES  *** */

/*        SEE (1) FOR FURTHER DISCUSSION OF THE ASSESSING AND MODEL */
/*     SWITCHING STRATEGIES.  WHILE NL2SOL CONSIDERS ONLY TWO MODELS, */
/*    DA7SST IS DESIGNED TO HANDLE ANY NUMBER OF MODELS. */

/*  ***  USAGE NOTES  *** */

/*        ON THE FIRST CALL OF AN ITERATION, ONLY THE I/O VARIABLES */
/*     STEP, X, IV(IRC), IV(MODEL), V(F), V(DSTNRM), V(GTSTEP), AND */
/*     V(PREDUC) NEED HAVE BEEN INITIALIZED.  BETWEEN CALLS, NO I/O */
/*     VALUES EXCEPT STEP, X, IV(MODEL), V(F) AND THE STOPPING TOLER- */
/*     ANCES SHOULD BE CHANGED. */
/*        AFTER A RETURN FOR CONVERGENCE OR FALSE CONVERGENCE, ONE CAN */
/*     CHANGE THE STOPPING TOLERANCES AND CALL DA7SST AGAIN, IN WHICH */
/*     CASE THE STOPPING TESTS WILL BE REPEATED. */

/*  ***  REFERENCES  *** */

/*     (1) DENNIS, J.E., JR., GAY, D.M., AND WELSCH, R.E. (1981), */
/*        AN ADAPTIVE NONLINEAR LEAST-SQUARES ALGORITHM, */
/*        ACM TRANS. MATH. SOFTWARE, VOL. 7, NO. 3. */

/*     (2) POWELL, M.J.D. (1970)  A FORTRAN SUBROUTINE FOR SOLVING */
/*        SYSTEMS OF NONLINEAR ALGEBRAIC EQUATIONS, IN NUMERICAL */
/*        METHODS FOR NONLINEAR ALGEBRAIC EQUATIONS, EDITED BY */
/*        P. RABINOWITZ, GORDON AND BREACH, LONDON. */

/*  ***  HISTORY  *** */

/*        JOHN DENNIS DESIGNED MUCH OF THIS ROUTINE, STARTING WITH */
/*     IDEAS IN (2). ROY WELSCH SUGGESTED THE MODEL SWITCHING STRATEGY. */
/*        DAVID GAY AND STEPHEN PETERS CAST THIS SUBROUTINE INTO A MORE */
/*     PORTABLE FORM (WINTER 1977), AND DAVID GAY CAST IT INTO ITS */
/*     PRESENT FORM (FALL 1978), WITH MINOR CHANGES TO THE SINGULAR */
/*     CONVERGENCE TEST IN MAY, 1984 (TO DEAL WITH FULL NEWTON STEPS). */

/*  ***  GENERAL  *** */

/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND */
/*     MCS-7906671. */

/* ------------------------  EXTERNAL QUANTITIES  ------------------------ */

/*  ***  NO EXTERNAL FUNCTIONS AND SUBROUTINES  *** */

/* --------------------------  LOCAL VARIABLES  -------------------------- */


/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  DATA INITIALIZATIONS  *** */

/* /6 */
/*     DATA HALF/0.5D+0/, ONE/1.D+0/, ONEP2/1.2D+0/, TWO/2.D+0/, */
/*    1     ZERO/0.D+0/ */
/* /7 */
/* / */

/* /6 */
/*     DATA IRC/29/, MLSTGD/32/, MODEL/5/, NFCALL/6/, NFGCAL/7/, */
/*    1     RADINC/8/, RESTOR/9/, STAGE/10/, STGLIM/11/, SWITCH/12/, */
/*    2     TOOBIG/2/, XIRC/13/ */
/* /7 */
/* / */
/* /6 */
/*     DATA AFCTOL/31/, DECFAC/22/, DSTNRM/2/, DST0/3/, DSTSAV/18/, */
/*    1     F/10/, FDIF/11/, FLSTGD/12/, F0/13/, GTSLST/14/, GTSTEP/4/, */
/*    2     INCFAC/23/, LMAXS/36/, NREDUC/6/, PLSTGD/15/, PREDUC/7/, */
/*    3     RADFAC/16/, RDFCMN/24/, RDFCMX/25/, RELDX/17/, RFCTOL/32/, */
/*    4     SCTOL/37/, STPPAR/5/, TUNER1/26/, TUNER2/27/, TUNER3/28/, */
/*    5     XCTOL/33/, XFTOL/34/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --v;

    /* Function Body */
    nfc = iv[6];
    iv[12] = 0;
    iv[9] = 0;
    rfac1 = 1.;
    goodx = TRUE_;
    i__ = iv[29];
    if (i__ >= 1 && i__ <= 12) {
	switch (i__) {
	    case 1:  goto L20;
	    case 2:  goto L30;
	    case 3:  goto L10;
	    case 4:  goto L10;
	    case 5:  goto L40;
	    case 6:  goto L280;
	    case 7:  goto L220;
	    case 8:  goto L220;
	    case 9:  goto L220;
	    case 10:  goto L220;
	    case 11:  goto L220;
	    case 12:  goto L170;
	}
    }
    iv[29] = 13;
    goto L999;

/*  ***  INITIALIZE FOR NEW ITERATION  *** */

L10:
    iv[10] = 1;
    iv[8] = 0;
    v[12] = v[13];
    if (iv[2] == 0) {
	goto L110;
    }
    iv[10] = -1;
    iv[13] = i__;
    goto L60;

/*  ***  STEP WAS RECOMPUTED WITH NEW MODEL OR SMALLER RADIUS  *** */
/*  ***  FIRST DECIDE WHICH  *** */

L20:
    if (iv[5] != iv[32]) {
	goto L30;
    }
/*        ***  OLD MODEL RETAINED, SMALLER RADIUS TRIED  *** */
/*        ***  DO NOT CONSIDER ANY MORE NEW MODELS THIS ITERATION  *** */
    iv[10] = iv[11];
    iv[8] = -1;
    goto L110;

/*  ***  A NEW MODEL IS BEING TRIED.  DECIDE WHETHER TO KEEP IT.  *** */

L30:
    ++iv[10];

/*     ***  NOW WE ADD THE POSSIBILITY THAT STEP WAS RECOMPUTED WITH  *** */
/*     ***  THE SAME MODEL, PERHAPS BECAUSE OF AN OVERSIZED STEP.     *** */

L40:
    if (iv[10] > 0) {
	goto L50;
    }

/*        ***  STEP WAS RECOMPUTED BECAUSE IT WAS TOO BIG.  *** */

    if (iv[2] != 0) {
	goto L60;
    }

/*        ***  RESTORE IV(STAGE) AND PICK UP WHERE WE LEFT OFF.  *** */

    iv[10] = -iv[10];
    i__ = iv[13];
    switch (i__) {
	case 1:  goto L20;
	case 2:  goto L30;
	case 3:  goto L110;
	case 4:  goto L110;
	case 5:  goto L70;
    }

L50:
    if (iv[2] == 0) {
	goto L70;
    }

/*  ***  HANDLE OVERSIZE STEP  *** */

    iv[2] = 0;
    if (iv[8] > 0) {
	goto L80;
    }
    iv[10] = -iv[10];
    iv[13] = iv[29];

L60:
    iv[2] = 0;
    v[16] = v[22];
    --iv[8];
    iv[29] = 5;
    iv[9] = 1;
    v[10] = v[12];
    goto L999;

L70:
    if (v[10] < v[12]) {
	goto L110;
    }

/*     *** THE NEW STEP IS A LOSER.  RESTORE OLD MODEL.  *** */

    if (iv[5] == iv[32]) {
	goto L80;
    }
    iv[5] = iv[32];
    iv[12] = 1;

/*     ***  RESTORE STEP, ETC. ONLY IF A PREVIOUS STEP DECREASED V(F). */

L80:
    if (v[12] >= v[13]) {
	goto L110;
    }
    if (iv[10] < iv[11]) {
	goodx = FALSE_;
    } else if (nfc < iv[7] + iv[11] + 2) {
	goodx = FALSE_;
    } else if (iv[12] != 0) {
	goodx = FALSE_;
    }
    iv[9] = 3;
    v[10] = v[12];
    v[7] = v[15];
    v[4] = v[14];
    if (iv[12] == 0) {
	rfac1 = v[2] / v[18];
    }
    v[2] = v[18];
    if (goodx) {

/*     ***  ACCEPT PREVIOUS SLIGHTLY REDUCING STEP *** */

	v[11] = v[13] - v[10];
	iv[29] = 4;
	v[16] = rfac1;
	goto L999;
    }
    nfc = iv[7];

L110:
    v[11] = v[13] - v[10];
    if (v[11] > v[27] * v[7]) {
	goto L140;
    }
    if (iv[8] > 0) {
	goto L140;
    }

/*        ***  NO (OR ONLY A TRIVIAL) FUNCTION DECREASE */
/*        ***  -- SO TRY NEW MODEL OR SMALLER RADIUS */

    if (v[10] < v[13]) {
	goto L120;
    }
    iv[32] = iv[5];
    v[12] = v[10];
    v[10] = v[13];
    iv[9] = 1;
    goto L130;
L120:
    iv[7] = nfc;
L130:
    iv[29] = 1;
    if (iv[10] < iv[11]) {
	goto L160;
    }
    iv[29] = 5;
    --iv[8];
    goto L160;

/*  ***  NONTRIVIAL FUNCTION DECREASE ACHIEVED  *** */

L140:
    iv[7] = nfc;
    rfac1 = 1.;
    v[18] = v[2];
    if (v[11] > v[7] * v[26]) {
	goto L190;
    }

/*  ***  DECREASE WAS MUCH LESS THAN PREDICTED -- EITHER CHANGE MODELS */
/*  ***  OR ACCEPT STEP WITH DECREASED RADIUS. */

    if (iv[10] >= iv[11]) {
	goto L150;
    }
/*        ***  CONSIDER SWITCHING MODELS  *** */
    iv[29] = 2;
    goto L160;

/*     ***  ACCEPT STEP WITH DECREASED RADIUS  *** */

L150:
    iv[29] = 4;

/*  ***  SET V(RADFAC) TO FLETCHER*S DECREASE FACTOR  *** */

L160:
    iv[13] = iv[29];
    emax = v[4] + v[11];
    v[16] = rfac1 * .5;
    if (emax < v[4]) {
/* Computing MAX */
	d__1 = v[24], d__2 = v[4] * .5 / emax;
	v[16] = rfac1 * max(d__1,d__2);
    }

/*  ***  DO FALSE CONVERGENCE TEST  *** */

L170:
    if (v[17] <= v[34]) {
	goto L180;
    }
    iv[29] = iv[13];
    if (v[10] < v[13]) {
	goto L200;
    }
    goto L230;

L180:
    iv[29] = 12;
    goto L240;

/*  ***  HANDLE GOOD FUNCTION DECREASE  *** */

L190:
    if (v[11] < -v[28] * v[4]) {
	goto L210;
    }

/*     ***  INCREASING RADIUS LOOKS WORTHWHILE.  SEE IF WE JUST */
/*     ***  RECOMPUTED STEP WITH A DECREASED RADIUS OR RESTORED STEP */
/*     ***  AFTER RECOMPUTING IT WITH A LARGER RADIUS. */

    if (iv[8] < 0) {
	goto L210;
    }
    if (iv[9] == 1) {
	goto L210;
    }
    if (iv[9] == 3) {
	goto L210;
    }

/*        ***  WE DID NOT.  TRY A LONGER STEP UNLESS THIS WAS A NEWTON */
/*        ***  STEP. */

    v[16] = v[25];
    gts = v[4];
    if (v[11] < (.5 / v[16] - 1.) * gts) {
/* Computing MAX */
	d__1 = v[23], d__2 = gts * .5 / (gts + v[11]);
	v[16] = max(d__1,d__2);
    }
    iv[29] = 4;
    if (v[5] == 0.) {
	goto L230;
    }
    if (v[3] >= 0. && (v[3] < v[2] * 2. || v[6] < v[11] * 1.2)) {
	goto L230;
    }
/*             ***  STEP WAS NOT A NEWTON STEP.  RECOMPUTE IT WITH */
/*             ***  A LARGER RADIUS. */
    iv[29] = 5;
    ++iv[8];

/*  ***  SAVE VALUES CORRESPONDING TO GOOD STEP  *** */

L200:
    v[12] = v[10];
    iv[32] = iv[5];
    if (iv[9] == 0) {
	iv[9] = 2;
    }
    v[18] = v[2];
    iv[7] = nfc;
    v[15] = v[7];
    v[14] = v[4];
    goto L230;

/*  ***  ACCEPT STEP WITH RADIUS UNCHANGED  *** */

L210:
    v[16] = 1.;
    iv[29] = 3;
    goto L230;

/*  ***  COME HERE FOR A RESTART AFTER CONVERGENCE  *** */

L220:
    iv[29] = iv[13];
    if (v[18] >= 0.) {
	goto L240;
    }
    iv[29] = 12;
    goto L240;

/*  ***  PERFORM CONVERGENCE TESTS  *** */

L230:
    iv[13] = iv[29];
L240:
    if (iv[9] == 1 && v[12] < v[13]) {
	iv[9] = 3;
    }
    if (abs(v[10]) < v[31]) {
	iv[29] = 10;
    }
    if (v[11] * .5 > v[7]) {
	goto L999;
    }
    emax = v[32] * abs(v[13]);
    emaxs = v[37] * abs(v[13]);
    if (v[7] <= emaxs && (v[2] > v[36] || v[5] == 0.)) {
	iv[29] = 11;
    }
    if (v[3] < 0.) {
	goto L250;
    }
    i__ = 0;
    if (v[6] > 0. && v[6] <= emax || v[6] == 0. && v[7] == 0.) {
	i__ = 2;
    }
    if (v[5] == 0. && v[17] <= v[33] && goodx) {
	++i__;
    }
    if (i__ > 0) {
	iv[29] = i__ + 6;
    }

/*  ***  CONSIDER RECOMPUTING STEP OF LENGTH V(LMAXS) FOR SINGULAR */
/*  ***  CONVERGENCE TEST. */

L250:
    if (iv[29] > 5 && iv[29] != 12) {
	goto L999;
    }
    if (v[5] == 0.) {
	goto L999;
    }
    if (v[2] > v[36]) {
	goto L260;
    }
    if (v[7] >= emaxs) {
	goto L999;
    }
    if (v[3] <= 0.) {
	goto L270;
    }
    if (v[3] * .5 <= v[36]) {
	goto L999;
    }
    goto L270;
L260:
    if (v[2] * .5 <= v[36]) {
	goto L999;
    }
    xmax = v[36] / v[2];
    if (xmax * (2. - xmax) * v[7] >= emaxs) {
	goto L999;
    }
L270:
    if (v[6] < 0.) {
	goto L290;
    }

/*  ***  RECOMPUTE V(PREDUC) FOR USE IN SINGULAR CONVERGENCE TEST  *** */

    v[14] = v[4];
    v[18] = v[2];
    if (iv[29] == 12) {
	v[18] = -v[18];
    }
    v[15] = v[7];
    i__ = iv[9];
    iv[9] = 2;
    if (i__ == 3) {
	iv[9] = 0;
    }
    iv[29] = 6;
    goto L999;

/*  ***  PERFORM SINGULAR CONVERGENCE TEST WITH RECOMPUTED V(PREDUC)  *** */

L280:
    v[4] = v[14];
    v[2] = abs(v[18]);
    iv[29] = iv[13];
    if (v[18] <= 0.) {
	iv[29] = 12;
    }
    v[6] = -v[7];
    v[7] = v[15];
    iv[9] = 3;
L290:
    if (-v[6] <= v[37] * abs(v[13])) {
	iv[29] = 11;
    }

L999:
    return 0;

/*  ***  LAST LINE OF DA7SST FOLLOWS  *** */
} /* da7sst_ */

/* Subroutine */ int dd7dgb_(doublereal *b, doublereal *d__, doublereal *dig, 
	doublereal *dst, doublereal *g, integer *ipiv, integer *ka, 
	doublereal *l, integer *lv, integer *p, integer *pc, doublereal *
	nwtst, doublereal *step, doublereal *td, doublereal *tg, doublereal *
	v, doublereal *w, doublereal *x0)
{
    /* Initialized data */

    static doublereal meps2 = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer p1;
    static doublereal t1, t2, ti, xi, x0i, rad;
    static integer p1m1;
    static doublereal nred, pred, gnorm;
    extern /* Subroutine */ int dd7dog_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dr7mdc_(integer *);
    extern /* Subroutine */ int dv7shf_(integer *, integer *, doublereal *), 
	    dl7ivm_(integer *, doublereal *, doublereal *, doublereal *);
    static doublereal gnorm0;
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int i7shft_(integer *, integer *, integer *), 
	    dl7vml_(integer *, doublereal *, doublereal *, doublereal *), 
	    dv7scp_(integer *, doublereal *, doublereal *);
    extern doublereal dv2nrm_(integer *, doublereal *);
    extern /* Subroutine */ int dl7itv_(integer *, doublereal *, doublereal *,
	     doublereal *), dq7rsh_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *), dv7ipr_(integer *, 
	    integer *, doublereal *), dv7cpy_(integer *, doublereal *, 
	    doublereal *), dl7tvm_(integer *, doublereal *, doublereal *, 
	    doublereal *), dv2axy_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dv7vmp_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal ghinvg, dnwtst;


/*  ***  COMPUTE DOUBLE-DOGLEG STEP, SUBJECT TO SIMPLE BOUNDS ON X  *** */


/*     DIMENSION L(P*(P+1)/2) */


/*  ***  LOCAL VARIABLES  *** */


/*  ***  V SUBSCRIPTS  *** */


/* /6 */
/*     DATA DGNORM/1/, DST0/3/, DSTNRM/2/, GRDFAC/45/, GTHG/44/, */
/*    1     GTSTEP/4/, NREDUC/6/, NWTFAC/46/, PREDUC/7/, RADIUS/8/, */
/*    2     STPPAR/5/ */
/* /7 */
/* / */
/* /6 */
/*     DATA HALF/0.5D+0/, ONE/1.D+0/, TWO/2.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --l;
    --v;
    --x0;
    --w;
    --tg;
    --td;
    --step;
    --nwtst;
    --ipiv;
    --g;
    --dst;
    --dig;
    --d__;
    b -= 3;

    /* Function Body */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    if (meps2 <= 0.) {
	meps2 = 2. * dr7mdc_(&c__3);
    }
    gnorm0 = v[1];
    v[2] = 0.;
    if (*ka < 0) {
	goto L10;
    }
    dnwtst = v[3];
    nred = v[6];
L10:
    pred = 0.;
    v[5] = 0.;
    rad = v[8];
    if (*pc > 0) {
	goto L20;
    }
    dnwtst = 0.;
    dv7scp_(p, &step[1], &c_b34);
    goto L140;

L20:
    p1 = *pc;
    dv7cpy_(p, &td[1], &d__[1]);
    dv7ipr_(p, &ipiv[1], &td[1]);
    dv7scp_(pc, &dst[1], &c_b34);
    dv7cpy_(p, &tg[1], &g[1]);
    dv7ipr_(p, &ipiv[1], &tg[1]);

L30:
    dl7ivm_(&p1, &nwtst[1], &l[1], &tg[1]);
    ghinvg = dd7tpr_(&p1, &nwtst[1], &nwtst[1]);
    v[6] = ghinvg * .5;
    dl7itv_(&p1, &nwtst[1], &l[1], &nwtst[1]);
    dv7vmp_(&p1, &step[1], &nwtst[1], &td[1], &c__1);
    v[3] = dv2nrm_(pc, &step[1]);
    if (*ka >= 0) {
	goto L40;
    }
    *ka = 0;
    dnwtst = v[3];
    nred = v[6];
L40:
    v[8] = rad - v[2];
    if (v[8] <= 0.) {
	goto L100;
    }
    dv7vmp_(&p1, &dig[1], &tg[1], &td[1], &c_n1);
    gnorm = dv2nrm_(&p1, &dig[1]);
    if (gnorm <= 0.) {
	goto L100;
    }
    v[1] = gnorm;
    dv7vmp_(&p1, &dig[1], &dig[1], &td[1], &c_n1);
    dl7tvm_(&p1, &w[1], &l[1], &dig[1]);
    v[44] = dv2nrm_(&p1, &w[1]);
    ++(*ka);
    dd7dog_(&dig[1], lv, &p1, &nwtst[1], &step[1], &v[1]);

/*     ***  FIND T SUCH THAT X - T*STEP IS STILL FEASIBLE. */

    t = 1.;
    k = 0;
    i__1 = p1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ipiv[i__];
	x0i = x0[j] + dst[i__] / td[i__];
	xi = x0i + step[i__];
	if (xi < b[(j << 1) + 1]) {
	    goto L50;
	}
	if (xi <= b[(j << 1) + 2]) {
	    goto L70;
	}
	ti = (b[(j << 1) + 2] - x0i) / step[i__];
	j = i__;
	goto L60;
L50:
	ti = (b[(j << 1) + 1] - x0i) / step[i__];
	j = -i__;
L60:
	if (t <= ti) {
	    goto L70;
	}
	k = j;
	t = ti;
L70:
	;
    }

/*  ***  UPDATE DST, TG, AND PRED  *** */

    dv7vmp_(&p1, &step[1], &step[1], &td[1], &c__1);
    dv2axy_(&p1, &dst[1], &t, &step[1], &dst[1]);
    v[2] = dv2nrm_(pc, &dst[1]);
    t1 = t * v[45];
    t2 = t * v[46];
/* Computing 2nd power */
    d__1 = v[44] * t1;
    pred = pred - t1 * gnorm * ((t2 + 1.) * gnorm) - t2 * (t2 * .5 + 1.) * 
	    ghinvg - d__1 * d__1 * .5;
    if (k == 0) {
	goto L100;
    }
    dl7vml_(&p1, &w[1], &l[1], &w[1]);
    t2 = 1. - t2;
    i__1 = p1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	tg[i__] = t2 * tg[i__] - t1 * w[i__];
    }

/*     ***  PERMUTE L, ETC. IF NECESSARY  *** */

    p1m1 = p1 - 1;
    j = abs(k);
    if (j == p1) {
	goto L90;
    }
    dq7rsh_(&j, &p1, &c_false, &tg[1], &l[1], &w[1]);
    i7shft_(&p1, &j, &ipiv[1]);
    dv7shf_(&p1, &j, &tg[1]);
    dv7shf_(&p1, &j, &td[1]);
    dv7shf_(&p1, &j, &dst[1]);
L90:
    if (k < 0) {
	ipiv[p1] = -ipiv[p1];
    }
    p1 = p1m1;
    if (p1 > 0) {
	goto L30;
    }

/*     ***  UNSCALE STEP, UPDATE X AND DIHDI  *** */

L100:
    dv7scp_(p, &step[1], &c_b34);
    i__1 = *pc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = (i__2 = ipiv[i__], abs(i__2));
	step[j] = dst[i__] / td[i__];
/* L110: */
    }

/*  ***  FUDGE STEP TO ENSURE THAT IT FORCES APPROPRIATE COMPONENTS */
/*  ***  TO THEIR BOUNDS  *** */

    if (p1 >= *pc) {
	goto L140;
    }
    dv2axy_(p, &td[1], &c_b52, &step[1], &x0[1]);
    k = p1 + 1;
    i__1 = *pc;
    for (i__ = k; i__ <= i__1; ++i__) {
	j = ipiv[i__];
	t = meps2;
	if (j > 0) {
	    goto L120;
	}
	t = -t;
	j = -j;
	ipiv[i__] = j;
L120:
/* Computing MAX */
	d__3 = (d__1 = td[j], abs(d__1)), d__4 = (d__2 = x0[j], abs(d__2));
	t *= max(d__3,d__4);
	step[j] += t;
/* L130: */
    }

L140:
    v[1] = gnorm0;
    v[6] = nred;
    v[7] = pred;
    v[8] = rad;
    v[3] = dnwtst;
    v[4] = dd7tpr_(p, &step[1], &g[1]);

/* L999: */
    return 0;
/*  ***  LAST LINE OF DD7DGB FOLLOWS  *** */
} /* dd7dgb_ */

/* Subroutine */ int dd7dog_(doublereal *dig, integer *lv, integer *n, 
	doublereal *nwtstp, doublereal *step, doublereal *v)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t, t1, t2, cfact, relax, cnorm, gnorm, rlambd, ghinvg, 
	    femnsq, ctrnwt, nwtnrm;


/*  ***  COMPUTE DOUBLE DOGLEG STEP  *** */

/*  ***  PARAMETER DECLARATIONS  *** */


/*  ***  PURPOSE  *** */

/*        THIS SUBROUTINE COMPUTES A CANDIDATE STEP (FOR USE IN AN UNCON- */
/*     STRAINED MINIMIZATION CODE) BY THE DOUBLE DOGLEG ALGORITHM OF */
/*     DENNIS AND MEI (REF. 1), WHICH IS A VARIATION ON POWELL*S DOGLEG */
/*     SCHEME (REF. 2, P. 95). */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/*    DIG (INPUT) DIAG(D)**-2 * G -- SEE ALGORITHM NOTES. */
/*      G (INPUT) THE CURRENT GRADIENT VECTOR. */
/*     LV (INPUT) LENGTH OF V. */
/*      N (INPUT) NUMBER OF COMPONENTS IN  DIG, G, NWTSTP,  AND  STEP. */
/* NWTSTP (INPUT) NEGATIVE NEWTON STEP -- SEE ALGORITHM NOTES. */
/*   STEP (OUTPUT) THE COMPUTED STEP. */
/*      V (I/O) VALUES ARRAY, THE FOLLOWING COMPONENTS OF WHICH ARE */
/*             USED HERE... */
/* V(BIAS)   (INPUT) BIAS FOR RELAXED NEWTON STEP, WHICH IS V(BIAS) OF */
/*             THE WAY FROM THE FULL NEWTON TO THE FULLY RELAXED NEWTON */
/*             STEP.  RECOMMENDED VALUE = 0.8 . */
/* V(DGNORM) (INPUT) 2-NORM OF DIAG(D)**-1 * G -- SEE ALGORITHM NOTES. */
/* V(DSTNRM) (OUTPUT) 2-NORM OF DIAG(D) * STEP, WHICH IS V(RADIUS) */
/*             UNLESS V(STPPAR) = 0 -- SEE ALGORITHM NOTES. */
/* V(DST0) (INPUT) 2-NORM OF DIAG(D) * NWTSTP -- SEE ALGORITHM NOTES. */
/* V(GRDFAC) (OUTPUT) THE COEFFICIENT OF  DIG  IN THE STEP RETURNED -- */
/*             STEP(I) = V(GRDFAC)*DIG(I) + V(NWTFAC)*NWTSTP(I). */
/* V(GTHG)   (INPUT) SQUARE-ROOT OF (DIG**T) * (HESSIAN) * DIG -- SEE */
/*             ALGORITHM NOTES. */
/* V(GTSTEP) (OUTPUT) INNER PRODUCT BETWEEN G AND STEP. */
/* V(NREDUC) (OUTPUT) FUNCTION REDUCTION PREDICTED FOR THE FULL NEWTON */
/*             STEP. */
/* V(NWTFAC) (OUTPUT) THE COEFFICIENT OF  NWTSTP  IN THE STEP RETURNED -- */
/*             SEE V(GRDFAC) ABOVE. */
/* V(PREDUC) (OUTPUT) FUNCTION REDUCTION PREDICTED FOR THE STEP RETURNED. */
/* V(RADIUS) (INPUT) THE TRUST REGION RADIUS.  D TIMES THE STEP RETURNED */
/*             HAS 2-NORM V(RADIUS) UNLESS V(STPPAR) = 0. */
/* V(STPPAR) (OUTPUT) CODE TELLING HOW STEP WAS COMPUTED... 0 MEANS A */
/*             FULL NEWTON STEP.  BETWEEN 0 AND 1 MEANS V(STPPAR) OF THE */
/*             WAY FROM THE NEWTON TO THE RELAXED NEWTON STEP.  BETWEEN */
/*             1 AND 2 MEANS A TRUE DOUBLE DOGLEG STEP, V(STPPAR) - 1 OF */
/*             THE WAY FROM THE RELAXED NEWTON TO THE CAUCHY STEP. */
/*             GREATER THAN 2 MEANS 1 / (V(STPPAR) - 1) TIMES THE CAUCHY */
/*             STEP. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  ALGORITHM NOTES  *** */

/*        LET  G  AND  H  BE THE CURRENT GRADIENT AND HESSIAN APPROXIMA- */
/*     TION RESPECTIVELY AND LET D BE THE CURRENT SCALE VECTOR.  THIS */
/*     ROUTINE ASSUMES DIG = DIAG(D)**-2 * G  AND  NWTSTP = H**-1 * G. */
/*     THE STEP COMPUTED IS THE SAME ONE WOULD GET BY REPLACING G AND H */
/*     BY  DIAG(D)**-1 * G  AND  DIAG(D)**-1 * H * DIAG(D)**-1, */
/*     COMPUTING STEP, AND TRANSLATING STEP BACK TO THE ORIGINAL */
/*     VARIABLES, I.E., PREMULTIPLYING IT BY DIAG(D)**-1. */

/*  ***  REFERENCES  *** */

/* 1.  DENNIS, J.E., AND MEI, H.H.W. (1979), TWO NEW UNCONSTRAINED OPTI- */
/*             MIZATION ALGORITHMS WHICH USE FUNCTION AND GRADIENT */
/*             VALUES, J. OPTIM. THEORY APPLIC. 28, PP. 453-482. */
/* 2. POWELL, M.J.D. (1970), A HYBRID METHOD FOR NON-LINEAR EQUATIONS, */
/*             IN NUMERICAL METHODS FOR NON-LINEAR EQUATIONS, EDITED BY */
/*             P. RABINOWITZ, GORDON AND BREACH, LONDON. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED */
/*     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND */
/*     MCS-7906671. */

/* ------------------------  EXTERNAL QUANTITIES  ------------------------ */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/* --------------------------  LOCAL VARIABLES  -------------------------- */


/*  ***  V SUBSCRIPTS  *** */


/*  ***  DATA INITIALIZATIONS  *** */

/* /6 */
/*     DATA HALF/0.5D+0/, ONE/1.D+0/, TWO/2.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* /6 */
/*     DATA BIAS/43/, DGNORM/1/, DSTNRM/2/, DST0/3/, GRDFAC/45/, */
/*    1     GTHG/44/, GTSTEP/4/, NREDUC/6/, NWTFAC/46/, PREDUC/7/, */
/*    2     RADIUS/8/, STPPAR/5/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --v;
    --step;
    --nwtstp;
    --dig;

    /* Function Body */
    nwtnrm = v[3];
    rlambd = 1.;
    if (nwtnrm > 0.) {
	rlambd = v[8] / nwtnrm;
    }
    gnorm = v[1];
    ghinvg = v[6] * 2.;
    v[45] = 0.;
    v[46] = 0.;
    if (rlambd < 1.) {
	goto L30;
    }

/*        ***  THE NEWTON STEP IS INSIDE THE TRUST REGION  *** */

    v[5] = 0.;
    v[2] = nwtnrm;
    v[4] = -ghinvg;
    v[7] = v[6];
    v[46] = -1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	step[i__] = -nwtstp[i__];
    }
    goto L999;

L30:
    v[2] = v[8];
/* Computing 2nd power */
    d__1 = gnorm / v[44];
    cfact = d__1 * d__1;
/*     ***  CAUCHY STEP = -CFACT * G. */
    cnorm = gnorm * cfact;
    relax = 1. - v[43] * (1. - gnorm * cnorm / ghinvg);
    if (rlambd < relax) {
	goto L50;
    }

/*        ***  STEP IS BETWEEN RELAXED NEWTON AND FULL NEWTON STEPS  *** */

    v[5] = 1. - (rlambd - relax) / (1. - relax);
    t = -rlambd;
    v[4] = t * ghinvg;
    v[7] = rlambd * (1. - rlambd * .5) * ghinvg;
    v[46] = t;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	step[i__] = t * nwtstp[i__];
    }
    goto L999;

L50:
    if (cnorm < v[8]) {
	goto L70;
    }

/*        ***  THE CAUCHY STEP LIES OUTSIDE THE TRUST REGION -- */
/*        ***  STEP = SCALED CAUCHY STEP  *** */

    t = -v[8] / gnorm;
    v[45] = t;
    v[5] = cnorm / v[8] + 1.;
    v[4] = -v[8] * gnorm;
/* Computing 2nd power */
    d__1 = v[44] / gnorm;
    v[7] = v[8] * (gnorm - v[8] * .5 * (d__1 * d__1));
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	step[i__] = t * dig[i__];
    }
    goto L999;

/*     ***  COMPUTE DOGLEG STEP BETWEEN CAUCHY AND RELAXED NEWTON  *** */
/*     ***  FEMUR = RELAXED NEWTON STEP MINUS CAUCHY STEP  *** */

L70:
    ctrnwt = cfact * relax * ghinvg / gnorm;
/*     *** CTRNWT = INNER PROD. OF CAUCHY AND RELAXED NEWTON STEPS, */
/*     *** SCALED BY GNORM**-1. */
/* Computing 2nd power */
    d__1 = cfact;
    t1 = ctrnwt - gnorm * (d__1 * d__1);
/*     ***  T1 = INNER PROD. OF FEMUR AND CAUCHY STEP, SCALED BY */
/*     ***  GNORM**-1. */
/* Computing 2nd power */
    d__1 = cfact;
    t2 = v[8] * (v[8] / gnorm) - gnorm * (d__1 * d__1);
    t = relax * nwtnrm;
    femnsq = t / gnorm * t - ctrnwt - t1;
/*     ***  FEMNSQ = SQUARE OF 2-NORM OF FEMUR, SCALED BY GNORM**-1. */
/* Computing 2nd power */
    d__1 = t1;
    t = t2 / (t1 + sqrt(d__1 * d__1 + femnsq * t2));
/*     ***  DOGLEG STEP  =  CAUCHY STEP  +  T * FEMUR. */
    t1 = (t - 1.) * cfact;
    v[45] = t1;
    t2 = -t * relax;
    v[46] = t2;
    v[5] = 2. - t;
/* Computing 2nd power */
    d__1 = gnorm;
    v[4] = t1 * (d__1 * d__1) + t2 * ghinvg;
/* Computing 2nd power */
    d__1 = v[44] * t1;
    v[7] = -t1 * gnorm * ((t2 + 1.) * gnorm) - t2 * (t2 * .5 + 1.) * ghinvg - 
	    d__1 * d__1 * .5;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	step[i__] = t1 * dig[i__] + t2 * nwtstp[i__];
    }

L999:
    return 0;
/*  ***  LAST LINE OF DD7DOG FOLLOWS  *** */
} /* dd7dog_ */

doublereal dd7tpr_(integer *p, doublereal *x, doublereal *y)
{
    /* Initialized data */

    static doublereal sqteta = 0.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__;
    static doublereal t;
    extern doublereal dr7mdc_(integer *);


/*  ***  RETURN THE INNER PRODUCT OF THE P-VECTORS X AND Y.  *** */



/*  ***  DR7MDC(2) RETURNS A MACHINE-DEPENDENT CONSTANT, SQTETA, WHICH */
/*  ***  IS SLIGHTLY LARGER THAN THE SMALLEST POSITIVE NUMBER THAT */
/*  ***  CAN BE SQUARED WITHOUT UNDERFLOWING. */

/* /6 */
/*     DATA ONE/1.D+0/, SQTETA/0.D+0/, ZERO/0.D+0/ */
/* /7 */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
/* / */

    ret_val = 0.;
    if (*p <= 0) {
	goto L999;
    }
    if (sqteta == 0.) {
	sqteta = dr7mdc_(&c__2);
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__3 = (d__1 = x[i__], abs(d__1)), d__4 = (d__2 = y[i__], abs(d__2));
	t = max(d__3,d__4);
	if (t > 1.) {
	    goto L10;
	}
	if (t < sqteta) {
	    goto L20;
	}
	t = x[i__] / sqteta * y[i__];
	if (abs(t) < sqteta) {
	    goto L20;
	}
L10:
	ret_val += x[i__] * y[i__];
L20:
	;
    }

L999:
    return ret_val;
/*  ***  LAST LINE OF DD7TPR FOLLOWS  *** */
} /* dd7tpr_ */

/* Subroutine */ int dh2rfa_(integer *n, doublereal *a, doublereal *b, 
	doublereal *x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal t;


/*  ***  APPLY 2X2 HOUSEHOLDER REFLECTION DETERMINED BY X, Y, Z TO */
/*  ***  N-VECTORS A, B  *** */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = a[i__] * *x + b[i__] * *y;
	a[i__] += t;
	b[i__] += t * *z__;
/* L10: */
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DH2RFA FOLLOWS  *** */
} /* dh2rfa_ */

doublereal dh2rfg_(doublereal *a, doublereal *b, doublereal *x, doublereal *y,
	 doublereal *z__)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__, t, a1, b1;


/*  ***  DETERMINE X, Y, Z SO  I + (1,Z)**T * (X,Y)  IS A 2X2 */
/*  ***  HOUSEHOLDER REFLECTION SENDING (A,B)**T INTO (C,0)**T, */
/*  ***  WHERE  C = -SIGN(A)*SQRT(A**2 + B**2)  IS THE VALUE DH2RFG */
/*  ***  RETURNS. */


/* /+ */
/* / */

/*  ***  BODY  *** */

    if (*b != zero) {
	goto L10;
    }
    *x = zero;
    *y = zero;
    *z__ = zero;
    ret_val = *a;
    goto L999;
L10:
    t = abs(*a) + abs(*b);
    a1 = *a / t;
    b1 = *b / t;
/* Computing 2nd power */
    d__1 = a1;
/* Computing 2nd power */
    d__2 = b1;
    c__ = sqrt(d__1 * d__1 + d__2 * d__2);
    if (a1 > zero) {
	c__ = -c__;
    }
    a1 -= c__;
    *z__ = b1 / a1;
    *x = a1 / c__;
    *y = b1 / c__;
    ret_val = t * c__;
L999:
    return ret_val;
/*  ***  LAST LINE OF DH2RFG FOLLOWS  *** */
} /* dh2rfg_ */

/* Subroutine */ int ditsum_(doublereal *d__, doublereal *g, integer *iv, 
	integer *liv, integer *lv, integer *p, doublereal *v, doublereal *x)
{
    /* Initialized data */

    static char model1[4*6+1] = "                  G   S ";
    static char model2[4*6+1] = " G   S  G-S S-G -S-G-G-S";

    /* Format strings */
    static char fmt_30[] = "(/\002   IT   NF\002,6x,\002F\002,7x,\002RELD\
F\002,3x,\002PRELDF\002,3x,\002RELDX\002,2x,\002MODEL  STPPAR\002)";
    static char fmt_40[] = "(/\002    IT   NF\002,7x,\002F\002,8x,\002RELD\
F\002,4x,\002PRELDF\002,4x,\002RELDX\002,3x,\002STPPAR\002)";
    static char fmt_100[] = "(i6,i5,d10.3,2d9.2,d8.1,a3,a4,2d8.1,d9.2)";
    static char fmt_110[] = "(i6,i5,d11.3,2d10.2,3d9.1,d10.2)";
    static char fmt_70[] = "(/\002    IT   NF\002,6x,\002F\002,7x,\002RELD\
F\002,3x,\002PRELDF\002,3x,\002RELDX\002,2x,\002MODEL  STPPAR\002,2x,\002D*S\
TEP\002,2x,\002NPRELDF\002)";
    static char fmt_80[] = "(/\002    IT   NF\002,7x,\002F\002,8x,\002RELD\
F\002,4x,\002PRELDF\002,4x,\002RELDX\002,3x,\002STPPAR\002,3x,\002D*STEP\002\
,3x,\002NPRELDF\002)";
    static char fmt_140[] = "(/\002 ***** X-CONVERGENCE *****\002)";
    static char fmt_160[] = "(/\002 ***** RELATIVE FUNCTION CONVERGENCE **\
***\002)";
    static char fmt_180[] = "(/\002 ***** X- AND RELATIVE FUNCTION CONVERGEN\
CE *****\002)";
    static char fmt_200[] = "(/\002 ***** ABSOLUTE FUNCTION CONVERGENCE **\
***\002)";
    static char fmt_220[] = "(/\002 ***** SINGULAR CONVERGENCE *****\002)";
    static char fmt_240[] = "(/\002 ***** FALSE CONVERGENCE *****\002)";
    static char fmt_260[] = "(/\002 ***** FUNCTION EVALUATION LIMIT *****\
\002)";
    static char fmt_280[] = "(/\002 ***** ITERATION LIMIT *****\002)";
    static char fmt_300[] = "(/\002 ***** STOPX *****\002)";
    static char fmt_320[] = "(/\002 ***** INITIAL F(X) CANNOT BE COMPUTED **\
***\002)";
    static char fmt_340[] = "(/\002 ***** BAD PARAMETERS TO ASSESS *****\002)"
	    ;
    static char fmt_360[] = "(/\002 ***** GRADIENT COULD NOT BE COMPUTED ***\
**\002)";
    static char fmt_380[] = "(/\002 ***** IV(1) =\002,i5,\002 *****\002)";
    static char fmt_400[] = "(/\002     I     INITIAL X(I)\002,8x,\002D(I\
)\002//(1x,i5,d17.6,d14.3))";
    static char fmt_410[] = "(/\002     0\002,i5,d10.3)";
    static char fmt_420[] = "(/\002     0\002,i5,d11.3)";
    static char fmt_450[] = "(/\002 FUNCTION\002,d17.6,\002   RELDX\002,d17.\
3/\002 FUNC. EVALS\002,i8,9x,\002GRAD. EVALS\002,i8/\002 PRELDF\002,d16.3,6x,\
\002NPRELDF\002,d15.3)";
    static char fmt_470[] = "(/\002     I      FINAL X(I)\002,8x,\002D(I)\
\002,10x,\002G(I)\002/)";
    static char fmt_490[] = "(1x,i5,d16.6,2d14.3)";
    static char fmt_510[] = "(/\002 INCONSISTENT DIMENSIONS\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, m, nf, ng, ol, pu, iv1, alg;
    static doublereal oldf, reldf, nreldf, preldf;

    /* Fortran I/O blocks */
    static cilist io___61 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___77 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___82 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___89 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_410, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_420, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_450, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_470, 0 };
    static cilist io___95 = { 0, 0, 0, fmt_490, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_510, 0 };



/*  ***  PRINT ITERATION SUMMARY FOR ***SOL (VERSION 2.3)  *** */

/*  ***  PARAMETER DECLARATIONS  *** */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */

/* /6S */
/*     REAL MODEL1(6), MODEL2(6) */
/* /7S */
/* / */

/*  ***  NO EXTERNAL FUNCTIONS OR SUBROUTINES  *** */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA ALGSAV/51/, NEEDHD/36/, NFCALL/6/, NFCOV/52/, NGCALL/30/, */
/*    1     NGCOV/53/, NITER/31/, OUTLEV/19/, PRNTIT/39/, PRUNIT/21/, */
/*    2     SOLPRT/22/, STATPR/23/, SUSED/64/, X0PRT/24/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DSTNRM/2/, F/10/, F0/13/, FDIF/11/, NREDUC/6/, PREDUC/7/, */
/*    1     RELDX/17/, STPPAR/5/ */
/* /7 */
/* / */

/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */
/* /6S */
/*     DATA MODEL1(1)/4H    /, MODEL1(2)/4H    /, MODEL1(3)/4H    /, */
/*    1     MODEL1(4)/4H    /, MODEL1(5)/4H  G /, MODEL1(6)/4H  S /, */
/*    2     MODEL2(1)/4H G  /, MODEL2(2)/4H S  /, MODEL2(3)/4HG-S /, */
/*    3     MODEL2(4)/4HS-G /, MODEL2(5)/4H-S-G/, MODEL2(6)/4H-G-S/ */
/* /7S */
    /* Parameter adjustments */
    --iv;
    --v;
    --x;
    --g;
    --d__;

    /* Function Body */
/* / */

/* -------------------------------  BODY  -------------------------------- */

    pu = iv[21];
    if (pu == 0) {
	goto L999;
    }
    iv1 = iv[1];
    if (iv1 > 62) {
	iv1 += -51;
    }
    ol = iv[19];
    alg = (iv[51] - 1) % 2 + 1;
    if (iv1 < 2 || iv1 > 15) {
	goto L370;
    }
    if (iv1 >= 12) {
	goto L120;
    }
    if (iv1 == 2 && iv[31] == 0) {
	goto L390;
    }
    if (ol == 0) {
	goto L120;
    }
    if (iv1 >= 10 && iv[39] == 0) {
	goto L120;
    }
    if (iv1 > 2) {
	goto L10;
    }
    ++iv[39];
    if (iv[39] < abs(ol)) {
	goto L999;
    }
L10:
    nf = iv[6] - abs(iv[52]);
    iv[39] = 0;
    reldf = 0.;
    preldf = 0.;
/* Computing MAX */
    d__1 = abs(v[13]), d__2 = abs(v[10]);
    oldf = max(d__1,d__2);
    if (oldf <= 0.) {
	goto L20;
    }
    reldf = v[11] / oldf;
    preldf = v[7] / oldf;
L20:
    if (ol > 0) {
	goto L60;
    }

/*        ***  PRINT SHORT SUMMARY LINE  *** */

    if (iv[36] == 1 && alg == 1) {
	io___61.ciunit = pu;
	s_wsfe(&io___61);
	e_wsfe();
    }
    if (iv[36] == 1 && alg == 2) {
	io___62.ciunit = pu;
	s_wsfe(&io___62);
	e_wsfe();
    }
    iv[36] = 0;
    if (alg == 2) {
	goto L50;
    }
    m = iv[64];
    io___64.ciunit = pu;
    s_wsfe(&io___64);
    do_fio(&c__1, (char *)&iv[31], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[17], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, model1 + (m - 1 << 2), (ftnlen)4);
    do_fio(&c__1, model2 + (m - 1 << 2), (ftnlen)4);
    do_fio(&c__1, (char *)&v[5], (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L120;

L50:
    io___65.ciunit = pu;
    s_wsfe(&io___65);
    do_fio(&c__1, (char *)&iv[31], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[17], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[5], (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L120;

/*     ***  PRINT LONG SUMMARY LINE  *** */

L60:
    if (iv[36] == 1 && alg == 1) {
	io___66.ciunit = pu;
	s_wsfe(&io___66);
	e_wsfe();
    }
    if (iv[36] == 1 && alg == 2) {
	io___67.ciunit = pu;
	s_wsfe(&io___67);
	e_wsfe();
    }
    iv[36] = 0;
    nreldf = 0.;
    if (oldf > 0.) {
	nreldf = v[6] / oldf;
    }
    if (alg == 2) {
	goto L90;
    }
    m = iv[64];
    io___69.ciunit = pu;
    s_wsfe(&io___69);
    do_fio(&c__1, (char *)&iv[31], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[17], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, model1 + (m - 1 << 2), (ftnlen)4);
    do_fio(&c__1, model2 + (m - 1 << 2), (ftnlen)4);
    do_fio(&c__1, (char *)&v[5], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[2], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nreldf, (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L120;

L90:
    io___70.ciunit = pu;
    s_wsfe(&io___70);
    do_fio(&c__1, (char *)&iv[31], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[17], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[5], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[2], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nreldf, (ftnlen)sizeof(doublereal));
    e_wsfe();

L120:
    if (iv1 <= 2) {
	goto L999;
    }
    i__ = iv[23];
    if (i__ == -1) {
	goto L460;
    }
    if (i__ + iv1 < 0) {
	goto L460;
    }
    switch (iv1) {
	case 1:  goto L999;
	case 2:  goto L999;
	case 3:  goto L130;
	case 4:  goto L150;
	case 5:  goto L170;
	case 6:  goto L190;
	case 7:  goto L210;
	case 8:  goto L230;
	case 9:  goto L250;
	case 10:  goto L270;
	case 11:  goto L290;
	case 12:  goto L310;
	case 13:  goto L330;
	case 14:  goto L350;
	case 15:  goto L500;
    }

L130:
    io___72.ciunit = pu;
    s_wsfe(&io___72);
    e_wsfe();
    goto L430;

L150:
    io___73.ciunit = pu;
    s_wsfe(&io___73);
    e_wsfe();
    goto L430;

L170:
    io___74.ciunit = pu;
    s_wsfe(&io___74);
    e_wsfe();
    goto L430;

L190:
    io___75.ciunit = pu;
    s_wsfe(&io___75);
    e_wsfe();
    goto L430;

L210:
    io___76.ciunit = pu;
    s_wsfe(&io___76);
    e_wsfe();
    goto L430;

L230:
    io___77.ciunit = pu;
    s_wsfe(&io___77);
    e_wsfe();
    goto L430;

L250:
    io___78.ciunit = pu;
    s_wsfe(&io___78);
    e_wsfe();
    goto L430;

L270:
    io___79.ciunit = pu;
    s_wsfe(&io___79);
    e_wsfe();
    goto L430;

L290:
    io___80.ciunit = pu;
    s_wsfe(&io___80);
    e_wsfe();
    goto L430;

L310:
    io___81.ciunit = pu;
    s_wsfe(&io___81);
    e_wsfe();

    goto L390;

L330:
    io___82.ciunit = pu;
    s_wsfe(&io___82);
    e_wsfe();
    goto L999;

L350:
    io___83.ciunit = pu;
    s_wsfe(&io___83);
    e_wsfe();
    if (iv[31] > 0) {
	goto L460;
    }
    goto L390;

L370:
    io___84.ciunit = pu;
    s_wsfe(&io___84);
    do_fio(&c__1, (char *)&iv[1], (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;

/*  ***  INITIAL CALL ON DITSUM  *** */

L390:
    if (iv[24] != 0) {
	io___85.ciunit = pu;
	s_wsfe(&io___85);
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&d__[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
/*     *** THE FOLLOWING ARE TO AVOID UNDEFINED VARIABLES WHEN THE */
/*     *** FUNCTION EVALUATION LIMIT IS 1... */
    v[2] = 0.;
    v[11] = 0.;
    v[6] = 0.;
    v[7] = 0.;
    v[17] = 0.;
    if (iv1 >= 12) {
	goto L999;
    }
    iv[36] = 0;
    iv[39] = 0;
    if (ol == 0) {
	goto L999;
    }
    if (ol < 0 && alg == 1) {
	io___86.ciunit = pu;
	s_wsfe(&io___86);
	e_wsfe();
    }
    if (ol < 0 && alg == 2) {
	io___87.ciunit = pu;
	s_wsfe(&io___87);
	e_wsfe();
    }
    if (ol > 0 && alg == 1) {
	io___88.ciunit = pu;
	s_wsfe(&io___88);
	e_wsfe();
    }
    if (ol > 0 && alg == 2) {
	io___89.ciunit = pu;
	s_wsfe(&io___89);
	e_wsfe();
    }
    if (alg == 1) {
	io___90.ciunit = pu;
	s_wsfe(&io___90);
	do_fio(&c__1, (char *)&iv[6], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (alg == 2) {
	io___91.ciunit = pu;
	s_wsfe(&io___91);
	do_fio(&c__1, (char *)&iv[6], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    goto L999;

/*  ***  PRINT VARIOUS INFORMATION REQUESTED ON SOLUTION  *** */

L430:
    iv[36] = 1;
    if (iv[23] <= 0) {
	goto L460;
    }
/* Computing MAX */
    d__1 = abs(v[13]), d__2 = abs(v[10]);
    oldf = max(d__1,d__2);
    preldf = 0.;
    nreldf = 0.;
    if (oldf <= 0.) {
	goto L440;
    }
    preldf = v[7] / oldf;
    nreldf = v[6] / oldf;
L440:
    nf = iv[6] - iv[52];
    ng = iv[30] - iv[53];
    io___93.ciunit = pu;
    s_wsfe(&io___93);
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[17], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ng, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nreldf, (ftnlen)sizeof(doublereal));
    e_wsfe();

L460:
    if (iv[22] == 0) {
	goto L999;
    }
    iv[36] = 1;
    if (iv[51] > 2) {
	goto L999;
    }
    io___94.ciunit = pu;
    s_wsfe(&io___94);
    e_wsfe();
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L480: */
	io___95.ciunit = pu;
	s_wsfe(&io___95);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&d__[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    goto L999;

L500:
    io___96.ciunit = pu;
    s_wsfe(&io___96);
    e_wsfe();
L999:
    return 0;
/*  ***  LAST CARD OF DITSUM FOLLOWS  *** */
} /* ditsum_ */

/* Subroutine */ int divset_(integer *alg, integer *iv, integer *liv, integer 
	*lv, doublereal *v)
{
    /* Initialized data */

    static integer miniv[4] = { 82,59,103,103 };
    static integer minv[4] = { 98,71,101,85 };

    static integer mv, miv, alg1;
    extern integer i7mdcn_(integer *);
    extern /* Subroutine */ int dv7dfl_(integer *, integer *, doublereal *);


/*  ***  SUPPLY ***SOL (VERSION 2.3) DEFAULT VALUES TO IV AND V  *** */

/*  ***  ALG = 1 MEANS REGRESSION CONSTANTS. */
/*  ***  ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS. */


/* I7MDCN... RETURNS MACHINE-DEPENDENT INTEGER CONSTANTS. */
/* DV7DFL.... PROVIDES DEFAULT VALUES TO V. */


/*  ***  SUBSCRIPTS FOR IV  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA ALGSAV/51/, COVPRT/14/, COVREQ/15/, DRADPR/101/, DTYPE/16/, */
/*    1     HC/71/, IERR/75/, INITH/25/, INITS/25/, IPIVOT/76/, */
/*    2     IVNEED/3/, LASTIV/44/, LASTV/45/, LMAT/42/, MXFCAL/17/, */
/*    3     MXITER/18/, NFCOV/52/, NGCOV/53/, NVDFLT/50/, NVSAVE/9/, */
/*    4     OUTLEV/19/, PARPRT/20/, PARSAV/49/, PERM/58/, PRUNIT/21/, */
/*    5     QRTYP/80/, RDREQ/57/, RMAT/78/, SOLPRT/22/, STATPR/23/, */
/*    6     VNEED/4/, VSAVE/60/, X0PRT/24/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --iv;
    --v;

    /* Function Body */

/* -------------------------------  BODY  -------------------------------- */

    if (21 <= *liv) {
	iv[21] = i7mdcn_(&c__1);
    }
    if (51 <= *liv) {
	iv[51] = *alg;
    }
    if (*alg < 1 || *alg > 4) {
	goto L40;
    }
    miv = miniv[*alg - 1];
    if (*liv < miv) {
	goto L20;
    }
    mv = minv[*alg - 1];
    if (*lv < mv) {
	goto L30;
    }
    alg1 = (*alg - 1) % 2 + 1;
    dv7dfl_(&alg1, lv, &v[1]);
    iv[1] = 12;
    if (*alg > 2) {
	iv[101] = 1;
    }
    iv[3] = 0;
    iv[44] = miv;
    iv[45] = mv;
    iv[42] = mv + 1;
    iv[17] = 200;
    iv[18] = 150;
    iv[19] = 1;
    iv[20] = 1;
    iv[58] = miv + 1;
    iv[22] = 1;
    iv[23] = 1;
    iv[4] = 0;
    iv[24] = 1;

    if (alg1 >= 2) {
	goto L10;
    }

/*  ***  REGRESSION  VALUES */

    iv[14] = 3;
    iv[15] = 1;
    iv[16] = 1;
    iv[71] = 0;
    iv[75] = 0;
    iv[25] = 0;
    iv[76] = 0;
    iv[50] = 32;
    iv[60] = 58;
    if (*alg > 2) {
	iv[60] += 3;
    }
    iv[49] = iv[60] + 9;
    iv[80] = 1;
    iv[57] = 3;
    iv[78] = 0;
    goto L999;

/*  ***  GENERAL OPTIMIZATION VALUES */

L10:
    iv[16] = 0;
    iv[25] = 1;
    iv[52] = 0;
    iv[53] = 0;
    iv[50] = 25;
    iv[49] = 47;
    if (*alg > 2) {
	iv[49] = 61;
    }
    goto L999;

L20:
    iv[1] = 15;
    goto L999;

L30:
    iv[1] = 16;
    goto L999;

L40:
    iv[1] = 67;

L999:
    return 0;
/*  ***  LAST CARD OF DIVSET FOLLOWS  *** */
} /* divset_ */

/* Subroutine */ int dl7itv_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, i0, ii, ij;
    static doublereal xi;
    static integer im1, np1;


/*  ***  SOLVE  (L**T)*X = Y,  WHERE  L  IS AN  N X N  LOWER TRIANGULAR */
/*  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME */
/*  ***  STORAGE.  *** */

/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = y[i__];
    }
    np1 = *n + 1;
    i0 = *n * (*n + 1) / 2;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	xi = x[i__] / l[i0];
	x[i__] = xi;
	if (i__ <= 1) {
	    goto L999;
	}
	i0 -= i__;
	if (xi == 0.) {
	    goto L30;
	}
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    ij = i0 + j;
	    x[j] -= xi * l[ij];
/* L20: */
	}
L30:
	;
    }
L999:
    return 0;
/*  ***  LAST CARD OF DL7ITV FOLLOWS  *** */
} /* dl7itv_ */

/* Subroutine */ int dl7ivm_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);


/*  ***  SOLVE  L*X = Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR */
/*  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME */
/*  ***  STORAGE.  *** */

/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (y[k] != 0.) {
	    goto L20;
	}
	x[k] = 0.;
/* L10: */
    }
    goto L999;
L20:
    j = k * (k + 1) / 2;
    x[k] = y[k] / l[j];
    if (k >= *n) {
	goto L999;
    }
    ++k;
    i__1 = *n;
    for (i__ = k; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	t = dd7tpr_(&i__2, &l[j + 1], &x[1]);
	j += i__;
	x[i__] = (y[i__] - t) / l[j];
/* L30: */
    }
L999:
    return 0;
/*  ***  LAST CARD OF DL7IVM FOLLOWS  *** */
} /* dl7ivm_ */

/* Subroutine */ int dl7tvm_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, i0, ij;
    static doublereal yi;


/*  ***  COMPUTE  X = (L**T)*Y, WHERE  L  IS AN  N X N  LOWER */
/*  ***  TRIANGULAR MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY */
/*  ***  OCCUPY THE SAME STORAGE.  *** */

/*     DIMENSION L(N*(N+1)/2) */
/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
    i0 = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yi = y[i__];
	x[i__] = 0.;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ij = i0 + j;
	    x[j] += yi * l[ij];
/* L10: */
	}
	i0 += i__;
/* L20: */
    }
/* L999: */
    return 0;
/*  ***  LAST CARD OF DL7TVM FOLLOWS  *** */
} /* dl7tvm_ */

/* Subroutine */ int dl7upd_(doublereal *beta, doublereal *gamma, doublereal *
	l, doublereal *lambda, doublereal *lplus, integer *n, doublereal *w, 
	doublereal *z__)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b;
    static integer i__, j, k;
    static doublereal s, bj, gj;
    static integer ij, jj;
    static doublereal lj, wj, nu, zj;
    static integer jp1, nm1, np1;
    static doublereal eta, lij, ljj, theta;


/*  ***  COMPUTE LPLUS = SECANT UPDATE OF L  *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION L(N*(N+1)/2), LPLUS(N*(N+1)/2) */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/*   BETA = SCRATCH VECTOR. */
/*  GAMMA = SCRATCH VECTOR. */
/*      L (INPUT) LOWER TRIANGULAR MATRIX, STORED ROWWISE. */
/* LAMBDA = SCRATCH VECTOR. */
/*  LPLUS (OUTPUT) LOWER TRIANGULAR MATRIX, STORED ROWWISE, WHICH MAY */
/*             OCCUPY THE SAME STORAGE AS  L. */
/*      N (INPUT) LENGTH OF VECTOR PARAMETERS AND ORDER OF MATRICES. */
/*      W (INPUT, DESTROYED ON OUTPUT) RIGHT SINGULAR VECTOR OF RANK 1 */
/*             CORRECTION TO  L. */
/*      Z (INPUT, DESTROYED ON OUTPUT) LEFT SINGULAR VECTOR OF RANK 1 */
/*             CORRECTION TO  L. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  APPLICATION AND USAGE RESTRICTIONS  *** */

/*        THIS ROUTINE UPDATES THE CHOLESKY FACTOR  L  OF A SYMMETRIC */
/*     POSITIVE DEFINITE MATRIX TO WHICH A SECANT UPDATE IS BEING */
/*     APPLIED -- IT COMPUTES A CHOLESKY FACTOR  LPLUS  OF */
/*     L * (I + Z*W**T) * (I + W*Z**T) * L**T.  IT IS ASSUMED THAT  W */
/*     AND  Z  HAVE BEEN CHOSEN SO THAT THE UPDATED MATRIX IS STRICTLY */
/*     POSITIVE DEFINITE. */

/*  ***  ALGORITHM NOTES  *** */

/*        THIS CODE USES RECURRENCE 3 OF REF. 1 (WITH D(J) = 1 FOR ALL J) */
/*     TO COMPUTE  LPLUS  OF THE FORM  L * (I + Z*W**T) * Q,  WHERE  Q */
/*     IS AN ORTHOGONAL MATRIX THAT MAKES THE RESULT LOWER TRIANGULAR. */
/*        LPLUS MAY HAVE SOME NEGATIVE DIAGONAL ELEMENTS. */

/*  ***  REFERENCES  *** */

/* 1.  GOLDFARB, D. (1976), FACTORIZED VARIABLE METRIC METHODS FOR UNCON- */
/*             STRAINED OPTIMIZATION, MATH. COMPUT. 30, PP. 796-811. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (FALL 1979). */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED */
/*     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND */
/*     MCS-7906671. */

/* ------------------------  EXTERNAL QUANTITIES  ------------------------ */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/* --------------------------  LOCAL VARIABLES  -------------------------- */


/*  ***  DATA INITIALIZATIONS  *** */

/* /6 */
/*     DATA ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --l;
    --lplus;
    --z__;
    --w;
    --lambda;
    --gamma;
    --beta;

    /* Function Body */
    nu = 1.;
    eta = 0.;
    if (*n <= 1) {
	goto L30;
    }
    nm1 = *n - 1;

/*  ***  TEMPORARILY STORE S(J) = SUM OVER K = J+1 TO N OF W(K)**2 IN */
/*  ***  LAMBDA(J). */

    s = 0.;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n - i__;
/* Computing 2nd power */
	d__1 = w[j + 1];
	s += d__1 * d__1;
	lambda[j] = s;
/* L10: */
    }

/*  ***  COMPUTE LAMBDA, GAMMA, AND BETA BY GOLDFARB*S RECURRENCE 3. */

    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	wj = w[j];
	a = nu * z__[j] - eta * wj;
	theta = a * wj + 1.;
	s = a * lambda[j];
/* Computing 2nd power */
	d__1 = theta;
	lj = sqrt(d__1 * d__1 + a * s);
	if (theta > 0.) {
	    lj = -lj;
	}
	lambda[j] = lj;
	b = theta * wj + s;
	gamma[j] = b * nu / lj;
	beta[j] = (a - b * eta) / lj;
	nu = -nu / lj;
/* Computing 2nd power */
	d__1 = a;
	eta = -(eta + d__1 * d__1 / (theta - lj)) / lj;
/* L20: */
    }
L30:
    lambda[*n] = (nu * z__[*n] - eta * w[*n]) * w[*n] + 1.;

/*  ***  UPDATE L, GRADUALLY OVERWRITING  W  AND  Z  WITH  L*W  AND  L*Z. */

    np1 = *n + 1;
    jj = *n * (*n + 1) / 2;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	j = np1 - k;
	lj = lambda[j];
	ljj = l[jj];
	lplus[jj] = lj * ljj;
	wj = w[j];
	w[j] = ljj * wj;
	zj = z__[j];
	z__[j] = ljj * zj;
	if (k == 1) {
	    goto L50;
	}
	bj = beta[j];
	gj = gamma[j];
	ij = jj + j;
	jp1 = j + 1;
	i__2 = *n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    lij = l[ij];
	    lplus[ij] = lj * lij + bj * w[i__] + gj * z__[i__];
	    w[i__] += lij * wj;
	    z__[i__] += lij * zj;
	    ij += i__;
/* L40: */
	}
L50:
	jj -= j;
/* L60: */
    }

/* L999: */
    return 0;
/*  ***  LAST CARD OF DL7UPD FOLLOWS  *** */
} /* dl7upd_ */

/* Subroutine */ int dl7vml_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer i0, ii, ij, np1;


/*  ***  COMPUTE  X = L*Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR */
/*  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME */
/*  ***  STORAGE.  *** */

/*     DIMENSION L(N*(N+1)/2) */
/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
    np1 = *n + 1;
    i0 = *n * (*n + 1) / 2;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	i0 -= i__;
	t = 0.;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ij = i0 + j;
	    t += l[ij] * y[j];
/* L10: */
	}
	x[i__] = t;
/* L20: */
    }
/* L999: */
    return 0;
/*  ***  LAST CARD OF DL7VML FOLLOWS  *** */
} /* dl7vml_ */

/* Subroutine */ int dparck_(integer *alg, doublereal *d__, integer *iv, 
	integer *liv, integer *lv, integer *n, doublereal *v)
{
    /* Initialized data */

    static doublereal big = 0.;
    static doublereal machep = -1.;
    static doublereal tiny = 1.;
    static doublereal zero = 0.;
    static char vn[4*2*34+1] = "EPSLON..PHMNFC..PHMXFC..DECFAC..INCFAC..RDFC\
MN..RDFCMX..TUNER1..TUNER2..TUNER3..TUNER4..TUNER5..AFCTOL..RFCTOL..XCTOL...\
XFTOL...LMAX0...LMAXS...SCTOL...DINIT...DTINIT..D0INIT..DFAC....DLTFDC..DLTF\
DJ..DELTA0..FUZZ....RLIMIT..COSMIN..HUBERC..RSPTOL..SIGMIN..ETA0....BIAS....";
    static doublereal vm[34] = { .001,-.99,.001,.01,1.2,.01,1.2,0.,0.,.001,
	    -1.,0.0,0.,0.0,0.,0.,0.0,0.0,0.,-10.,0.,0.,0.,0.0,0.0,0.0,1.01,
	    1e10,0.0,0.,0.,0.,0.0,0. };
    static doublereal vx[34] = { .9,-.001,10.,.8,100.,.8,100.,.5,.5,1.,1.,0.0,
	    0.0,.1,1.,1.,0.0,0.0,1.,0.0,0.0,0.0,1.,1.,1.,1.,1e10,0.0,1.,0.0,
	    1.,1.,1.,1. };
    static char varnm[1*2+1] = "PP";
    static char sh[1*2+1] = "SH";
    static char cngd[4*3+1] = "---CHANGED V";
    static char dflt[4*3+1] = "NONDEFAULT V";
    static integer ijmp = 33;
    static integer jlim[4] = { 0,24,0,24 };
    static integer ndflt[4] = { 32,25,32,25 };
    static integer miniv[4] = { 82,59,103,103 };

    /* Format strings */
    static char fmt_10[] = "(/\002 THE FIRST PARAMETER TO DIVSET SHOULD B\
E\002,i3,\002 RATHER THAN\002,i3)";
    static char fmt_40[] = "(/\002 /// BAD\002,a1,\002 =\002,i5)";
    static char fmt_70[] = "(/\002 /// \002,1a1,\002 CHANGED FROM \002,i5\
,\002 TO \002,i5)";
    static char fmt_90[] = "(/\002 ///  IV(1) =\002,i5,\002 SHOULD BE BETWEE\
N 0 AND 14.\002)";
    static char fmt_130[] = "(/\002 ///  \002,2a4,\002.. V(\002,i2,\002) \
=\002,d11.3,\002 SHOULD\002,\002 BE BETWEEN\002,d11.3,\002 AND\002,d11.3)";
    static char fmt_160[] = "(/\002 IV(NVDFLT) =\002,i5,\002 RATHER THAN \
\002,i5)";
    static char fmt_180[] = "(/\002 ///  D(\002,i3,\002) =\002,d11.3,\002 SH\
OULD BE POSITIVE\002)";
    static char fmt_220[] = "(/\002 NONDEFAULT VALUES....\002/\002 INIT\002,\
a1,\002..... IV(25) =\002,i3)";
    static char fmt_260[] = "(/\002 \002,3a4,\002ALUES....\002/)";
    static char fmt_240[] = "(\002 DTYPE..... IV(16) =\002,i3)";
    static char fmt_270[] = "(1x,2a4,\002.. V(\002,i2,\002) =\002,d15.7)";
    static char fmt_310[] = "(/\002 /// LIV =\002,i5,\002 MUST BE AT LEAS\
T\002,i5)";
    static char fmt_330[] = "(/\002 /// LV =\002,i5,\002 MUST BE AT LEAST\
\002,i5)";
    static char fmt_350[] = "(/\002 /// ALG =\002,i5,\002 MUST BE 1 2, 3, OR\
 4\002)";
    static char fmt_370[] = "(/\002 /// LIV =\002,i5,\002 MUST BE AT LEAS\
T\002,i5,\002 TO COMPUTE TRUE MIN. LIV AND MIN. LV\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, l, m, ii;
    static doublereal vk;
    static integer pu, iv1, alg1, miv1, miv2;
    static char which[4*3];
    extern doublereal dr7mdc_(integer *);
    extern /* Subroutine */ int dv7dfl_(integer *, integer *, doublereal *), 
	    dv7cpy_(integer *, doublereal *, doublereal *);
    static integer parsv1, ndfalt;
    extern /* Subroutine */ int divset_(integer *, integer *, integer *, 
	    integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___163 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___168 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___171 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___172 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___179 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___180 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___181 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___182 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___183 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___184 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___186 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___187 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___189 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___190 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___191 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___192 = { 0, 0, 0, fmt_370, 0 };



/*  ***  CHECK ***SOL (VERSION 2.3) PARAMETERS, PRINT CHANGED VALUES  *** */

/*  ***  ALG = 1 FOR REGRESSION, ALG = 2 FOR GENERAL UNCONSTRAINED OPT. */


/* DIVSET  -- SUPPLIES DEFAULT VALUES TO BOTH IV AND V. */
/* DR7MDC -- RETURNS MACHINE-DEPENDENT CONSTANTS. */
/* DV7CPY  -- COPIES ONE VECTOR TO ANOTHER. */
/* DV7DFL  -- SUPPLIES DEFAULT PARAMETER VALUES TO V ALONE. */

/*  ***  LOCAL VARIABLES  *** */

/* /6S */
/*     INTEGER VARNM(2), SH(2) */
/*     REAL CNGD(3), DFLT(3), VN(2,34), WHICH(3) */
/* /7S */
/* / */

/*  ***  IV AND V SUBSCRIPTS  *** */



/* /6 */
/*     DATA ALGSAV/51/, DINIT/38/, DTYPE/16/, DTYPE0/54/, EPSLON/19/, */
/*    1     INITS/25/, IVNEED/3/, LASTIV/44/, LASTV/45/, LMAT/42/, */
/*    2     NEXTIV/46/, NEXTV/47/, NVDFLT/50/, OLDN/38/, PARPRT/20/, */
/*    3     PARSAV/49/, PERM/58/, PRUNIT/21/, VNEED/4/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --iv;
    --v;
    --d__;

    /* Function Body */
/* /6S */
/*     DATA VN(1,1),VN(2,1)/4HEPSL,4HON../ */
/*     DATA VN(1,2),VN(2,2)/4HPHMN,4HFC../ */
/*     DATA VN(1,3),VN(2,3)/4HPHMX,4HFC../ */
/*     DATA VN(1,4),VN(2,4)/4HDECF,4HAC../ */
/*     DATA VN(1,5),VN(2,5)/4HINCF,4HAC../ */
/*     DATA VN(1,6),VN(2,6)/4HRDFC,4HMN../ */
/*     DATA VN(1,7),VN(2,7)/4HRDFC,4HMX../ */
/*     DATA VN(1,8),VN(2,8)/4HTUNE,4HR1../ */
/*     DATA VN(1,9),VN(2,9)/4HTUNE,4HR2../ */
/*     DATA VN(1,10),VN(2,10)/4HTUNE,4HR3../ */
/*     DATA VN(1,11),VN(2,11)/4HTUNE,4HR4../ */
/*     DATA VN(1,12),VN(2,12)/4HTUNE,4HR5../ */
/*     DATA VN(1,13),VN(2,13)/4HAFCT,4HOL../ */
/*     DATA VN(1,14),VN(2,14)/4HRFCT,4HOL../ */
/*     DATA VN(1,15),VN(2,15)/4HXCTO,4HL.../ */
/*     DATA VN(1,16),VN(2,16)/4HXFTO,4HL.../ */
/*     DATA VN(1,17),VN(2,17)/4HLMAX,4H0.../ */
/*     DATA VN(1,18),VN(2,18)/4HLMAX,4HS.../ */
/*     DATA VN(1,19),VN(2,19)/4HSCTO,4HL.../ */
/*     DATA VN(1,20),VN(2,20)/4HDINI,4HT.../ */
/*     DATA VN(1,21),VN(2,21)/4HDTIN,4HIT../ */
/*     DATA VN(1,22),VN(2,22)/4HD0IN,4HIT../ */
/*     DATA VN(1,23),VN(2,23)/4HDFAC,4H..../ */
/*     DATA VN(1,24),VN(2,24)/4HDLTF,4HDC../ */
/*     DATA VN(1,25),VN(2,25)/4HDLTF,4HDJ../ */
/*     DATA VN(1,26),VN(2,26)/4HDELT,4HA0../ */
/*     DATA VN(1,27),VN(2,27)/4HFUZZ,4H..../ */
/*     DATA VN(1,28),VN(2,28)/4HRLIM,4HIT../ */
/*     DATA VN(1,29),VN(2,29)/4HCOSM,4HIN../ */
/*     DATA VN(1,30),VN(2,30)/4HHUBE,4HRC../ */
/*     DATA VN(1,31),VN(2,31)/4HRSPT,4HOL../ */
/*     DATA VN(1,32),VN(2,32)/4HSIGM,4HIN../ */
/*     DATA VN(1,33),VN(2,33)/4HETA0,4H..../ */
/*     DATA VN(1,34),VN(2,34)/4HBIAS,4H..../ */
/* /7S */
/* / */


/* /6S */
/*     DATA VARNM(1)/1HP/, VARNM(2)/1HP/, SH(1)/1HS/, SH(2)/1HH/ */
/*     DATA CNGD(1),CNGD(2),CNGD(3)/4H---C,4HHANG,4HED V/, */
/*    1     DFLT(1),DFLT(2),DFLT(3)/4HNOND,4HEFAU,4HLT V/ */
/* /7S */
/* / */

/* ...............................  BODY  ................................ */

    pu = 0;
    if (21 <= *liv) {
	pu = iv[21];
    }
    if (51 > *liv) {
	goto L20;
    }
    if (*alg == iv[51]) {
	goto L20;
    }
    if (pu != 0) {
	io___163.ciunit = pu;
	s_wsfe(&io___163);
	do_fio(&c__1, (char *)&(*alg), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iv[51], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    iv[1] = 67;
    goto L999;
L20:
    if (*alg < 1 || *alg > 4) {
	goto L340;
    }
    miv1 = miniv[*alg - 1];
    if (iv[1] == 15) {
	goto L360;
    }
    alg1 = (*alg - 1) % 2 + 1;
    if (iv[1] == 0) {
	divset_(alg, &iv[1], liv, lv, &v[1]);
    }
    iv1 = iv[1];
    if (iv1 != 13 && iv1 != 12) {
	goto L30;
    }
    if (58 <= *liv) {
/* Computing MAX */
	i__1 = miv1, i__2 = iv[58] - 1;
	miv1 = max(i__1,i__2);
    }
    if (3 <= *liv) {
	miv2 = miv1 + max(iv[3],0);
    }
    if (44 <= *liv) {
	iv[44] = miv2;
    }
    if (*liv < miv1) {
	goto L300;
    }
    iv[3] = 0;
    iv[45] = max(iv[4],0) + iv[42] - 1;
    iv[4] = 0;
    if (*liv < miv2) {
	goto L300;
    }
    if (*lv < iv[45]) {
	goto L320;
    }
L30:
    if (iv1 < 12 || iv1 > 14) {
	goto L60;
    }
    if (*n >= 1) {
	goto L50;
    }
    iv[1] = 81;
    if (pu == 0) {
	goto L999;
    }
    io___168.ciunit = pu;
    s_wsfe(&io___168);
    do_fio(&c__1, varnm + (alg1 - 1), (ftnlen)1);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;
L50:
    if (iv1 != 14) {
	iv[46] = iv[58];
    }
    if (iv1 != 14) {
	iv[47] = iv[42];
    }
    if (iv1 == 13) {
	goto L999;
    }
    k = iv[49] - 19;
    i__1 = *lv - k;
    dv7dfl_(&alg1, &i__1, &v[k + 1]);
    iv[54] = 2 - alg1;
    iv[38] = *n;
    s_copy(which, dflt, (ftnlen)4, (ftnlen)4);
    s_copy(which + 4, dflt + 4, (ftnlen)4, (ftnlen)4);
    s_copy(which + 8, dflt + 8, (ftnlen)4, (ftnlen)4);
    goto L110;
L60:
    if (*n == iv[38]) {
	goto L80;
    }
    iv[1] = 17;
    if (pu == 0) {
	goto L999;
    }
    io___171.ciunit = pu;
    s_wsfe(&io___171);
    do_fio(&c__1, varnm + (alg1 - 1), (ftnlen)1);
    do_fio(&c__1, (char *)&iv[38], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;

L80:
    if (iv1 <= 11 && iv1 >= 1) {
	goto L100;
    }
    iv[1] = 80;
    if (pu != 0) {
	io___172.ciunit = pu;
	s_wsfe(&io___172);
	do_fio(&c__1, (char *)&iv1, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L999;

L100:
    s_copy(which, cngd, (ftnlen)4, (ftnlen)4);
    s_copy(which + 4, cngd + 4, (ftnlen)4, (ftnlen)4);
    s_copy(which + 8, cngd + 8, (ftnlen)4, (ftnlen)4);

L110:
    if (iv1 == 14) {
	iv1 = 12;
    }
    if (big > tiny) {
	goto L120;
    }
    tiny = dr7mdc_(&c__1);
    machep = dr7mdc_(&c__3);
    big = dr7mdc_(&c__6);
    vm[11] = machep;
    vx[11] = big;
    vx[12] = big;
    vm[13] = machep;
    vm[16] = tiny;
    vx[16] = big;
    vm[17] = tiny;
    vx[17] = big;
    vx[19] = big;
    vx[20] = big;
    vx[21] = big;
    vm[23] = machep;
    vm[24] = machep;
    vm[25] = machep;
    vx[27] = dr7mdc_(&c__5);
    vm[28] = machep;
    vx[29] = big;
    vm[32] = machep;
L120:
    m = 0;
    i__ = 1;
    j = jlim[alg1 - 1];
    k = 19;
    ndfalt = ndflt[alg1 - 1];
    i__1 = ndfalt;
    for (l = 1; l <= i__1; ++l) {
	vk = v[k];
	if (vk >= vm[i__ - 1] && vk <= vx[i__ - 1]) {
	    goto L140;
	}
	m = k;
	if (pu != 0) {
	    io___179.ciunit = pu;
	    s_wsfe(&io___179);
	    do_fio(&c__1, vn + ((i__ << 1) - 2 << 2), (ftnlen)4);
	    do_fio(&c__1, vn + ((i__ << 1) - 1 << 2), (ftnlen)4);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&vk, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&vm[i__ - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&vx[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
L140:
	++k;
	++i__;
	if (i__ == j) {
	    i__ = ijmp;
	}
/* L150: */
    }

    if (iv[50] == ndfalt) {
	goto L170;
    }
    iv[1] = 51;
    if (pu == 0) {
	goto L999;
    }
    io___180.ciunit = pu;
    s_wsfe(&io___180);
    do_fio(&c__1, (char *)&iv[50], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ndfalt, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;
L170:
    if ((iv[16] > 0 || v[38] > zero) && iv1 == 12) {
	goto L200;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] > zero) {
	    goto L190;
	}
	m = 18;
	if (pu != 0) {
	    io___181.ciunit = pu;
	    s_wsfe(&io___181);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&d__[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
L190:
	;
    }
L200:
    if (m == 0) {
	goto L210;
    }
    iv[1] = m;
    goto L999;

L210:
    if (pu == 0 || iv[20] == 0) {
	goto L999;
    }
    if (iv1 != 12 || iv[25] == alg1 - 1) {
	goto L230;
    }
    m = 1;
    io___182.ciunit = pu;
    s_wsfe(&io___182);
    do_fio(&c__1, sh + (alg1 - 1), (ftnlen)1);
    do_fio(&c__1, (char *)&iv[25], (ftnlen)sizeof(integer));
    e_wsfe();
L230:
    if (iv[16] == iv[54]) {
	goto L250;
    }
    if (m == 0) {
	io___183.ciunit = pu;
	s_wsfe(&io___183);
	do_fio(&c__3, which, (ftnlen)4);
	e_wsfe();
    }
    m = 1;
    io___184.ciunit = pu;
    s_wsfe(&io___184);
    do_fio(&c__1, (char *)&iv[16], (ftnlen)sizeof(integer));
    e_wsfe();
L250:
    i__ = 1;
    j = jlim[alg1 - 1];
    k = 19;
    l = iv[49];
    ndfalt = ndflt[alg1 - 1];
    i__1 = ndfalt;
    for (ii = 1; ii <= i__1; ++ii) {
	if (v[k] == v[l]) {
	    goto L280;
	}
	if (m == 0) {
	    io___186.ciunit = pu;
	    s_wsfe(&io___186);
	    do_fio(&c__3, which, (ftnlen)4);
	    e_wsfe();
	}
	m = 1;
	io___187.ciunit = pu;
	s_wsfe(&io___187);
	do_fio(&c__1, vn + ((i__ << 1) - 2 << 2), (ftnlen)4);
	do_fio(&c__1, vn + ((i__ << 1) - 1 << 2), (ftnlen)4);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&v[k], (ftnlen)sizeof(doublereal));
	e_wsfe();
L280:
	++k;
	++l;
	++i__;
	if (i__ == j) {
	    i__ = ijmp;
	}
/* L290: */
    }

    iv[54] = iv[16];
    parsv1 = iv[49];
    dv7cpy_(&iv[50], &v[parsv1], &v[19]);
    goto L999;

L300:
    iv[1] = 15;
    if (pu == 0) {
	goto L999;
    }
    io___189.ciunit = pu;
    s_wsfe(&io___189);
    do_fio(&c__1, (char *)&(*liv), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&miv2, (ftnlen)sizeof(integer));
    e_wsfe();
    if (*liv < miv1) {
	goto L999;
    }
    if (*lv < iv[45]) {
	goto L320;
    }
    goto L999;

L320:
    iv[1] = 16;
    if (pu != 0) {
	io___190.ciunit = pu;
	s_wsfe(&io___190);
	do_fio(&c__1, (char *)&(*lv), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iv[45], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L999;

L340:
    iv[1] = 67;
    if (pu != 0) {
	io___191.ciunit = pu;
	s_wsfe(&io___191);
	do_fio(&c__1, (char *)&(*alg), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L999;
L360:
    if (pu != 0) {
	io___192.ciunit = pu;
	s_wsfe(&io___192);
	do_fio(&c__1, (char *)&(*liv), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&miv1, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (44 <= *liv) {
	iv[44] = miv1;
    }
    if (45 <= *liv) {
	iv[45] = 0;
    }

L999:
    return 0;
/*  ***  LAST LINE OF DPARCK FOLLOWS  *** */
} /* dparck_ */

/* Subroutine */ int dq7rsh_(integer *k, integer *p, logical *havqtr, 
	doublereal *qtr, doublereal *r__, doublereal *w)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a, b;
    static integer i__, j;
    static doublereal t, x, y, z__;
    static integer i1, j1, k1;
    static doublereal wj;
    static integer jm1, km1, jp1, pm1;
    extern /* Subroutine */ int dh2rfa_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    extern doublereal dh2rfg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern /* Subroutine */ int dv7cpy_(integer *, doublereal *, doublereal *)
	    ;


/*  ***  PERMUTE COLUMN K OF R TO COLUMN P, MODIFY QTR ACCORDINGLY  *** */

/*     DIMSNSION R(P*(P+1)/2) */


/*  ***  LOCAL VARIABLES  *** */


    /* Parameter adjustments */
    --w;
    --qtr;
    --r__;

    /* Function Body */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    if (*k >= *p) {
	goto L999;
    }
    km1 = *k - 1;
    k1 = *k * km1 / 2;
    dv7cpy_(k, &w[1], &r__[k1 + 1]);
    wj = w[*k];
    pm1 = *p - 1;
    j1 = k1 + km1;
    i__1 = pm1;
    for (j = *k; j <= i__1; ++j) {
	jm1 = j - 1;
	jp1 = j + 1;
	if (jm1 > 0) {
	    dv7cpy_(&jm1, &r__[k1 + 1], &r__[j1 + 2]);
	}
	j1 += jp1;
	k1 += j;
	a = r__[j1];
	b = r__[j1 + 1];
	if (b != zero) {
	    goto L10;
	}
	r__[k1] = a;
	x = zero;
	z__ = zero;
	goto L40;
L10:
	r__[k1] = dh2rfg_(&a, &b, &x, &y, &z__);
	if (j == pm1) {
	    goto L30;
	}
	i1 = j1;
	i__2 = pm1;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    i1 += i__;
	    dh2rfa_(&c__1, &r__[i1], &r__[i1 + 1], &x, &y, &z__);
/* L20: */
	}
L30:
	if (*havqtr) {
	    dh2rfa_(&c__1, &qtr[j], &qtr[jp1], &x, &y, &z__);
	}
L40:
	t = x * wj;
	w[j] = wj + t;
	wj = t * z__;
/* L50: */
    }
    w[*p] = wj;
    dv7cpy_(p, &r__[k1 + 1], &w[1]);
L999:
    return 0;
} /* dq7rsh_ */

doublereal dr7mdc_(integer *k)
{
    /* Initialized data */

    static doublereal big = 0.;
    static doublereal eta = 0.;
    static doublereal machep = 0.;
    static doublereal zero = 0.;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern doublereal d1mach_(integer *);


/*  ***  RETURN MACHINE DEPENDENT CONSTANTS USED BY NL2SOL  *** */


/*  ***  THE CONSTANT RETURNED DEPENDS ON K... */

/*  ***        K = 1... SMALLEST POS. ETA SUCH THAT -ETA EXISTS. */
/*  ***        K = 2... SQUARE ROOT OF ETA. */
/*  ***        K = 3... UNIT ROUNDOFF = SMALLEST POS. NO. MACHEP SUCH */
/*  ***                 THAT 1 + MACHEP .GT. 1 .AND. 1 - MACHEP .LT. 1. */
/*  ***        K = 4... SQUARE ROOT OF MACHEP. */
/*  ***        K = 5... SQUARE ROOT OF BIG (SEE K = 6). */
/*  ***        K = 6... LARGEST MACHINE NO. BIG SUCH THAT -BIG EXISTS. */

/* /+ */
/* / */

    if (big > zero) {
	goto L1;
    }
    big = d1mach_(&c__2);
    eta = d1mach_(&c__1);
    machep = d1mach_(&c__4);
L1:

/* -------------------------------  BODY  -------------------------------- */

    switch (*k) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	case 5:  goto L50;
	case 6:  goto L60;
    }

L10:
    ret_val = eta;
    goto L999;

L20:
    ret_val = sqrt(eta * 256.) / 16.;
    goto L999;

L30:
    ret_val = machep;
    goto L999;

L40:
    ret_val = sqrt(machep);
    goto L999;

L50:
    ret_val = sqrt(big / 256.) * 16.;
    goto L999;

L60:
    ret_val = big;

L999:
    return ret_val;
/*  ***  LAST CARD OF DR7MDC FOLLOWS  *** */
} /* dr7mdc_ */

doublereal drldst_(integer *p, doublereal *d__, doublereal *x, doublereal *x0)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal t, emax, xmax;


/*  ***  COMPUTE AND RETURN RELATIVE DIFFERENCE BETWEEN X AND X0  *** */
/*  ***  NL2SOL VERSION 2.2  *** */


/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

/*  ***  BODY  *** */

    /* Parameter adjustments */
    --x0;
    --x;
    --d__;

    /* Function Body */
    emax = 0.;
    xmax = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = (d__1 = d__[i__] * (x[i__] - x0[i__]), abs(d__1));
	if (emax < t) {
	    emax = t;
	}
	t = d__[i__] * ((d__1 = x[i__], abs(d__1)) + (d__2 = x0[i__], abs(
		d__2)));
	if (xmax < t) {
	    xmax = t;
	}
/* L10: */
    }
    ret_val = 0.;
    if (xmax > 0.) {
	ret_val = emax / xmax;
    }
/* L999: */
    return ret_val;
/*  ***  LAST CARD OF DRLDST FOLLOWS  *** */
} /* drldst_ */

/* Subroutine */ int drmngb_(doublereal *b, doublereal *d__, doublereal *fx, 
	doublereal *g, integer *iv, integer *liv, integer *lv, integer *n, 
	doublereal *v, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t;
    static integer z__, i1, n1, w1, g01;
    static doublereal gi;
    static integer x01;
    static doublereal xi;
    static integer dg1, td1, tg1, np1, ipi, ipn, temp0, temp1, step1, dummy;
    extern /* Subroutine */ int dd7dgb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern logical stopx_(integer *);
    extern /* Subroutine */ int dl7upd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *);
    static integer dstep1;
    extern /* Subroutine */ int dw7zbf_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int da7sst_(integer *, integer *, integer *, 
	    doublereal *), i7shft_(integer *, integer *, integer *), dl7vml_(
	    integer *, doublereal *, doublereal *, doublereal *);
    extern doublereal dv2nrm_(integer *, doublereal *);
    extern /* Subroutine */ int dq7rsh_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *), dv7scp_(integer *, 
	    doublereal *, doublereal *), dv7ipr_(integer *, integer *, 
	    doublereal *), dv7cpy_(integer *, doublereal *, doublereal *), 
	    dl7tvm_(integer *, doublereal *, doublereal *, doublereal *), 
	    dv2axy_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dv7vmp_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static integer nwtst1;
    extern /* Subroutine */ int dparck_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *);
    extern doublereal drldst_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int divset_(integer *, integer *, integer *, 
	    integer *, doublereal *), ditsum_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *);
    static integer lstgst, rstrst;


/*  ***  CARRY OUT  DMNGB (SIMPLY BOUNDED MINIMIZATION) ITERATIONS, */
/*  ***  USING DOUBLE-DOGLEG/BFGS STEPS. */

/*  ***  PARAMETER DECLARATIONS  *** */


/* --------------------------  PARAMETER USAGE  -------------------------- */

/* B.... VECTOR OF LOWER AND UPPER BOUNDS ON X. */
/* D.... SCALE VECTOR. */
/* FX... FUNCTION VALUE. */
/* G.... GRADIENT VECTOR. */
/* IV... INTEGER VALUE ARRAY. */
/* LIV.. LENGTH OF IV (AT LEAST 59) + N. */
/* LV... LENGTH OF V (AT LEAST 71 + N*(N+19)/2). */
/* N.... NUMBER OF VARIABLES (COMPONENTS IN X AND G). */
/* V.... FLOATING-POINT VALUE ARRAY. */
/* X.... VECTOR OF PARAMETERS TO BE OPTIMIZED. */

/*  ***  DISCUSSION  *** */

/*        PARAMETERS IV, N, V, AND X ARE THE SAME AS THE CORRESPONDING */
/*     ONES TO  DMNGB (WHICH SEE), EXCEPT THAT V CAN BE SHORTER (SINCE */
/*     THE PART OF V THAT  DMNGB USES FOR STORING G IS NOT NEEDED). */
/*     MOREOVER, COMPARED WITH  DMNGB, IV(1) MAY HAVE THE TWO ADDITIONAL */
/*     OUTPUT VALUES 1 AND 2, WHICH ARE EXPLAINED BELOW, AS IS THE USE */
/*     OF IV(TOOBIG) AND IV(NFGCAL).  THE VALUE IV(G), WHICH IS AN */
/*     OUTPUT VALUE FROM  DMNGB (AND SMSNOB), IS NOT REFERENCED BY */
/*     DRMNGB OR THE SUBROUTINES IT CALLS. */
/*        FX AND G NEED NOT HAVE BEEN INITIALIZED WHEN DRMNGB IS CALLED */
/*     WITH IV(1) = 12, 13, OR 14. */

/* IV(1) = 1 MEANS THE CALLER SHOULD SET FX TO F(X), THE FUNCTION VALUE */
/*             AT X, AND CALL DRMNGB AGAIN, HAVING CHANGED NONE OF THE */
/*             OTHER PARAMETERS.  AN EXCEPTION OCCURS IF F(X) CANNOT BE */
/*             (E.G. IF OVERFLOW WOULD OCCUR), WHICH MAY HAPPEN BECAUSE */
/*             OF AN OVERSIZED STEP.  IN THIS CASE THE CALLER SHOULD SET */
/*             IV(TOOBIG) = IV(2) TO 1, WHICH WILL CAUSE DRMNGB TO IG- */
/*             NORE FX AND TRY A SMALLER STEP.  THE PARAMETER NF THAT */
/*              DMNGB PASSES TO CALCF (FOR POSSIBLE USE BY CALCG) IS A */
/*             COPY OF IV(NFCALL) = IV(6). */
/* IV(1) = 2 MEANS THE CALLER SHOULD SET G TO G(X), THE GRADIENT VECTOR */
/*             OF F AT X, AND CALL DRMNGB AGAIN, HAVING CHANGED NONE OF */
/*             THE OTHER PARAMETERS EXCEPT POSSIBLY THE SCALE VECTOR D */
/*             WHEN IV(DTYPE) = 0.  THE PARAMETER NF THAT  DMNGB PASSES */
/*             TO CALCG IS IV(NFGCAL) = IV(7).  IF G(X) CANNOT BE */
/*             EVALUATED, THEN THE CALLER MAY SET IV(NFGCAL) TO 0, IN */
/*             WHICH CASE DRMNGB WILL RETURN WITH IV(1) = 65. */
/* . */
/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (DECEMBER 1979).  REVISED SEPT. 1982. */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED */
/*     IN PART BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324 AND MCS-7906671. */

/*        (SEE  DMNG FOR REFERENCES.) */

/* +++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */


/*     ***  CONSTANTS  *** */


/*  ***  NO INTRINSIC FUNCTIONS  *** */

/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* DA7SST.... ASSESSES CANDIDATE STEP. */
/* DD7DGB... COMPUTES SIMPLY BOUNDED DOUBLE-DOGLEG (CANDIDATE) STEP. */
/* DIVSET.... SUPPLIES DEFAULT IV AND V INPUT COMPONENTS. */
/* DD7TPR... RETURNS INNER PRODUCT OF TWO VECTORS. */
/* I7SHFT... CYCLICALLLY SHIFTS AN ARRAY OF INTEGERS. */
/* DITSUM.... PRINTS ITERATION SUMMARY AND INFO ON INITIAL AND FINAL X. */
/* DL7TVM... MULTIPLIES TRANSPOSE OF LOWER TRIANGLE TIMES VECTOR. */
/* LUPDT.... UPDATES CHOLESKY FACTOR OF HESSIAN APPROXIMATION. */
/* DL7VML.... MULTIPLIES LOWER TRIANGLE TIMES VECTOR. */
/* DPARCK.... CHECKS VALIDITY OF INPUT IV AND V VALUES. */
/* DQ7RSH... CYCLICALLY SHIFTS CHOLESKY FACTOR. */
/* DRLDST... COMPUTES V(RELDX) = RELATIVE STEP SIZE. */
/* STOPX.... RETURNS .TRUE. IF THE BREAK KEY HAS BEEN PRESSED. */
/* DV2NRM... RETURNS THE 2-NORM OF A VECTOR. */
/* DV2AXY.... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER. */
/* DV7CPY.... COPIES ONE VECTOR TO ANOTHER. */
/* DV7IPR... CYCLICALLY SHIFTS A FLOATING-POINT ARRAY. */
/* DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */
/* DV7VMP... MULTIPLIES VECTOR BY VECTOR RAISED TO POWER (COMPONENTWISE). */
/* DW7ZBF... COMPUTES W AND Z FOR DL7UPD CORRESPONDING TO BFGS UPDATE. */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/*  ***  (NOTE THAT NC IS STORED IN IV(G0)) *** */

/* /6 */
/*     DATA CNVCOD/55/, DG/37/, INITH/25/, IRC/29/, IVNEED/3/, KAGQT/33/, */
/*    1     MODE/35/, MODEL/5/, MXFCAL/17/, MXITER/18/, NC/48/, */
/*    2     NEXTIV/46/, NEXTV/47/, NFCALL/6/, NFGCAL/7/, NGCALL/30/, */
/*    3     NITER/31/, NWTSTP/34/, PERM/58/, RADINC/8/, RESTOR/9/, */
/*    4     STEP/40/, STGLIM/11/, STLSTG/41/, TOOBIG/2/, XIRC/13/, X0/43/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DGNORM/1/, DINIT/38/, DSTNRM/2/, F/10/, F0/13/, FDIF/11/, */
/*    1     GTSTEP/4/, INCFAC/23/, LMAT/42/, LMAX0/35/, LMAXS/36/, */
/*    2     PREDUC/7/, RADFAC/16/, RADIUS/8/, RAD0/9/, RELDX/17/, */
/*    3     TUNER4/29/, TUNER5/30/, VNEED/4/ */
/* /7 */
/* / */

/* /6 */
/*     DATA NEGONE/-1.D+0/, ONE/1.D+0/, ONEP2/1.2D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --v;
    --x;
    --g;
    --d__;
    b -= 3;

    /* Function Body */
    i__ = iv[1];
    if (i__ == 1) {
	goto L70;
    }
    if (i__ == 2) {
	goto L80;
    }

/*  ***  CHECK VALIDITY OF IV AND V INPUT VALUES  *** */

    if (iv[1] == 0) {
	divset_(&c__2, &iv[1], liv, lv, &v[1]);
    }
    if (iv[1] < 12) {
	goto L10;
    }
    if (iv[1] > 13) {
	goto L10;
    }
    iv[4] += *n * (*n + 19) / 2;
    iv[3] += *n;
L10:
    dparck_(&c__2, &d__[1], &iv[1], liv, lv, n, &v[1]);
    i__ = iv[1] - 2;
    if (i__ > 12) {
	goto L999;
    }
    switch (i__) {
	case 1:  goto L250;
	case 2:  goto L250;
	case 3:  goto L250;
	case 4:  goto L250;
	case 5:  goto L250;
	case 6:  goto L250;
	case 7:  goto L190;
	case 8:  goto L150;
	case 9:  goto L190;
	case 10:  goto L20;
	case 11:  goto L20;
	case 12:  goto L30;
    }

/*  ***  STORAGE ALLOCATION  *** */

L20:
    l = iv[42];
    iv[43] = l + *n * (*n + 1) / 2;
    iv[40] = iv[43] + (*n << 1);
    iv[41] = iv[40] + (*n << 1);
    iv[34] = iv[41] + *n;
    iv[37] = iv[34] + (*n << 1);
    iv[47] = iv[37] + (*n << 1);
    iv[46] = iv[58] + *n;
    if (iv[1] != 13) {
	goto L30;
    }
    iv[1] = 14;
    goto L999;

/*  ***  INITIALIZATION  *** */

L30:
    iv[31] = 0;
    iv[6] = 1;
    iv[30] = 1;
    iv[7] = 1;
    iv[35] = -1;
    iv[5] = 1;
    iv[11] = 1;
    iv[2] = 0;
    iv[55] = 0;
    iv[8] = 0;
    iv[48] = *n;
    v[9] = 0.;

/*  ***  CHECK CONSISTENCY OF B AND INITIALIZE IP ARRAY  *** */

    ipi = iv[58];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iv[ipi] = i__;
	++ipi;
	if (b[(i__ << 1) + 1] > b[(i__ << 1) + 2]) {
	    goto L410;
	}
/* L40: */
    }

    if (v[38] >= 0.) {
	dv7scp_(n, &d__[1], &v[38]);
    }
    if (iv[25] != 1) {
	goto L60;
    }

/*     ***  SET THE INITIAL HESSIAN APPROXIMATION TO DIAG(D)**-2  *** */

    l = iv[42];
    i__1 = *n * (*n + 1) / 2;
    dv7scp_(&i__1, &v[l], &c_b34);
    k = l - 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k += i__;
	t = d__[i__];
	if (t <= 0.) {
	    t = 1.;
	}
	v[k] = t;
/* L50: */
    }

/*  ***  GET INITIAL FUNCTION VALUE  *** */

L60:
    iv[1] = 1;
    goto L440;

L70:
    v[10] = *fx;
    if (iv[35] >= 0) {
	goto L250;
    }
    v[13] = *fx;
    iv[1] = 2;
    if (iv[2] == 0) {
	goto L999;
    }
    iv[1] = 63;
    goto L430;

/*  ***  MAKE SURE GRADIENT COULD BE COMPUTED  *** */

L80:
    if (iv[2] == 0) {
	goto L90;
    }
    iv[1] = 65;
    goto L430;

/*  ***  CHOOSE INITIAL PERMUTATION  *** */

L90:
    ipi = iv[58];
    ipn = ipi + *n;
    n1 = *n;
    np1 = *n + 1;
    l = iv[42];
    w1 = iv[34] + *n;
    k = *n - iv[48];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--ipn;
	j = iv[ipn];
	if (b[(j << 1) + 1] >= b[(j << 1) + 2]) {
	    goto L100;
	}
	xi = x[j];
	gi = g[j];
	if (xi <= b[(j << 1) + 1] && gi > 0.) {
	    goto L100;
	}
	if (xi >= b[(j << 1) + 2] && gi < 0.) {
	    goto L100;
	}
/*           *** DISALLOW CONVERGENCE IF X(J) HAS JUST BEEN FREED *** */
	if (i__ <= k) {
	    iv[55] = 0;
	}
	goto L120;
L100:
	i1 = np1 - i__;
	if (i1 >= n1) {
	    goto L110;
	}
	i7shft_(&n1, &i1, &iv[ipi]);
	dq7rsh_(&i1, &n1, &c_false, &g[1], &v[l], &v[w1]);
L110:
	--n1;
L120:
	;
    }

    iv[48] = n1;
    v[1] = 0.;
    if (n1 <= 0) {
	goto L130;
    }
    dg1 = iv[37];
    dv7vmp_(n, &v[dg1], &g[1], &d__[1], &c_n1);
    dv7ipr_(n, &iv[ipi], &v[dg1]);
    v[1] = dv2nrm_(&n1, &v[dg1]);
L130:
    if (iv[55] != 0) {
	goto L420;
    }
    if (iv[35] == 0) {
	goto L370;
    }

/*  ***  ALLOW FIRST STEP TO HAVE SCALED 2-NORM AT MOST V(LMAX0)  *** */

    v[8] = v[35];

    iv[35] = 0;


/* -----------------------------  MAIN LOOP  ----------------------------- */


/*  ***  PRINT ITERATION SUMMARY, CHECK ITERATION LIMIT  *** */

L140:
    ditsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);
L150:
    k = iv[31];
    if (k < iv[18]) {
	goto L160;
    }
    iv[1] = 10;
    goto L430;

/*  ***  UPDATE RADIUS  *** */

L160:
    iv[31] = k + 1;
    if (k == 0) {
	goto L170;
    }
    t = v[16] * v[2];
    if (v[16] < 1. || t > v[8]) {
	v[8] = t;
    }

/*  ***  INITIALIZE FOR START OF NEXT ITERATION  *** */

L170:
    x01 = iv[43];
    v[13] = v[10];
    iv[29] = 4;
    iv[33] = -1;

/*     ***  COPY X TO X0  *** */

    dv7cpy_(n, &v[x01], &x[1]);

/*  ***  CHECK STOPX AND FUNCTION EVALUATION LIMIT  *** */

L180:
    if (! stopx_(&dummy)) {
	goto L200;
    }
    iv[1] = 11;
    goto L210;

/*     ***  COME HERE WHEN RESTARTING AFTER FUNC. EVAL. LIMIT OR STOPX. */

L190:
    if (v[10] >= v[13]) {
	goto L200;
    }
    v[16] = 1.;
    k = iv[31];
    goto L160;

L200:
    if (iv[6] < iv[17]) {
	goto L220;
    }
    iv[1] = 9;
L210:
    if (v[10] >= v[13]) {
	goto L430;
    }

/*        ***  IN CASE OF STOPX OR FUNCTION EVALUATION LIMIT WITH */
/*        ***  IMPROVED V(F), EVALUATE THE GRADIENT AT X. */

    iv[55] = iv[1];
    goto L360;

/* . . . . . . . . . . . . .  COMPUTE CANDIDATE STEP  . . . . . . . . . . */

L220:
    step1 = iv[40];
    dg1 = iv[37];
    nwtst1 = iv[34];
    w1 = nwtst1 + *n;
    dstep1 = step1 + *n;
    ipi = iv[58];
    l = iv[42];
    tg1 = dg1 + *n;
    x01 = iv[43];
    td1 = x01 + *n;
    dd7dgb_(&b[3], &d__[1], &v[dg1], &v[dstep1], &g[1], &iv[ipi], &iv[33], &v[
	    l], lv, n, &iv[48], &v[nwtst1], &v[step1], &v[td1], &v[tg1], &v[1]
	    , &v[w1], &v[x01]);
    if (iv[29] != 6) {
	goto L230;
    }
    if (iv[9] != 2) {
	goto L250;
    }
    rstrst = 2;
    goto L260;

/*  ***  CHECK WHETHER EVALUATING F(X0 + STEP) LOOKS WORTHWHILE  *** */

L230:
    iv[2] = 0;
    if (v[2] <= 0.) {
	goto L250;
    }
    if (iv[29] != 5) {
	goto L240;
    }
    if (v[16] <= 1.) {
	goto L240;
    }
    if (v[7] > v[11] * 1.2) {
	goto L240;
    }
    if (iv[9] != 2) {
	goto L250;
    }
    rstrst = 0;
    goto L260;

/*  ***  COMPUTE F(X0 + STEP)  *** */

L240:
    dv2axy_(n, &x[1], &c_b52, &v[step1], &v[x01]);
    ++iv[6];
    iv[1] = 1;
    goto L440;

/* . . . . . . . . . . . . .  ASSESS CANDIDATE STEP  . . . . . . . . . . . */

L250:
    rstrst = 3;
L260:
    x01 = iv[43];
    v[17] = drldst_(n, &d__[1], &x[1], &v[x01]);
    da7sst_(&iv[1], liv, lv, &v[1]);
    step1 = iv[40];
    lstgst = iv[41];
    i__ = iv[9] + 1;
    switch (i__) {
	case 1:  goto L300;
	case 2:  goto L270;
	case 3:  goto L280;
	case 4:  goto L290;
    }
L270:
    dv7cpy_(n, &x[1], &v[x01]);
    goto L300;
L280:
    dv7cpy_(n, &v[lstgst], &x[1]);
    goto L300;
L290:
    dv7cpy_(n, &x[1], &v[lstgst]);
    dv2axy_(n, &v[step1], &c_b403, &v[x01], &x[1]);
    v[17] = drldst_(n, &d__[1], &x[1], &v[x01]);
    iv[9] = rstrst;

L300:
    k = iv[29];
    switch (k) {
	case 1:  goto L310;
	case 2:  goto L340;
	case 3:  goto L340;
	case 4:  goto L340;
	case 5:  goto L310;
	case 6:  goto L320;
	case 7:  goto L330;
	case 8:  goto L330;
	case 9:  goto L330;
	case 10:  goto L330;
	case 11:  goto L330;
	case 12:  goto L330;
	case 13:  goto L400;
	case 14:  goto L370;
    }

/*     ***  RECOMPUTE STEP WITH CHANGED RADIUS  *** */

L310:
    v[8] = v[16] * v[2];
    goto L180;

/*  ***  COMPUTE STEP OF LENGTH V(LMAXS) FOR SINGULAR CONVERGENCE TEST. */

L320:
    v[8] = v[36];
    goto L220;

/*  ***  CONVERGENCE OR FALSE CONVERGENCE  *** */

L330:
    iv[55] = k - 4;
    if (v[10] >= v[13]) {
	goto L420;
    }
    if (iv[13] == 14) {
	goto L420;
    }
    iv[13] = 14;

/* . . . . . . . . . . . .  PROCESS ACCEPTABLE STEP  . . . . . . . . . . . */

L340:
    x01 = iv[43];
    step1 = iv[40];
    dv2axy_(n, &v[step1], &c_b403, &v[x01], &x[1]);
    if (iv[29] != 3) {
	goto L360;
    }

/*     ***  SET  TEMP1 = HESSIAN * STEP  FOR USE IN GRADIENT TESTS  *** */

/*     ***  USE X0 AS TEMPORARY... */

    ipi = iv[58];
    dv7cpy_(n, &v[x01], &v[step1]);
    dv7ipr_(n, &iv[ipi], &v[x01]);
    l = iv[42];
    dl7tvm_(n, &v[x01], &v[l], &v[x01]);
    dl7vml_(n, &v[x01], &v[l], &v[x01]);

/*        *** UNPERMUTE X0 INTO TEMP1 *** */

    temp1 = iv[41];
    temp0 = temp1 - 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iv[ipi];
	++ipi;
	k = temp0 + j;
	v[k] = v[x01];
	++x01;
/* L350: */
    }

/*  ***  SAVE OLD GRADIENT, COMPUTE NEW ONE  *** */

L360:
    g01 = iv[34] + *n;
    dv7cpy_(n, &v[g01], &g[1]);
    ++iv[30];
    iv[2] = 0;
    iv[1] = 2;
    goto L999;

/*  ***  INITIALIZATIONS -- G0 = G - G0, ETC.  *** */

L370:
    g01 = iv[34] + *n;
    dv2axy_(n, &v[g01], &c_b403, &v[g01], &g[1]);
    step1 = iv[40];
    temp1 = iv[41];
    if (iv[29] != 3) {
	goto L390;
    }

/*  ***  SET V(RADFAC) BY GRADIENT TESTS  *** */

/*     ***  SET  TEMP1 = DIAG(D)**-1 * (HESSIAN*STEP + (G(X0)-G(X)))  *** */

    dv2axy_(n, &v[temp1], &c_b403, &v[g01], &v[temp1]);
    dv7vmp_(n, &v[temp1], &v[temp1], &d__[1], &c_n1);

/*        ***  DO GRADIENT TESTS  *** */

    if (dv2nrm_(n, &v[temp1]) <= v[1] * v[29]) {
	goto L380;
    }
    if (dd7tpr_(n, &g[1], &v[step1]) >= v[4] * v[30]) {
	goto L390;
    }
L380:
    v[16] = v[23];

/*  ***  UPDATE H, LOOP  *** */

L390:
    w1 = iv[34];
    z__ = iv[43];
    l = iv[42];
    ipi = iv[58];
    dv7ipr_(n, &iv[ipi], &v[step1]);
    dv7ipr_(n, &iv[ipi], &v[g01]);
    dw7zbf_(&v[l], n, &v[step1], &v[w1], &v[g01], &v[z__]);

/*     ** USE THE N-VECTORS STARTING AT V(STEP1) AND V(G01) FOR SCRATCH.. */
    dl7upd_(&v[temp1], &v[step1], &v[l], &v[g01], &v[l], n, &v[w1], &v[z__]);
    iv[1] = 2;
    goto L140;

/* . . . . . . . . . . . . . .  MISC. DETAILS  . . . . . . . . . . . . . . */

/*  ***  BAD PARAMETERS TO ASSESS  *** */

L400:
    iv[1] = 64;
    goto L430;

/*  ***  INCONSISTENT B  *** */

L410:
    iv[1] = 82;
    goto L430;

/*  ***  PRINT SUMMARY OF FINAL ITERATION AND OTHER REQUESTED ITEMS  *** */

L420:
    iv[1] = iv[55];
    iv[55] = 0;
L430:
    ditsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);
    goto L999;

/*  ***  PROJECT X INTO FEASIBLE REGION (PRIOR TO COMPUTING F OR G)  *** */

L440:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] < b[(i__ << 1) + 1]) {
	    x[i__] = b[(i__ << 1) + 1];
	}
	if (x[i__] > b[(i__ << 1) + 2]) {
	    x[i__] = b[(i__ << 1) + 2];
	}
/* L450: */
    }

L999:
    return 0;

/*  ***  LAST CARD OF DRMNGB FOLLOWS  *** */
} /* drmngb_ */

/* Subroutine */ int dv2axy_(integer *p, doublereal *w, doublereal *a, 
	doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET W = A*X + Y  --  W, X, Y = P-VECTORS, A = SCALAR  *** */



    /* Parameter adjustments */
    --y;
    --x;
    --w;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	w[i__] = *a * x[i__] + y[i__];
    }
    return 0;
} /* dv2axy_ */

doublereal dv2nrm_(integer *p, doublereal *x)
{
    /* Initialized data */

    static doublereal sqteta = 0.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal r__, t, xi, scale;
    extern doublereal dr7mdc_(integer *);


/*  ***  RETURN THE 2-NORM OF THE P-VECTOR X, TAKING  *** */
/*  ***  CARE TO AVOID THE MOST LIKELY UNDERFLOWS.    *** */


/* /+ */
/* / */

/* /6 */
/*     DATA ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --x;

    /* Function Body */

    if (*p > 0) {
	goto L10;
    }
    ret_val = 0.;
    goto L999;
L10:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] != 0.) {
	    goto L30;
	}
/* L20: */
    }
    ret_val = 0.;
    goto L999;

L30:
    scale = (d__1 = x[i__], abs(d__1));
    if (i__ < *p) {
	goto L40;
    }
    ret_val = scale;
    goto L999;
L40:
    t = 1.;
    if (sqteta == 0.) {
	sqteta = dr7mdc_(&c__2);
    }

/*     ***  SQTETA IS (SLIGHTLY LARGER THAN) THE SQUARE ROOT OF THE */
/*     ***  SMALLEST POSITIVE FLOATING POINT NUMBER ON THE MACHINE. */
/*     ***  THE TESTS INVOLVING SQTETA ARE DONE TO PREVENT UNDERFLOWS. */

    j = i__ + 1;
    i__1 = *p;
    for (i__ = j; i__ <= i__1; ++i__) {
	xi = (d__1 = x[i__], abs(d__1));
	if (xi > scale) {
	    goto L50;
	}
	r__ = xi / scale;
	if (r__ > sqteta) {
	    t += r__ * r__;
	}
	goto L60;
L50:
	r__ = scale / xi;
	if (r__ <= sqteta) {
	    r__ = 0.;
	}
	t = t * r__ * r__ + 1.;
	scale = xi;
L60:
	;
    }

    ret_val = scale * sqrt(t);
L999:
    return ret_val;
/*  ***  LAST LINE OF DV2NRM FOLLOWS  *** */
} /* dv2nrm_ */

/* Subroutine */ int dv7cpy_(integer *p, doublereal *y, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET Y = X, WHERE X AND Y ARE P-VECTORS  *** */



    /* Parameter adjustments */
    --x;
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = x[i__];
    }
    return 0;
} /* dv7cpy_ */

/* Subroutine */ int dv7dfl_(integer *alg, integer *lv, doublereal *v)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal dr7mdc_(integer *);
    static doublereal machep, mepcrt, sqteps;


/*  ***  SUPPLY ***SOL (VERSION 2.3) DEFAULT VALUES TO V  *** */

/*  ***  ALG = 1 MEANS REGRESSION CONSTANTS. */
/*  ***  ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS. */


/* DR7MDC... RETURNS MACHINE-DEPENDENT CONSTANTS */


/*  ***  SUBSCRIPTS FOR V  *** */


/* /6 */
/*     DATA ONE/1.D+0/, THREE/3.D+0/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA AFCTOL/31/, BIAS/43/, COSMIN/47/, DECFAC/22/, DELTA0/44/, */
/*    1     DFAC/41/, DINIT/38/, DLTFDC/42/, DLTFDJ/43/, DTINIT/39/, */
/*    2     D0INIT/40/, EPSLON/19/, ETA0/42/, FUZZ/45/, HUBERC/48/, */
/*    3     INCFAC/23/, LMAX0/35/, LMAXS/36/, PHMNFC/20/, PHMXFC/21/, */
/*    4     RDFCMN/24/, RDFCMX/25/, RFCTOL/32/, RLIMIT/46/, RSPTOL/49/, */
/*    5     SCTOL/37/, SIGMIN/50/, TUNER1/26/, TUNER2/27/, TUNER3/28/, */
/*    6     TUNER4/29/, TUNER5/30/, XCTOL/33/, XFTOL/34/ */
/* /7 */
/* / */

/* -------------------------------  BODY  -------------------------------- */

    /* Parameter adjustments */
    --v;

    /* Function Body */
    machep = dr7mdc_(&c__3);
    v[31] = 1e-20;
    if (machep > 1e-10) {
/* Computing 2nd power */
	d__1 = machep;
	v[31] = d__1 * d__1;
    }
    v[22] = .5;
    sqteps = dr7mdc_(&c__4);
    v[41] = .6;
    v[39] = 1e-6;
    mepcrt = pow_dd(&machep, &c_b433);
    v[40] = 1.;
    v[19] = .1;
    v[23] = 2.;
    v[35] = 1.;
    v[36] = 1.;
    v[20] = -.1;
    v[21] = .1;
    v[24] = .1;
    v[25] = 4.;
/* Computing MAX */
/* Computing 2nd power */
    d__3 = mepcrt;
    d__1 = 1e-10, d__2 = d__3 * d__3;
    v[32] = max(d__1,d__2);
    v[37] = v[32];
    v[26] = .1;
    v[27] = 1e-4;
    v[28] = .75;
    v[29] = .5;
    v[30] = .75;
    v[33] = sqteps;
    v[34] = machep * 100.;

    if (*alg >= 2) {
	goto L10;
    }

/*  ***  REGRESSION  VALUES */

/* Computing MAX */
    d__1 = 1e-6, d__2 = machep * 100.;
    v[47] = max(d__1,d__2);
    v[38] = 0.;
    v[44] = sqteps;
    v[42] = mepcrt;
    v[43] = sqteps;
    v[45] = 1.5;
    v[48] = .7;
    v[46] = dr7mdc_(&c__5);
    v[49] = .001;
    v[50] = 1e-4;
    goto L999;

/*  ***  GENERAL OPTIMIZATION VALUES */

L10:
    v[43] = .8;
    v[38] = -1.;
    v[42] = machep * 1e3;

L999:
    return 0;
/*  ***  LAST CARD OF DV7DFL FOLLOWS  *** */
} /* dv7dfl_ */

/* Subroutine */ int dv7ipr_(integer *n, integer *ip, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;


/*     PERMUTE X SO THAT X.OUTPUT(I) = X.INPUT(IP(I)). */
/*     IP IS UNCHANGED ON OUTPUT. */


    /* Parameter adjustments */
    --x;
    --ip;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ip[i__];
	if (j == i__) {
	    goto L30;
	}
	if (j > 0) {
	    goto L10;
	}
	ip[i__] = -j;
	goto L30;
L10:
	t = x[i__];
	k = i__;
L20:
	x[k] = x[j];
	k = j;
	j = ip[k];
	ip[k] = -j;
	if (j > i__) {
	    goto L20;
	}
	x[k] = t;
L30:
	;
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DV7IPR FOLLOWS  *** */
} /* dv7ipr_ */

/* Subroutine */ int dv7scp_(integer *p, doublereal *y, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET P-VECTOR Y TO SCALAR S  *** */



    /* Parameter adjustments */
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = *s;
    }
    return 0;
} /* dv7scp_ */

/* Subroutine */ int dv7shf_(integer *n, integer *k, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal t;
    static integer nm1;


/*  ***  SHIFT X(K),...,X(N) LEFT CIRCULARLY ONE POSITION  *** */



    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*k >= *n) {
	goto L999;
    }
    nm1 = *n - 1;
    t = x[*k];
    i__1 = nm1;
    for (i__ = *k; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = x[i__ + 1];
    }
    x[*n] = t;
L999:
    return 0;
} /* dv7shf_ */

/* Subroutine */ int dv7vmp_(integer *n, doublereal *x, doublereal *y, 
	doublereal *z__, integer *k)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* ***  SET X(I) = Y(I) * Z(I)**K, 1 .LE. I .LE. N (FOR K = 1 OR -1)  *** */


    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    if (*k >= 0) {
	goto L20;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = y[i__] / z__[i__];
    }
    goto L999;

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	x[i__] = y[i__] * z__[i__];
    }
L999:
    return 0;
/*  ***  LAST CARD OF DV7VMP FOLLOWS  *** */
} /* dv7vmp_ */

/* Subroutine */ int dw7zbf_(doublereal *l, integer *n, doublereal *s, 
	doublereal *w, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal cs, cy, ys, shs, theta, epsrt;
    extern /* Subroutine */ int dl7ivm_(integer *, doublereal *, doublereal *,
	     doublereal *);
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int dl7tvm_(integer *, doublereal *, doublereal *,
	     doublereal *);


/*  ***  COMPUTE  Y  AND  Z  FOR  DL7UPD  CORRESPONDING TO BFGS UPDATE. */

/*     DIMENSION L(N*(N+1)/2) */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/* L (I/O) CHOLESKY FACTOR OF HESSIAN, A LOWER TRIANG. MATRIX STORED */
/*             COMPACTLY BY ROWS. */
/* N (INPUT) ORDER OF  L  AND LENGTH OF  S,  W,  Y,  Z. */
/* S (INPUT) THE STEP JUST TAKEN. */
/* W (OUTPUT) RIGHT SINGULAR VECTOR OF RANK 1 CORRECTION TO L. */
/* Y (INPUT) CHANGE IN GRADIENTS CORRESPONDING TO S. */
/* Z (OUTPUT) LEFT SINGULAR VECTOR OF RANK 1 CORRECTION TO L. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  ALGORITHM NOTES  *** */

/*        WHEN  S  IS COMPUTED IN CERTAIN WAYS, E.G. BY  GQTSTP  OR */
/*     DBLDOG,  IT IS POSSIBLE TO SAVE N**2/2 OPERATIONS SINCE  (L**T)*S */
/*     OR  L*(L**T)*S IS THEN KNOWN. */
/*        IF THE BFGS UPDATE TO L*(L**T) WOULD REDUCE ITS DETERMINANT TO */
/*     LESS THAN EPS TIMES ITS OLD VALUE, THEN THIS ROUTINE IN EFFECT */
/*     REPLACES  Y  BY  THETA*Y + (1 - THETA)*L*(L**T)*S,  WHERE  THETA */
/*     (BETWEEN 0 AND 1) IS CHOSEN TO MAKE THE REDUCTION FACTOR = EPS. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (FALL 1979). */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED */
/*     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND */
/*     MCS-7906671. */

/* ------------------------  EXTERNAL QUANTITIES  ------------------------ */

/*  ***  FUNCTIONS AND SUBROUTINES CALLED  *** */

/* DD7TPR RETURNS INNER PRODUCT OF TWO VECTORS. */
/* DL7IVM MULTIPLIES L**-1 TIMES A VECTOR. */
/* DL7TVM MULTIPLIES L**T TIMES A VECTOR. */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/* --------------------------  LOCAL VARIABLES  -------------------------- */


/*  ***  DATA INITIALIZATIONS  *** */

/* /6 */
/*     DATA EPS/0.1D+0/, ONE/1.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --l;
    --z__;
    --y;
    --w;
    --s;

    /* Function Body */
    dl7tvm_(n, &w[1], &l[1], &s[1]);
    shs = dd7tpr_(n, &w[1], &w[1]);
    ys = dd7tpr_(n, &y[1], &s[1]);
    if (ys >= shs * .1) {
	goto L10;
    }
    theta = shs * .90000000000000002 / (shs - ys);
    epsrt = sqrt(.1);
    cy = theta / (shs * epsrt);
    cs = ((theta - 1.) / epsrt + 1.) / shs;
    goto L20;
L10:
    cy = 1. / (sqrt(ys) * sqrt(shs));
    cs = 1. / shs;
L20:
    dl7ivm_(n, &z__[1], &l[1], &y[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	z__[i__] = cy * z__[i__] - cs * w[i__];
    }

/* L999: */
    return 0;
/*  ***  LAST CARD OF DW7ZBF FOLLOWS  *** */
} /* dw7zbf_ */

integer i7mdcn_(integer *k)
{
    /* Initialized data */

    static integer mdperm[3] = { 2,4,1 };

    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer i1mach_(integer *);



/*  ***  RETURN INTEGER MACHINE-DEPENDENT CONSTANTS  *** */

/*     ***  K = 1 MEANS RETURN STANDARD OUTPUT UNIT NUMBER.   *** */
/*     ***  K = 2 MEANS RETURN ALTERNATE OUTPUT UNIT NUMBER.  *** */
/*     ***  K = 3 MEANS RETURN  INPUT UNIT NUMBER.            *** */
/*          (NOTE -- K = 2, 3 ARE USED ONLY BY TEST PROGRAMS.) */

/*  +++  PORT VERSION FOLLOWS... */
    ret_val = i1mach_(&mdperm[(0 + (0 + (*k - 1 << 2))) / 4]);
/*  +++  END OF PORT VERSION  +++ */

/*  +++  NON-PORT VERSION FOLLOWS... */
/*     INTEGER MDCON(3) */
/*     DATA MDCON(1)/6/, MDCON(2)/8/, MDCON(3)/5/ */
/*     I7MDCN = MDCON(K) */
/*  +++  END OF NON-PORT VERSION  +++ */

/* L999: */
    return ret_val;
/*  ***  LAST CARD OF I7MDCN FOLLOWS  *** */
} /* i7mdcn_ */

/* Subroutine */ int i7shft_(integer *n, integer *k, integer *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, t, k1, ii, nm1;


/*  ***  SHIFT X(K),...,X(N) LEFT CIRCULARLY ONE POSITION IF K .GT. 0. */
/*  ***  SHIFT X(-K),...,X(N) RIGHT CIRCULARLY ONE POSITION IF K .LT. 0. */



    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*k < 0) {
	goto L20;
    }
    if (*k >= *n) {
	goto L999;
    }
    nm1 = *n - 1;
    t = x[*k];
    i__1 = nm1;
    for (i__ = *k; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = x[i__ + 1];
    }
    x[*n] = t;
    goto L999;

L20:
    k1 = -(*k);
    if (k1 >= *n) {
	goto L999;
    }
    t = x[*n];
    nm1 = *n - k1;
    i__1 = nm1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n - ii;
	x[i__ + 1] = x[i__];
/* L30: */
    }
    x[k1] = t;
L999:
    return 0;
/*  ***  LAST LINE OF I7SHFT FOLLOWS  *** */
} /* i7shft_ */

logical stopx_(integer *idummy)
{
    /* System generated locals */
    logical ret_val;

/*     *****PARAMETERS... */

/*     .................................................................. */

/*     *****PURPOSE... */
/*     THIS FUNCTION MAY SERVE AS THE STOPX (ASYNCHRONOUS INTERRUPTION) */
/*     FUNCTION FOR THE NL2SOL (NONLINEAR LEAST-SQUARES) PACKAGE AT */
/*     THOSE INSTALLATIONS WHICH DO NOT WISH TO IMPLEMENT A */
/*     DYNAMIC STOPX. */

/*     *****ALGORITHM NOTES... */
/*     AT INSTALLATIONS WHERE THE NL2SOL SYSTEM IS USED */
/*     INTERACTIVELY, THIS DUMMY STOPX SHOULD BE REPLACED BY A */
/*     FUNCTION THAT RETURNS .TRUE. IF AND ONLY IF THE INTERRUPT */
/*     (BREAK) KEY HAS BEEN PRESSED SINCE THE LAST CALL ON STOPX. */

/*     .................................................................. */

    ret_val = FALSE_;
    return ret_val;
} /* stopx_ */

#ifdef __cplusplus
	}
#endif
