/* dpptri.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c/f2c.h"

/* Table of constant values */

static doublereal c_b8 = 1.;
static integer c__1 = 1;

/* Subroutine */ int dpptri_(char *uplo, integer *n, doublereal *ap, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, jc, jj;
    static doublereal ajj;
    static integer jjn;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dtptri_(
	    char *, char *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);


/*  -- LAPACK routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     March 31, 1993 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DPPTRI computes the inverse of a real symmetric positive definite */
/*  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T */
/*  computed by DPPTRF. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  Upper triangular factor is stored in AP; */
/*          = 'L':  Lower triangular factor is stored in AP. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/*          On entry, the triangular factor U or L from the Cholesky */
/*          factorization A = U**T*U or A = L*L**T, packed columnwise as */
/*          a linear array.  The j-th column of U or L is stored in the */
/*          array AP as follows: */
/*          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; */
/*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */

/*          On exit, the upper or lower triangle of the (symmetric) */
/*          inverse of A, overwriting the input factor U or L. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, the (i,i) element of the factor U or L is */
/*                zero, and the inverse could not be computed. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --ap;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPPTRI", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Invert the triangular Cholesky factor U or L. */

    dtptri_(uplo, "Non-unit", n, &ap[1], info, (ftnlen)1, (ftnlen)8);
    if (*info > 0) {
	return 0;
    }

    if (upper) {

/*        Compute the product inv(U) * inv(U)'. */

	jj = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    jc = jj + 1;
	    jj += j;
	    if (j > 1) {
		i__2 = j - 1;
		dspr_("Upper", &i__2, &c_b8, &ap[jc], &c__1, &ap[1], (ftnlen)
			5);
	    }
	    ajj = ap[jj];
	    dscal_(&j, &ajj, &ap[jc], &c__1);
/* L10: */
	}

    } else {

/*        Compute the product inv(L)' * inv(L). */

	jj = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    jjn = jj + *n - j + 1;
	    i__2 = *n - j + 1;
	    ap[jj] = ddot_(&i__2, &ap[jj], &c__1, &ap[jj], &c__1);
	    if (j < *n) {
		i__2 = *n - j;
		dtpmv_("Lower", "Transpose", "Non-unit", &i__2, &ap[jjn], &ap[
			jj + 1], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
	    }
	    jj = jjn;
/* L20: */
	}
    }

    return 0;

/*     End of DPPTRI */

} /* dpptri_ */

/* Subroutine */ int dtptri_(char *uplo, char *diag, integer *n, doublereal *
	ap, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, jc, jj;
    static doublereal ajj;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer jclast;
    static logical nounit;


/*  -- LAPACK routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTPTRI computes the inverse of a real upper or lower triangular */
/*  matrix A stored in packed format. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  A is upper triangular; */
/*          = 'L':  A is lower triangular. */

/*  DIAG    (input) CHARACTER*1 */
/*          = 'N':  A is non-unit triangular; */
/*          = 'U':  A is unit triangular. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/*          On entry, the upper or lower triangular matrix A, stored */
/*          columnwise in a linear array.  The j-th column of A is stored */
/*          in the array AP as follows: */
/*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/*          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n. */
/*          See below for further details. */
/*          On exit, the (triangular) inverse of the original matrix, in */
/*          the same packed storage format. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular */
/*                matrix is singular and its inverse can not be computed. */

/*  Further Details */
/*  =============== */

/*  A triangular matrix A can be transferred to packed storage using one */
/*  of the following program segments: */

/*  UPLO = 'U':                      UPLO = 'L': */

/*        JC = 1                           JC = 1 */
/*        DO 2 J = 1, N                    DO 2 J = 1, N */
/*           DO 1 I = 1, J                    DO 1 I = J, N */
/*              AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J) */
/*      1    CONTINUE                    1    CONTINUE */
/*           JC = JC + J                      JC = JC + N - J + 1 */
/*      2 CONTINUE                       2 CONTINUE */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --ap;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTPTRI", &i__1, (ftnlen)6);
	return 0;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
	if (upper) {
	    jj = 0;
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		jj += *info;
		if (ap[jj] == 0.) {
		    return 0;
		}
/* L10: */
	    }
	} else {
	    jj = 1;
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		if (ap[jj] == 0.) {
		    return 0;
		}
		jj = jj + *n - *info + 1;
/* L20: */
	    }
	}
	*info = 0;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

	jc = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (nounit) {
		ap[jc + j - 1] = 1. / ap[jc + j - 1];
		ajj = -ap[jc + j - 1];
	    } else {
		ajj = -1.;
	    }

/*           Compute elements 1:j-1 of j-th column. */

	    i__2 = j - 1;
	    dtpmv_("Upper", "No transpose", diag, &i__2, &ap[1], &ap[jc], &
		    c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
	    i__2 = j - 1;
	    dscal_(&i__2, &ajj, &ap[jc], &c__1);
	    jc += j;
/* L30: */
	}

    } else {

/*        Compute inverse of lower triangular matrix. */

	jc = *n * (*n + 1) / 2;
	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		ap[jc] = 1. / ap[jc];
		ajj = -ap[jc];
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

		i__1 = *n - j;
		dtpmv_("Lower", "No transpose", diag, &i__1, &ap[jclast], &ap[
			jc + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
		i__1 = *n - j;
		dscal_(&i__1, &ajj, &ap[jc + 1], &c__1);
	    }
	    jclast = jc;
	    jc = jc - *n + j - 2;
/* L40: */
	}
    }

    return 0;

/*     End of DTPTRI */

} /* dtptri_ */

logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static integer inta, intb, zcode;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/*  case. */

/*  Arguments */
/*  ========= */

/*  CA      (input) CHARACTER*1 */
/*  CB      (input) CHARACTER*1 */
/*          CA and CB specify the single characters to be compared. */

/* ===================================================================== */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test if the characters are equal */

    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower or */
/*        upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or */
/*        upper case 'Z'. */

	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code */
/*        plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

    return ret_val;
} /* lsame_ */

/* Subroutine */ int xerbla_(char *srname, integer *info, ftnlen srname_len)
{
    /* Format strings */
    static char fmt_9999[] = "(\002 ** On entry to \002,a6,\002 parameter nu\
mber \002,i2,\002 had \002,\002an illegal value\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___17 = { 0, 6, 0, fmt_9999, 0 };



/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  XERBLA  is an error handler for the LAPACK routines. */
/*  It is called by an LAPACK routine if an input parameter has an */
/*  invalid value.  A message is printed and execution stops. */

/*  Installers may consider modifying the STOP statement in order to */
/*  call system-specific exception-handling facilities. */

/*  Arguments */
/*  ========= */

/*  SRNAME  (input) CHARACTER*6 */
/*          The name of the routine which called XERBLA. */

/*  INFO    (input) INTEGER */
/*          The position of the invalid parameter in the parameter list */
/*          of the calling routine. */

/* ===================================================================== */

/*     .. Executable Statements .. */

    s_wsfe(&io___17);
    do_fio(&c__1, srname, (ftnlen)6);
    do_fio(&c__1, (char *)&(*info), (ftnlen)sizeof(integer));
    e_wsfe();

    s_stop("", (ftnlen)0);


/*     End of XERBLA */

    return 0;
} /* xerbla_ */

#ifdef __cplusplus
	}
#endif
