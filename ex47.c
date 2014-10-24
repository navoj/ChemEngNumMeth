/* ex47.f -- translated by f2c (version 19950110).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/**************************  ABSTRACT  ************************************/

/*   THIS PROGRAM USES THE EXPLICIT EULER METHOD TO INTEGRATE A SYSTEM OF */
/*  FIRST ORDER ODE'S.  THE USER SUPPLIES THE INITIAL CONDITIONS AND THE */
/*  FINAL VALUE OF X, AS WELL AS THE EQUATIONS FOR F(X,Y). */

/* ****************************  NOMENCLATURE  ************************** */

/*  DX-  THE STEP SIZE */
/*  DXPNT- THE PRINT INTERVAL */
/*  DYDX(I)- THE DERIVATIVE OF Y(I) WITH RESPECT TO X */
/*  FUNCT- THE NAME OF THE SUBROUTINE THAT CALCULATES F(X,Y) */
/*  IPNT- PRINT FLAG (=0 NO PRINT; =1 PRINT) */
/*  N- THE NUMBER OF DEPENDENT VARIABLES */
/*  X-  THE INDEPENDENT VARIABLE */
/*  XMAX-  THE FINAL VALUE OF X */
/*  Y(I)-  THE ITH DEPENDENT VARIABLE */

/* ********************************************************************** */

/* Main program */ MAIN__()
{
    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    static doublereal dydx[2];
    static integer ipnt;
    static doublereal xmax;
    static integer n;
    static doublereal x, y[2];
    extern /* Subroutine */ int euler_();
    extern /* Subroutine */ int funct_();
    static doublereal dxpnt, dx;

/*  EXTERNAL STATEMENT OF DERIVATIVE SUBROUTINE */
/*  SPECIFY INPUT DATA */
    n = 2;
    x = (float)0.;
    y[0] = (float)1.;
    y[1] = (float)300.;
    xmax = (float)100.;
    dx = (float).02;
    ipnt = 1;
    dxpnt = (float)10.;
/*  CALL EULER METHOD */
    euler_(funct_, &n, &dx, &x, y, &xmax, &ipnt, &dxpnt, dydx);
    s_stop("", 0L);
} /* MAIN__ */


/* *****************************  ABSTRACT  ***************************** */

/*     THIS SUBROUTINE CALLS SUBROUTINE E1 AT EACH PRINT INTERVAL */
/*  AND PRINTS OUT THE INTERMEDIATE RESULTS IF IPNT=1. */

/* ****************************  NOMENCLATURE  ************************** */

/*  DX-  THE STEP SIZE */
/*  DXPNT- THE PRINT INTERVAL */
/*  DXP- THE INTERNALLY CALCULATED PRINT INTERVAL */
/*  DYDX(I)- THE DERIVATIVE OF Y(I) WITH RESPECT TO X */
/*  F- THE NAME OF THE SUBROUTINE THAT CALCULATES F(X,Y) */
/*  IPNT- PRINT FLAG (=0 NO PRINT; =1 PRINT) */
/*  N- THE NUMBER OF DEPENDENT VARIABLES */
/*  NS- THE NUMBER OF PRINT INTERVAL STEPS */
/*  PRO- THE PRINT INTERVAL NEAR XMAX */
/*  X- THE INDEPENDENT VARIABLE */
/*  XF- THE VALUE OF X AT THE END OF THE PRINT INTERVAL */
/*  XO- THE VALUE OF X AT THE BEGINNING OF THE PRINT INTERVAL */
/*  XMAX-  THE FINAL VALUE OF X */
/*  Y(I)-  THE ITH DEPENDENT VARIABLE */

/* ********************************************************************** */

/* Subroutine */ int euler_(f, n, dx, x, y, xmax, ipnt, dxpnt, dydx)
/* Subroutine */ int (*f) ();
integer *n;
doublereal *dx, *x, *y, *xmax;
integer *ipnt;
doublereal *dxpnt, *dydx;
{
    /* Format strings */
    static char fmt_88[] = "(\002 DX HAS BEEN SPECIFIED TO BE MORE THAN DX\
PNT\002)";
    static char fmt_10[] = "(10x,\002 X=\002,d12.5,3x,\002 Y=\002,4(d14.7,3x\
))";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe();
    /* Subroutine */ int s_stop();
    integer do_fio();

    /* Local variables */
    extern /* Subroutine */ int e_();
    static integer i, j;
    static doublereal xf;
    static integer ns;
    static doublereal xo, dxp, pro;

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 6, 0, fmt_88, 0 };
    static cilist io___10 = { 0, 6, 0, fmt_10, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_10, 0 };


/*  CHECK DXPNT */
    /* Parameter adjustments */
    --dydx;
    --y;

    /* Function Body */
    if (*dx > *dxpnt) {
	s_wsfe(&io___9);
	e_wsfe();
    }
    if (*dx > *dxpnt) {
	s_stop("", 0L);
    }
/*  PRINT INITIAL CONDITIONS IF IPNT = 1 */
    if (*ipnt == 1) {
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&y[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
/*  CALCULATE THE NUMBER OF STEPS PER DXPNT */
    dxp = *dxpnt;
    ns = (integer) (*xmax / *dxpnt) + 1;
    xo = *x;
    pro = *xmax - *dxpnt * (doublereal) (ns - 1);

/*  CALL SUBROUTINE E1 TIL X=XMAX */

    i__1 = ns;
    for (i = 1; i <= i__1; ++i) {
	if (i == ns) {
	    dxp = pro;
	}
	if (i == ns && pro < (float)1e-6) {
	    goto L1000;
	}
	xf = xo + dxp;
/*  CALL SUBROUTINE E1 */
	e_(f, n, dx, &xo, &xf, &y[1], &dydx[1]);
/*  PRINT INTERMEDIATE RESULTS IF IPNT=1 */
	if (*ipnt == 1) {
	    s_wsfe(&io___18);
	    do_fio(&c__1, (char *)&xf, (ftnlen)sizeof(doublereal));
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&y[j], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
	}
/*  INCREMENT X */
L1000:
	xo += dxp;
    }
    return 0;
} /* euler_ */


/* *****************************  ABSTRACT  ***************************** */

/*     THIS SUBROUTINE APPLIES THE EULER METHOD FOR THE INTEGRATION OF */
/*  A SYSTEM OF FIRST ORDER ODE'S FROM X=XO TO X=XF. */

/* ****************************  NOMENCLATURE  ************************* */

/*  DX-  THE STEP SIZE */
/*  DXX- INTERNALLY USED STEP SIZE */
/*  DYDX(I)- THE DERIVATIVE OF Y(I) WITH RESPECT TO X */
/*  F- THE NAME OF THE SUBROUTINE THAT CALCULATES F(X,Y) */
/*  N- THE NUMBER OF DEPENDENT VARIABLES */
/*  NS- THE NUMBER OF DX INTERVAL STEPS */
/*  PRO- THE VALUE OF DX NEAR XF */
/*  X- THE INDEPENDENT VARIABLE */
/*  XF- THE VALUE OF X AT THE END OF THE PRINT INTERVAL */
/*  XO- THE VALUE OF X AT THE BEGINNING OF THE PRINT INTERVAL */
/*  XMAX-  THE FINAL VALUE OF X */
/*  Y(I)-  THE ITH DEPENDENT VARIABLE */

/**************************************************************************/


/* Subroutine */ int e_(f, n, dx, xo, xf, y, dydx)
/* Subroutine */ int (*f) ();
integer *n;
doublereal *dx, *xo, *xf, *y, *dydx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, k;
    static doublereal x, xd;
    static integer ns;
    static doublereal pro, dxx;

    /* Parameter adjustments */
    --dydx;
    --y;

    /* Function Body */
    x = *xo;
/*  CALCULATE THE NUMBER OF STEPS TO REACH XF */
    xd = *xf - *xo;
    ns = (integer) (xd / *dx) + 1;
    dxx = *dx;
    pro = xd - *dx * (doublereal) (ns - 1);

/*  APPLY EULER APPROXIMATION FROM X=XO TO X=XF */

    i__1 = ns;
    for (i = 1; i <= i__1; ++i) {
/*  CALCULATE DX FOR END OF PRINT INTERVAL */
	if (i == ns) {
	    dxx = pro;
	}
/*  CALCULATE DYDX(I) */
	(*f)(n, &x, &y[1], &dydx[1]);
	if (*xf - x < dxx) {
	    dxx = *xf - x;
	}
	x += dxx;
/*  APPLY EULER APPROXIMATION */
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
/* L10: */
	    y[k] += dxx * dydx[k];
	}
/* L100: */
    }
    return 0;
} /* e_ */


/* *****************************  ABSTRACT  ***************************** */

/*     THIS SUBROUTINE CALCULATES THE DERIVATIVE OF Y WITH RESPECT TO X, */
/*  DYDX.  THIS EQUATION IS USER SUPPLIED. */

/* ********************************************************************** */
/* Subroutine */ int funct_(n, x, y, dydx)
integer *n;
doublereal *x, *y, *dydx;
{
    /* Builtin functions */
    double exp();

    /* Parameter adjustments */
    --dydx;
    --y;

    /* Function Body */
    dydx[1] = y[1] * (float)-.1 * exp((float)-300. / y[2]);
    dydx[2] = y[1] * (float)1. * exp((float)-300. / y[2]);
    return 0;
} /* funct_ */

