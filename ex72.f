
C
C************************  ABSTRACT  *******************************
C
C     THIS PROGRAM USES LINEAR REGRESSION TO DETERMINE THE RATE CONSTANT
C  FROM CONCENTRATION VERSUS TIME DATA FOR A FIRST ORDER REACTION.
C
C***********************  NOMENCLATURE  ****************************
C
C  CA(I)-  CONCENTRATION OF THE REACTANT (GMOLES/SEC)
C  N-  THE NUMBER OF DATA POINTS
C  SEI-  THE STANDARD ERROR OF THE Y-INTERCEPT
C  SELR-  THE STANDARD ERROR OF THE LINEAR REGRESSION
C  SES-  THE STANDARD ERROR OF THE SLOPE
C  R-  THE CORRELATION COEFICIENT
C  RSELR- THE NORMALIZED STANDARD ERROR OF THE LINEAR REGRESSION
C  TIME(I)-  THE VALUE OF THE TIME FOR THE ITH DATA POINT  (SEC)
C  X(I)-  THE ITH VALUE OF THE INDEPENDENT VARIABLE(DATA)
C  Y(I)-  THE ITH VALUE OF THE DEPENDENT VARIABLE (DATA)
C  YI- THE Y-INTERCEPT OF THE LINEAR REGRESSION
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(100),Y(100)
      DIMENSION CA(20),TIME(20)
C
C  INPUT DATA CONCENTRATION VERSUS TIME DATA
C
      CA(1)=45.01
      CA(2)=45.81
      CA(3)=46.65
      CA(4)=47.52
      CA(5)=48.44
      CA(6)=49.39
      TIME(1)=320.0
      TIME(2)=330.0
      TIME(3)=340.0
      TIME(4)=350.0
      TIME(5)=360.0
      TIME(6)=370.0
      N=6
C
C  FORM (X,Y) DATA FROM LINEARIZATION OF CONCENTRATION/TIME DATA
C
      DO 10 I=1,N
      Y(I)=CA(I)
   10 X(I)=TIME(I)
C
C  CALL LINEAR REGRESSION SUBROUTINE
C
      CALL LNRGS(N,X,Y,S,YI,SES,SEYI,SELR,R,RSELR)
C
C  PRINT OUT RATE CONSTANT AND ITS UNCERTAINTY
C
      WRITE(6,20)S
   20 FORMAT( 34H THE FIRST ORDER RATE CONSTANT IS ,D14.7,7HSECONDS)
      WRITE(6,25)SES
   25 FORMAT( 43H THE STANDARD ERROR OF THE RATE CONSTANT IS,D14.7)
      STOP
      END
C
C******************************  ABSTRACT  ******************************
C
C     THIS SUBROUTINE PERFORMS A LINEAR REGRESSION ON A SET OF (X/Y)
C  DATA.  IN ADDITION TO DETERMINING THE SLOPE AND INTERCEPT OF THE
C  DATA, THE STANDARD ERROR OF THE REGRESSION, THE SLOPE, AND THE
C  INTERCEPT ARE DETERMINED.  ALSO, THE CORRELATION COEFICIENT AND THE
C  NORMALIZED ERROR OF THE REGRESSION ARE CALCULATED.
C
C****************************  NOMENCLATURE  ****************************
C
C  R-  THE CORRELATION COEFICIENT
C  RSELR-  THE NORMALIZED STANDARD ERROR OF THE LINEAR REGRESSION
C  S-  THE SLOPE OF THE LINEAR REGRESSION
C  SB-  THE SUM OF THE SQUARES  BETWEEN THE DATA AND THE AVERAGE Y
C  SELR-  THE STANDARD ERROR OF THE LINERAR REGRESSION
C  SES-  THE STANDARD ERROR OF THE SLOPE
C  SEYI- THE STANDARD ERROR OF THE INTERCEPT
C  S2-  THE SUM OF THE SQUARES OF THE ERROR
C  X(I)-  THE ITH VALUE OF THE DEPENDENT VARIABLE
C  XAVG-  THE AVERAGE VALUE OF X
C  XT-  THE SUM OF ALL X(I)
C  XY-  THE SUM OF ALL X(I)*Y(I)
C  X2-  THE SUM OF ALL X(I)*X(I)
C  Y(I)-  THE ITH VALUE OF THE DEPENDENT VARIABLE
C  YAVG-  THE AVERAGE VALUE OF Y
C  YI-  THE Y-INTERCEPT OF THE REGRESSION LINE
C  YT-  THE SUM OF ALL Y(I)
C  Y2-  THE SUM OF ALL Y(I)*Y(I)
C
C************************************************************************
C
C
      SUBROUTINE LNRGS(N,X,Y,S,YI,SES,SEYI,SELR,R,RSELR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(1),Y(1)
C
C  DETERMINE XT, YT, X2, XY, AND Y2
C
      XT=0.0
      YT=0.0
      X2=0.0
      XY=0.0
      Y2=0.0
      DO 10 I=1,N
      XT=XT+X(I)
      YT=YT+Y(I)
      X2=X2+X(I)*X(I)
      Y2=Y2+Y(I)*Y(I)
   10 XY=XY+X(I)*Y(I)
C
C  FROM XT, YT, XY, AND X2 CALCULATE THE SLOPE AND Y-INTERCEPT
C
      XN=FLOAT(N)
      S=(XN*XY-XT*YT)/(XN*X2-XT*XT)
      XAVG=XT/XN
      YAVG=YT/XN
      YI=YAVG-S*XAVG
C
C  CALCULATE S2 AND SB
C
      S2=0.0
      SB=0.0
      DO 100 I=1,N
      SB=SB+(Y(I)-YAVG)**2
  100 S2=S2+(YI+S*X(I)-Y(I))**2
C
C  CALCULATE SELR AND THE CORRELATION FACTOR
C
      SELR=DSQRT(S2/(XN-2))
      R=DSQRT((SB-S2)/SB)
C
C  CALCULATE SES AND SEYI
C
      XD=0.0
      DO 200 I=1,N
  200 XD=XD+(X(I)-XAVG)**2
      SES=DSQRT(S2/((XN-2.)*XD))
      SEYI=DSQRT(S2*X2/(XN*(XN-2)*XD))
C
C  CALCULATE RSELR
C
      YMAX=Y(1)
      YMIN=Y(1)
      DO 400 I=2,N
      IF(Y(I).LT.YMIN)YMIN=Y(I)
  400 IF(Y(I).GT.YMAX)YMAX=Y(I)
      RSELR=SELR/DABS(YMAX-YMIN)
C
C  PRINT OUT RESULTS
C
      WRITE(6,19)
   19 FORMAT( ///)
      WRITE(6,20)YI,S
   20 FORMAT( 19H THE Y-INTERCEPT IS,D14.7,5X,12H THE SLOPE =,D14.7)
      WRITE(6,21)SELR
   21 FORMAT( 45H THE STANDARD ERROR OF THE REGRESSION LINE = ,D14.7)
      WRITE(6,22)R
   22 FORMAT( 30H THE CORRELATION COEFICIENT = ,D14.7)
      WRITE(6,23)SES
   23 FORMAT( 35H THE STANDARD ERROR OF THE SLOPE = ,D14.7)
      WRITE(6,24)SEYI
   24 FORMAT( 40H THE STANDARD ERROR OF THE Y-INTERCEPT =,D14.7)
      WRITE(6,25)RSELR
   25 FORMAT( 45H NORMALIZED STANDARD ERROR OF THE REGRESSION=,D14.7)
      WRITE(6,19)
      RETURN
      END
