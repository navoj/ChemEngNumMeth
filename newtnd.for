C
C****************************  ABSTRACT   *******************************
C
C      THIS SUBROUTINE EMPLOYES NEWTON'S METHOD IN ORDER TO SOLVE A
C  SET OF N NONLINEAR EQUATIONS CONTAINING N UNKNOWNS.
C  THIS SUBROUTINE IS CALLED BY THE MAIN PROGRAM AND IS SUPPLIED
C  THE VALUES OF THE INITIAL GUESS FOR X(I)'S AS WELL AS THE VALUE OF
C  N.  THIS SUBROUTINE USES THE VALUES OF THE FUNCTION FROM FUNC AND
C  THE VALUES OF THE PARTIAL DERIVATIVES OF THE FUNCTION IN ORDER TO
C  DETERMINE THE SOLUTION. THIS METHOD USES THE LIBRARY ROUTINE LINPAC
C  TO SOLVE THE SYSTEM OF LINEAR EQUATION USED BY NEWTON'S METHOD.
C
C************************************************************************
C
C
      SUBROUTINE NEWTND(FUNC,DER,N,X,FX,FX1,ERLIM,A,B,XX,IPVT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,1),B(1),XX(1),IPVT(1),X(1),FX(1),FX1(1)
      EXTERNAL FUNC,DER
    1 CONTINUE
      ITEST=0
C
C  MAKE FUNCTION EVALUATIONS
C
      CALL FUNC(N,X,FX)
      DO 3 I=1,N
    3 B(I)=-FX(I)
C
C  EVALUATE JACOBIAN MATRIX
C
      CALL DER(FUNC,N,X,A,FX,FX1)
C
C  CALL LINEAR EQUATION SOLVER
C
      CALL LINPAC(N,A,B,XX,IPVT)
C
C  MAKE AN IMPROVED VALUE FOR X(I)
C
      ITEST=0
      DO 5 I=1,N
C  ANTICIPATE ZERO ROOTS
      TS=DABS(X(I))
      IF(TS.LT.1.D-10)GO TO 5
      RAT=XX(I)/X(I)
      IF(DABS(RAT).GT.ERLIM)ITEST=ITEST+1
    5 X(I)=X(I)+XX(I)
C
C  CHECK FOR CONVERGENCE
C
      IF(ITEST.NE.0)GO TO 1
      RETURN
      END
C
C*************************   ABSTRACT  **********************************
C
C     THIS SUBROUTINE CALCULATES THE PARTIAL DERIVATIVES OF THE
C  THE FUNCTIONS WITH RESPECT TO THE INDEPENDENT VARIABLES.
C  A(I,J) REPRESENTS THE PARTIAL OF THE ITH FUNCTION WITH RESPECT
C  TO THE JTH VARIABLE.
C
C************************************************************************
C
C
      SUBROUTINE DER(FUNC,N,X,A,FXB,FXD)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,1),X(1),FXB(1),FXD(1)
      REAL K1,K2,K3,K4
      COMMON /ONE/K1,K2,K3,K4,VR,Q,CAO
      EXTERNAL FUNC
      DELTA=.01
C  CALCULATE FUNCTION VALUE AT X(I)
      CALL FUNC(N,X,FXB)
C
C  CALCULATE NUMERICAL DERIVATIVES USING FINITE DIFFERENCE  EQUATIONS
C
      DO 1 I=1,N
C  INCREMENT X(I) BY DELTA
      X(I)=X(I)*(1.+DELTA)
C  CONSIDER CASES FOR VERY SMALL VALUES OF X(I)
      IF(DABS(X(I)).LT.1.D-10)XS=X(I)
      IF(DABS(X(I)).LT.1.D-10)X(I)=1.D-3
C  CALL FUNCTION VALUE AT INCREMENTED VALUE OF X(I)
      CALL FUNC(N,X,FXD)
      DEL=DELTA*X(I)/(1.+DELTA)
C  SET DEL FOR CASE IN WHICH X(I) IS SMALL
      IF(X(I).EQ.1.D-3)DEL=1.D-3
      IF(X(I).EQ.1.D-3)X(I)=XS
C  CORRECT X(I) BACK TO ORIGINAL VALUE
      X(I)=X(I)/(1.+DELTA)
C  CALCULATE DERIVATIVES WITH RESPECT TO X(I)
      DO 2 J=1,N
    2 A(J,I)=(FXD(J)-FXB(J))/DEL
    1 CONTINUE
      RETURN
      END
