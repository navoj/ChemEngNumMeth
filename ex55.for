C
C*****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM SOLVES THE BVP POSED IN EX 5.3 USING FINITE DIFFER-
C  ENCE APPROXIMATIONS WITH AN DIRECT SOLUTION PROCEDURE.
C
C***************************  NOMENCLATURE  *****************************
C
C  B(I)-  THE CONSTANT TERM IN THE ITH EQUATION
C  C(I)-  THE COEFFICIENT OF THE TERM TO THE LEFT TO THE ITH DIAGONAL
C         ELEMENT
C  COND-  THE THERMAL CONDUCTIVITY OF TIN (KJ/M-SEC-DEG C)
C  D(I)-  THE ITH DIAGONAL ELEMENT
C  DX-  THE SEPARATION BETWEEN NODE POINTS (M)
C  E(I)-  THE COEFFICIENT OF THE TERM TO THE RIGHT OF THE ITH DIAGONAL
C         ELEMENT
C  H-  THE HEAT TRANSFER COEFFICIENT (KJ/M**2-SEC-DEG C)
C  N-  THE NUMBER OF INDEPENDENT EQUATIONS
C  TB-  THE AMBIENT AIR TEMPERATURE (DEG C)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(50),D(50),E(50),B(50),X(50),BETA(50),GAM(50)
C  SET THE PARAMETERS OF THE PROBLEM
      N=19
      DX=.05/FLOAT(N+1)
      H=4180.
      TB=25.
      COND=.602
C  FORM THE NON-ZERO COEFICIENTS OF THE LINEAR EQUATIONS FOR THE
C  INTERIOR NODES
      NM=N-1
      DO 1 I=2,NM
      C(I)=1.
      D(I)=-2.
      E(I)=1.
    1 B(I)=-59340.*DX*DX
C  CALCULATE THE COEFICIENTS OF THE NODES ADJACENT TO THE BOUNDARY NODES
      Y=H+3.*COND/2./DX
      D(1)=2.*COND/DX/Y-2.
      E(1)=1.-COND/2./DX/Y
      B(1)=-59340.*DX*DX-H*TB/Y
      C(N)=2./3.
      D(N)=-2./3.
      B(N)=-59340.*DX*DX
C  CALL THOMAS METHOD FOR SOLUTION
      CALL TM(N,C,D,E,B,X,BETA,GAM)
      X(N+1)=(4.*X(N)-X(N-1))/3.
      NP=N+1
C  PRINT OUT RESULTS
      WRITE(6,90)
   90 FORMAT( 29H       X(M)       TEMP(DEG C))
      DO 100 I=1,NP
      Z=DX*FLOAT(I)
  100 WRITE(6,110)Z,X(I)
  110 FORMAT( 6X,F6.4,7X,F7.3)
      STOP
      END
C
C************************   ABSTRACT   *****************************
C
C      THIS SUBROUTINE CALCULATES THE SOLUTION OF A SYSTEM OF LINEAR
C   EQUATIONS WHICH HAVE A TRIDIAGONAL COEFICIENT MATRIX USING THE
C   THE THOMAS METHOD
C
C************************************************************************
C
      SUBROUTINE TM(N,C,D,E,B,X,BETA,GAM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(1),D(1),E(1),B(1),X(1),BETA(1),GAM(1)
      BETA(1)=D(1)
      GAM(1)=B(1)/BETA(1)
      DO 10 I=2,N
      BETA(I)=D(I)-C(I)*E(I-1)/BETA(I-1)
   10 GAM(I)=(B(I)-C(I)*GAM(I-1))/BETA(I)
      X(N)=GAM(N)
      DO 20 I=2,N
      J=N-I+1
   20 X(J)=GAM(J)-E(J)*X(J+1)/BETA(J)
      RETURN
      END
