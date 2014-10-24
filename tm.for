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
