C
C*****************************  ABSTRACT  *******************************
C
C     THIS SUBROUTINE CALCULATES THE SECOND DERIVATIVES AT EACH NODE
C  POINT.
C
C***************************  NOMENCLATURE  *****************************
      SUBROUTINE SECDER(ND,X,FD,FPP,C,D,E,B,FX,BETA,GAM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(1),FD(1),FPP(1),C(1),D(1),E(1),B(1),FX(1)
      DIMENSION BETA(1),GAM(1)
C  FORM LINEAR EQUATION THAT DETERMINE THE SECOND DERIVATIVES
      N=ND-2
      DO 1 I=1,N
      C(I)=X(I+1)-X(I)
      D(I)=2.*(X(I+2)-X(I))
      E(I)=X(I+2)-X(I+1)
    1 B(I)=6.*(FD(I+2)-FD(I+1))/(X(I+2)-X(I+1))+6.*(FD(I)-FD(I+1))/(X(I
     *+1)-X(I))
C  CALL THE THOMAS METHOD
      CALL TM(N,C,D,E,B,FX,BETA,GAM)
C  ALLOCATE THE VALUES OF THE SECOND DERIVATIVES TO THE PROPER LOCATION
      DO 5 I=1,N
    5 FPP(I+1)=FX(I)
      FPP(1)=0.0
      FPP(N+2)=0.0
      RETURN
      END
C
C****************************  ABSTRACT  ********************************
C
C     THIS SUBROUTINE SOLVES A SYSTEM OF LINEAR EQUATIONS WITH A TRI-
C  DIAGONAL COEFFICIENT MATRIX USING THE THOMAS METHOD
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
