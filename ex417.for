C
C****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM APPLIES THE CRANK-NICHOLSON METHOD TO THE INTEGRATION
C  OF AN INITIAL VALUED PARTIAL DIFFERENTIAL EQUATION.  THE IV-PDE
C  FOR THIS PROBLEM REPRESENTS UNSTEADY STATE HEAT CONDUCTION IN A
C  SLAB.
C
C***********************************************************************
C
C
C************************  NOMENCLATURE  *******************************
C
C   ALFA-  THERMAL DIFFUSIVITY (M**2/SEC)
C   DX-  SPATIAL DESCRETIZATION (M)
C   DT-  TIME STEP (SEC)
C   NTMAX-  THE NUMBER  OF TIME STEPS
C   NX- THE NUMBER OF SPATIAL NODES
C   T(I,J)-  TEMPERATURE OF THE ITH SPATIAL NODE AND THE JTH TIME NODE
C   TI- INSIDE TEMPERATURE (DEG K)
C   TO- OUTSIDE TEMPERATURE (DEG K)
C   XL-  THE THICKNESS OF THE WALL (M)
C
C**********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(15,15),C(5),D(5),E(5),B(5),X(5),BETA(5),GAM(5)
C  SET SYSTEM PARAMETERS
      DT=1.
      DX=.125
      ALFA=.010
      TO=0.0
      TI=25.
      TMAX=10.
      XL=.50
C  CALCULATE GRID MATRIX
      NX=IFIX(XL/DX)+1
      NXM=NX-1
      NTMAX=IFIX(TMAX/DT)+1
C  SET INITIAL CONDITIONS
      DO 10 I=1,NX
   10 T(I,1)=TO
C  SET BOUNDARY CONDITIONS
      DO 20 J=1,NTMAX
      T(1,J)=TO
   20 T(NX,J)=TI
C  GENERATE SYSTEM OF LINEAR EQUATIONS WITH TRIDIAGONAL COEFICIENT MATRIX
      DO 100 J=2,NTMAX
      DO 200 I=2,NXM
      C(I-1)=-1.
      D(I-1)=2.*(1.+DX*DX/DT/ALFA)
      E(I-1)=-1.
  200 B(I-1)=T(I-1,J-1)+2.*(DX*DX/DT/ALFA-1.)*T(I,J-1)+T(I+1,J-1)
      B(1)=B(1)+T(1,J)
      B(NXM-1)=B(NXM-1)+T(NX,J)
      NXMM=NX-2
C  CALL THOMAS METHOD
      CALL TM(NXMM,C,D,E,B,X,BETA,GAM)
      DO 35 K=1,NXMM
   35 T(K+1,J)=X(K)
  100 CONTINUE
C  PRINT OUT RESULTS
      WRITE(6,16)
   16 FORMAT( '                      TEMPERATURE  (DEG C)  ')
      WRITE(6,15)
   15 FORMAT( '         ____________________________________________')
   14 FORMAT( '  TIME        X=0  X=.125   X=.25  X=.375   X=.5       ')
      WRITE(6,14)
      DO 800 J=1,NTMAX
  800 WRITE(6,17)J,(T(I,J),I=1,NX)
   17 FORMAT( 3H T=,I3,5X,6(F6.3,2X))
      END
C
C*****************************  ABSTRACT  *****************************
C
C     THIS SUBROUTINE USES THE THOMAS METHOD TO SOLVE A SYSTEM OF LINEAR
C   EQUATIONS WITH A TRIADIAGONAL COEFFICIENT MATRIX.  SEE CHAPTER 2 FOR
C   DEFINITION OF TERMS.
C
C**********************************************************************
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

