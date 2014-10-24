C
C****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM APPLIES THE METHOD OF LINE TO THE SOLUTION OF AN
C  INITIAL VALUED PARTIAL DIFFERENTIAL EQUATION.  NOMENCLATURE FOR THIS
C  PROGRAM CAN BE FOUND IN FILE LSDP1.
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FEX,JAC
      DIMENSION Y(90),ATOL(90),RWORK(999),IWORK(99),TC(10,5)
C  INITIAL CONDITIONS
      NEQ=3
      DO 111 I=1,NEQ
  111 Y(I)=0.0
      TIN=0.0
      NT=NEQ+2
      DO 43 I=1,NT
   43 TC(1,I)=0.0
      DO 44 I=1,10
      TC(I,1)=0.0
   44 TC(I,NT)=25.
C  SPECIFY LSODE PARAMETERS
      RTOL=1.E-4
      ITOL=2
      DO 112 I=1,NEQ
  112 ATOL(I)=1.E-6
      ITASK=1
      ISTATE=1
      IOPT=0
      LRW=22+NEQ*13
      LIW=20+NEQ
      MF=25
      IWORK(1)=1
      IWORK(2)=1
C  INTEGRATE IN EQUALLY SPACED SEGMENTS USING A DO-LOOP
      DO 10 K=1,10
      TOUT=FLOAT(K)
      CALL LSODE(FEX,NEQ,Y,TIN,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,
     1RWORK,LRW,IWORK,LIW,JAC,MF)
C  PRINT OUT RESULTS
      T=TOUT
      DO 20 J=1,NEQ
   20 TC(K,J+1)=Y(J)
   10 CONTINUE
c
C  PRINT OUT RESULTS
C
      WRITE(6,16)
   16 FORMAT( 10X,'                      TEMPERATURE ')
      WRITE(6,15)
   15 FORMAT( 10X,'          _________________________________________')
   14 FORMAT( 10X,'   TIME       X=0  X=.125   X=.25  X=.375   X=.5')
      WRITE(6,14)
      DO 1000 I=1,10
      WRITE(6,17)I,(TC(I,J),J=1,5)
   17 FORMAT( 10X,3H T=,I3,5X,6(F6.3,2X))
 1000 CONTINUE
      STOP
      END
C
C*********************  SUBROUTINE FCN  ******************************
C
C     THIS SUBROUTINE CALCULATES THE DERIVATIVE OF EACH DEPENDENT
C  VARIABLE WITH RESPECT TO X, YDOT(I).
C
C*********************************************************************
C
      SUBROUTINE FEX(N,X,Y,DYDX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(1),DYDX(1)
      DYDX(1)=-2.*Y(1)+Y(2)
      DYDX(2)= Y(1)-2.*Y(2)+Y(3)
      DYDX(3)=Y(2)-2.*Y(3)+25.
      RETURN
      END
C
C***************************  SUBROUTINE FCNJ  ************************
C
C     THIS IS A DUMMY SUBROUTINE SINCE MF IS EQUAL TO 22 AND THE
C  JACOBIAN IS ESTIMATED NUMERICALLY.
C
C*********************************************************************
C
      SUBROUTINE JAC(N,X,Y,ML,MU,PD,NRPD)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(1),PD(NRPD,1)
      RETURN
      END
