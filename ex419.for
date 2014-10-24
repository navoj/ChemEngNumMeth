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
      DIMENSION Y(90),ATOL(90),RWORK(999),IWORK(99)
C  INITIAL CONDITIONS
      NEQ=9
      DO 111 I=1,NEQ
  111 Y(I)=1.0
      TIN=0.0
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
      TOUT=FLOAT(K)/10.
      CALL LSODE(FEX,NEQ,Y,TIN,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,
     1RWORK,LRW,IWORK,LIW,JAC,MF)
C  PRINT OUT RESULTS
      TIN=TOUT
      WRITE(6,631)TIN
  631 FORMAT( //10X,' TIME=',F4.1)
   10 WRITE(6,12)(Y(IX),IX=1,NEQ)
   12 FORMAT( 5X,' Y=',5(D11.4,3X))
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
      SUBROUTINE FEX(N,X,Y,YDOT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(9),YDOT(9)
      DX=1./DFLOAT(N+1)
      YDOT(1)=Y(1)*(1.-2.*Y(1)+Y(2))/DX/DX
      YDOT(N)=Y(N)*(Y(N-1)-2.*Y(N))/DX/DX
      NM=N-1
      DO 10 I=2,NM
   10 YDOT(I)=Y(I)*(Y(I-1)-2.*Y(I)+Y(I+1))/DX/DX
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
