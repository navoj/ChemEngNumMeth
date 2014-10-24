C
C****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM APPLIES THE  LIBRARY ROUTINE LSODE TO THE
C  INTEGRATION ON A NON-STIFF SET OF ODE'S.  NOMENCLATURE FOR THIS
C  PROGRAM CAN BE FOUND IN FILE LSDP1.
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(3),ATOL(3),RWORK(68),IWORK(20)
      EXTERNAL FEX,JAC
C  INITIAL CONDITIONS
      Y(1)=1.
      Y(2)=0.0
      Y(3)=0.0
      TIN=0.0
C  DGEAR PARAMETER SPECIFICATIONS
      NEQ=3
      RTOL=1.E-6
      ITOL=1
      ATOL(1)=1.E-6
      ITASK=1
      ISTATE=1
      IOPT=0
      LRW=20+NEQ*16
      LIW=20
      MF=10
C  INTEGRATE IN EQUALLY SPACED SEGMENTS USING A DO-LOOP
      DO 10 K=1,10
      TOUT=FLOAT(K)
      CALL LSODE(FEX,NEQ,Y,TIN,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,
     1RWORK,LRW,IWORK,LIW,JAC,MF)
C  PRINT OUT RESULTS
   10 WRITE(6,12)TOUT,Y(1),Y(2),Y(3)
   12 FORMAT( 14X,2HX=,F4.1,2X,5HY(1)=,D11.4,2X,5HY(2)=,D11.4,2X,5HY(3)=
     *,D11.4)
      STOP
      END
C
C*********************  SUBROUTINE FCN  ******************************
C
C     THIS SUBROUTINE CALCULATES THE DERIVATIVE OF EACH DEPENDENT
C  VARIABLE WITH RESPECT TO X, YPRIME(I).
C
C*********************************************************************
C
      SUBROUTINE FEX(N,X,Y,YPRIME)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N
      DIMENSION Y(3),YPRIME(3)
      YPRIME(1)=-.04*Y(1)+1.00*Y(2)*Y(3)
      YPRIME(2)=.04*Y(1)-1.00*Y(2)*Y(3)-3.00*Y(2)*Y(2)
      YPRIME(3)=3.00*Y(2)*Y(2)
      RETURN
      END
C
C***************************  SUBROUTINE FCNJ  ************************
C
C     THIS IS A DUMMY SUBROUTINE SINCE MF IS EQUAL TO 10 AND THE
C  JACOBIAN IS ESTIMATED NUMERICALLY.
C
C*********************************************************************
C
      SUBROUTINE JAC(N,X,Y,ML,MU,PD,NRPD)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N
      REAL Y(3),PD(NRPD,3),X
      RETURN
      END
