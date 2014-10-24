C
C****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM APPLIES THE  LIBRARY ROUTINE LSODE TO THE
C  INTEGRATION ON A STIFF SET OF ODE'S.  NOMENCLATURE FOR THIS
C  PROGRAM CAN BE FOUND IN FILE LSDP1.
C
C***********************************************************************
C
      EXTERNAL FEX,JAC
      DOUBLE PRECISION ATOL,RWORK,RTOL,T,TOUT,Y
      DIMENSION Y(3),ATOL(3),RWORK(58),IWORK(23)
C  INITIAL CONDITIONS
      Y(1)=1.
      Y(2)=0.0
      Y(3)=0.0
      T=0.0
      TOUT=1.0
C  LSODE PARAMETER SPECIFICATIONS
      NEQ=3
      RTOL=1.D-4
      ITOL=2
      ATOL(1)=1.D-6
      ATOL(2)=1.D-10
      ATOL(3)=1.D-6
      ITASK=1
      ISTATE=1
      IOPT=0
      LRW=58
      LIW=23
      MF=22
C  INTEGRATE IN EQUALLY SPACED SEGMENTS USING A DO-LOOP
      DO 10 K=1,10
      CALL LSODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     1    IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
      T=TOUT
      TOUT=T+1.
C  PRINT OUT RESULTS
   10 WRITE(6,12)T,Y(1),Y(2),Y(3)
   12 FORMAT( 14X,2HT=,F4.1,2X,5HY(1)=,D11.4,2X,5HY(2)=,D11.4,2X,5HY(3)=
     *,D11.4)
      WRITE(6,899)
      WRITE(6,898)IWORK(11),IWORK(12)
  898 FORMAT( 13X,12H NO. STEPS =,I4,18H   NO. FEX CALLS =,I4)
  899 FORMAT( //)
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
      SUBROUTINE FEX(NEQ,T,Y,YDOT)
      DOUBLE PRECISION T,Y,YDOT
      DIMENSION Y(3),YDOT(3)
      YDOT(1)=-.04*Y(1)+1.D4*Y(2)*Y(3)
      YDOT(2)=.04*Y(1)-1.D4*Y(2)*Y(3)-3.D7*Y(2)*Y(2)
      YDOT(3)=3.D7*Y(2)*Y(2)
      RETURN
      END
C
C***************************  SUBROUTINE FCNJ  ************************
C
C     THIS IS A DUMMY SUBROUTINE SINCE MF IS EQUAL TO 10
C
C*********************************************************************
C
      SUBROUTINE JAC(NEQ,T,Y,ML,MU,PD,NRPD)
      DOUBLE PRECISION PD,T,Y
      DIMENSION Y(3),PD(NRPD,3)
      RETURN
      END
