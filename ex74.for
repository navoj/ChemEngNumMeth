C
C*****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM DETERMINES THE ANTOINE CONSTANTS FOR CARBON DIOXIDE
C  FROM VAPOR PRESSURE DATA USING NON-LINEAR REGRESSION.  THE NELDER-
C  MEAD PATTERN SEARCH OPTIMIZATION METHOD.
C
C****************************  NOMENCLATURE  ****************************
C
C  C(I)-  THE ANTOINE CONSTANTS (UNKNOWN PARAMETERS)
C  H-  THE INITIAL SIMPLEX SIZE
C  N-  THE NUMBER OF DATA POINTS USED
C  P(I)-  THE ITH VAPOR PRESSUREOF CO2 (MM-HG)
C  T(I)-  THE ITH TEMPERATURE (DEG-C)
C
C***********************************************************************
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 C(10)
      COMMON /DATA/T(100),P(100),N
C
C   SET THE INITIAL SIMPLEX SIZE AND SELECT PRINT OPTION
C
      H=.5
      IPRINT=0
C
C   INPUT THE VAPOR PRESSURE/TEMPERATURE DATA
C
      N=32
      DO 1 I=1,32
    1 T(I)=FLOAT(I-1)
      P(1)=26142
      P(2)=26840
      P(3)=27552
      P(4)=28277
      P(5)=29017
      P(6)=29771
      P(7)=30539
      P(8)=31323
      P(9)=32121
      P(10)=32934
      P(11)=33763
      P(12)=34607
      P(13)=35467
      P(14)=36343
      P(15)=37236
      P(16)=38146
      P(17)=39073
      P(18)=40017
      P(19)=40980
      P(20)=41960
      P(21)=42959
      P(22)=43977
      P(23)=45014
      P(24)=46072
      P(25)=47150
      P(26)=48250
      P(27)=49370
      P(28)=50514
      P(29)=51680
      P(30)=52871
      P(31)=54086
      P(32)=55327
C
C   INPUT INITIAL GUESSES FOR THE ANTOINE CONSTANTS
C
      C(1)=20.
      C(2)=1000.
      C(3)=150.
C
C  CALL NELDER-MEAD OPTIMIZER
C
      CALL NMEAD(C,3,H,IPRINT)
C
C  PRINT OUT RESULTS
C
      WRITE(6,22)
      WRITE(6,20)
   20 FORMAT( 35H THE ANTOINE CONSTANTS FOR CO2 ARE:)
      WRITE(6,21)C(1),C(2),C(3)
   21 FORMAT( 3H A=,D12.5,4X,3H B=,D12.5,4X,3H C=,D12.5)
      WRITE(6,22)
   22 FORMAT( ///)
      STOP
      END
C
C
C****************************  ABSTRACT  ********************************
C
C     THIS SUBROUTINE CALCULATES THE SUM OF THE SQUARES OF THE ERROR
C  GIVEN THE ANTOINE CONSTANTS, C(I).
C
C************************************************************************
C
      SUBROUTINE NSOLV(X,F)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(3),DP(100)
      COMMON /DATA/T(100),P(100),N
      F=0.
      DO 10 I=1,N
C  NOTE THE EXPONENTIAL FORM OF THE ANTOINE EQUATION IS USED HERE
      DP(I)=P(I)-DEXP(X(1)-X(2)/(T(I)+273.15+X(3)))
   10 F=F+DP(I)*DP(I)
      RETURN
      END
