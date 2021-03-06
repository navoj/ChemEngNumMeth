
C
C
C****************************  ABSTRACT  *******************************
C
C
C     THIS PROGRAM SOLVES EXAMPLE 6.12 USING THE IMSL LIBRARY ROUTINE
C  ZX3LP WHICH SOLVES LINEAR PROGRAMING PROBLEM USING THE REVISED
C  SIMPLEX METHOD (EASY TO USE VERSION).
C
C**************************  NOMENCLATURE  ****************************
C
C
C  A(I,J)-  THE COEFICIENT OF THE JTHVARIABLE IN THE ITH CONSTRAINT
C  B(I)-  THE CONSTANT FOR THE ITH CONSTRAINT
C  C(I)-  THE COEFICIENT OF THE ITH VARIABLE IN THE OBJECTIVE FUNCTION
C  IA-  EQUAL TO M1+M2+1
C  IW-  WORK VECTOR
C  M1-  NUMBER OF INEQUALITY CONSTRAINTS
C  M2-  NUMBER OF EQUALITY CONSTRAINTS
C  N-  NUMBER OF INDEPENDENT VARIABLES
C  PSOL(I)- THE VALUE OF THE ITH VARIABLE AT THE OPTIMUM VERTEX
C  RW-  WORK VECTOR
C
C***********************************************************************
      DIMENSION A(10,4),B(10),C(4),RW(128),IW(28),PSOL(20)
      DIMENSION DSOL(10)
C
C  INPUT DATA FOR LP PROBLEM
C
      N=4
      IA=10
      M1=8
      M2=0
C  ZERO A MATRIX
      DO 1 I=1,N
      DO 1 J=1,IA
    1 A(J,I)=0.0
      DO 2 I=1,N
    2 A(I,I)=1.0
C  SPECIFY NON-ZERO VALUES OF A MATRIX
      A(5,1)=.1
      A(5,2)=.6
      A(5,3)=.35
      A(5,4)=.25
      A(6,1)=.3
      A(6,2)=.2
      A(6,3)=.25
      A(6,4)=.29
      A(7,1)=.2
      A(7,2)=.2
      A(7,3)=.28
      A(7,4)=.2
      A(8,1)=.25
      A(8,4)=.15
C  SPECIFY B VECTOR
      B(1)=25.
      B(2)=20.
      B(3)=30.
      B(4)=35.
      B(5)=35.
      B(6)=20.
      B(7)=25.
      B(8)=4.
C  SPECIFY C VECTOR
      C(1)=4100.
      C(2)=4800.
      C(3)=4180.
      C(4)=4520.
C
C  CALL IMSL SUBROUTINE ZX3LP
C
      CALL ZX3LP(A,IA,B,C,N,M1,M2,S,PSOL,DSOL,RW,IW,IER)
C
C  PRINT OUT RESULTS
C
      WRITE(6,11)S
   11 FORMAT( 22H THE MAXIMUM PROFIT IS,F10.2,2X,15HDOLLARS PER DAY)
      WRITE(6,12)
   12 FORMAT( 36H THE OPTIMAL PROCESSING SCHELDULE IS)
   13 FORMAT( 5X,9H CRUDE NO,I2,15H THE FLOW RATE=,E11.4,21H 1000 BARREL
     *S PER DAY)
      DO 20 I=1,N
   20 WRITE(6,13)I,PSOL(I)
      END
