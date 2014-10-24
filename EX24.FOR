C*****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM SOLVES THE MATERIAL BALANCE PROBLEM POSED IN EXAMPLE
C  2.4 USING THE LINEAR EQUATION SOLVER LINPAC.
C
C****************************  NOMENCLATURE  ****************************
C
C  A(I,J)-  THE COEFFICIENT OF THE JTH VARIABLE IN THE ITH EQUATION
C  B(I)-  THE CONSTANT TERM IN THE ITH EQUATION
C  N- THE NUMBER OF INDEPENDENT EQUATIONS
C  X(I)-  THE VECTOR CONTAINING THE SOLUTION TO THE SET OF LINEAR
C         EQUATIONS
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(15,15),B(15),X(15),IPVT(15)
C*** NOTE THAT A,B,X,IPVT MUST BE DOUBLE PRECISIONED AND DIMENSIONED BY N
C
C  INPUT DATA FOR THE PROBLEM
C
      N=15
C  ZERO B(I) AND A(I,J)
      DO 1 I=1,N
      B(I)=0.0
      DO 1 J=1,N
    1 A(I,J)=0.0
C  SET THE NONZERO VALUES OF A(I,J) AND B(I)
      A(1,1)=1.
      A(1,2)=1.
      A(2,1)=1.
      A(2,3)=-1.
      A(2,4)=-1.
      A(2,5)=-1.
      A(2,6)=-1.
      A(3,1)=.65
      A(3,6)=-1.
      A(4,1)=.05
      A(4,5)=-1.
      A(5,4)=1.
      A(5,8)=-1.
      A(6,5)=1.
      A(6,9)=-1.
      A(7,6)=1.
      A(7,10)=-1.
      A(8,7)=1.
      A(8,11)=-1.
      A(9,2)=.3
      A(9,8)=1.
      A(9,13)=-1.
      A(10,2)=.65
      A(10,10)=1.
      A(10,12)=-1.
      A(11,2)=.05
      A(11,9)=1.
      A(11,14)=-1.
      A(12,11)=1.
      A(12,15)=-1.
      A(14,1)=-.2
      A(14,3)=1.
      A(13,12)=-.03
      A(13,13)=-.03
      A(13,14)=-.03
      A(13,15)=.97
      A(15,4)=-.5
      A(15,7)=1.
      B(1)=1000.
C
C  CALL LINPAC
C
      CALL LINPAC(N,A,B,X,IPVT)
C
C  PRINT OUT RESULTS
C
      WRITE(6,88)
   88 FORMAT( ///)
      DO 20 I=1,N
   20 WRITE(6,10)I,X(I)
      WRITE(6,88)
   10 FORMAT( 10X,3H X(,I2,3H)= ,D14.7)
      STOP
      END
