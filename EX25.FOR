C
C****************************  ABSTRACT  *******************************
C
C    THIS PROGRAM CALL SUBROUTINE SVD TO DETERMINE THE SINGULAR VALUES
C AND THE CONDITION NUMBER OF A COEFICIENT MATRIX.
C
C**************************  NOMENCLATURE  *****************************
C
C  A(I,J)- THE COEFICIENT MATRIX
C  CN- THE CONDITION NUMBER
C  M- THE NUMBER OF ROWS IN THE COEFICIENT MATRIX
C  N- THE NUMBER OF COLUMNS IN THE COEFICIENT MATRIX
C  RV1(I)- A WORK VECTOR USED BY SVD
C  U(I,J)- A RESULTING MATRIX FROM SVD
C  V(I,J)- A RESULTING MATRIX FROM SVD
C  W(I)- THE SINGULAR VALUES OF THE COEFICIENT MATRIX
C
C********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL A(3,3),W(3),U(3,3),V(3,3),RV1(3)  
C
C  SPECIFY NUMBER OF ROWS AND COLUMNS OF THE COEFICIENT MATRIX
C
      M = 3
      N = 3
C
C  SPECIFY THE ELEMENT OF THE COEFICIENT MATRIX
C
      A(1,1) = 3
      A(1,2) = -2
      A(1,3) = 3
      A(2,1) = 4
      A(2,2) = 2
      A(2,3) = -2
      A(3,1) = 1
      A(3,2) = 1
      A(3,3) = 2
C
C  CALL SVD
C
      CALL SVD(M,N,A,W,U,V,RV1,CN)
C
C  WRITE OUT RESULTS
C
      DO 10 I=1,N
   10 WRITE(6,11)I,W(I)
   11 FORMAT( 10X,' SINGULAR VALUE NUMBER',I3,'   = ',F6.2)
      WRITE(6,22)CN
   22 FORMAT( 10X,'     THE CONDITION NUMBER = ',F6.2)
      STOP
      END
