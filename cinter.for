C
C*****************************  ABSTRACT  *******************************
C
C     THIS SUBROUTINE APPLIES THE CUBIC SPLINE FORMULA USING THE SECOND
C  DERIVATIVES AND THE INPUT DATA.
C
C***********************************************************************
C
      SUBROUTINE CINTER(ND,X,FD,FPP,XEV,FEV)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(1),FD(1),FPP(1)
C  PERFORM THE INTERPOLATION
      DO 3 J=1,ND
      IF(XEV.LE.X(J))GO TO 4
    3 CONTINUE
    4 FEV=FPP(J-1)/6./(X(J)-X(J-1))*(X(J)-XEV)**3+FPP(J)/6./(
     1X(J)-X(J-1))*(XEV-X(J-1))**3+(FD(J-1)/(X(J)-X(J-1))-FPP(J-1)
     2*(X(J)-X(J-1))/6.)*(X(J)-XEV)+(FD(J)/(X(J)-X(J-1))-FPP(J)*(X
     3(J)-X(J-1))/6.)*(XEV-X(J-1))
      RETURN
      END
