C
C
C*****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM USES THE LIBRARY ROUTINE NMEAD TO OPTIMIZE A FORM OF
C  THE ROSENBROCK FUNCTION.  SINCE THE NEALDER-MEAD SIMPLEX METHOD IS
C  A PATTERN SEARCH METHOD, IT REQUIRES ONLY FUNCTIONS VALUES WHICH ARE
C  SUPPLIED BY SUBROUTINE NSOLV.
C
C****************************  NOMENCLATURE  ****************************
C
C  IPRINT- THE PRINT OPTION FOR THE NEALDER-NEAD METHOD
C  H-  THE INITIAL SIMPLEX SIZE
C  N-  THE NUMBER OF UNKNOWNS
C  X(I)-  THE UNKNOWNS
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /TWO/NFUN
      DIMENSION X(10)
C
C   SET THE INITIAL SIMPLEX SIZE AND SELECT PRINT OPTION
C
      H=.1
      IPRINT=0
C
C   INPUT INITIAL GUESSES
C
      X(1)=10.
      X(2)=10.
C
C  CALL NELDER-MEAD OPTIMIZER
C
      CALL NMEAD(X,2,H,IPRINT)
      WRITE(6,11)NFUN
   11 FORMAT( 10X,' NUMBER OF FUNCTION EVALUATIONS = ',I4)
      STOP
      END
C
C
C****************************  ABSTRACT  ********************************
C
C     THIS SUBROUTINE CALCULATES THE VALUE OF F GIVEN THE VALUES OF X(I)
C
C************************************************************************
C
      SUBROUTINE NSOLV(X,F)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /TWO/NFUN
      DIMENSION X(10)
      F=10.*(X(2)-X(1)*X(1))**2+(1.-X(1))**2
      NFUN=NFUN+1
      RETURN
      END
