C
C*****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM USES THE LIBRARY ROUTINE CONGRAD TO OPTIMIZE A FORM OF
C  THE ROSENBROCK FUNCTION.  CONGRAD USES BOTH THE FUNCTION VALUES AND
C  THE GRADIENT TO FIND THE OPTIMUM.
C
C****************************  NOMENCLATURE  ****************************
C
C  IPRINT- THE PRINT OPTION FOR CONGRAD
C  N-  THE NUMBER OF UNKNOWNS
C  X(I)-  THE UNKNOWNS
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(10),GRAD(10)
      EXTERNAL FUN,DER
C
C   SET THE INPUT DATA
C
      N=2
      IPRINT=0
C
C   INPUT INITIAL GUESSES
C
      X(1)=10.
      X(2)=10.
C
C  CALL CONGRAD
C
      CALL CONGRAD(FUN,DER,N,X,GRAD,NFUN,NDEV,IPRINT,FY)
C
C  PRINT OUT RESULTS
C
      WRITE(6,2080)NFUN,NDEV
 2080 FORMAT( 10X,' NFUN=',I4,5X,' NDER=',I4)
      WRITE(6,2222)FY,(X(I),I=1,N)
 2222 FORMAT( 5X,' OPTIMUM FUNCTION VALUE=',D12.5,5X,' X=',10D12.5)
      STOP
      END
C
C****************************  ABSTRACT  ********************************
C
C     THIS SUBROUTINE CALCULATES THE VALUE OF F GIVEN THE VALUES OF X(I)
C
C************************************************************************
C
      SUBROUTINE FUN(N,X,F,NFUN)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(10)
      F=10.*(X(2)-X(1)*X(1))**2+(1.-X(1))**2
      NFUN=NFUN+1
      RETURN
      END
C
C****************************  ABSTRACT  ********************************
C
C     THIS SUBROUTINE CALCULATES THE GRADIENT OF F GIVEN THE VALUES OF X(I)
C
C************************************************************************
C
      SUBROUTINE DER(N,X,GRAD,NDER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(10),GRAD(10)
      GRAD(1)=-40.*X(1)*(X(2)-X(1)**2)-2.*(1.-X(1))
      GRAD(2)=20.*(X(2)-X(1)**2)
      NDER=NDER+1
      RETURN
      END
