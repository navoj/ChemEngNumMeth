C
C*****************************  ABSTRACT  *****************************
C
C     THIS SUBROUTINE CALLS SUBROUTINE RK1 AT EACH PRINT INTERVAL
C  AND PRINTS OUT THE INTERMEDIATE RESULTS IF IPNT=1.
C
C****************************  NOMENCLATURE  **************************
C
C  DX-  THE STEP SIZE
C  DXPNT- THE PRINT INTERVAL
C  DXP- THE INTERNALLY CALCULATED PRINT INTERVAL
C  F- THE NAME OF THE SUBROUTINE THAT CALCULATES F(X,Y)
C  IPNT- PRINT FLAG (=0 NO PRINT; =1 PRINT)
C  NS- THE NUMBER OF PRINT INTERVAL STEPS
C  PRO- THE PRINT INTERVAL NEAR XMAX
C  X- THE INDEPENDENT VARIABLE
C  XF- THE VALUE OF X AT THE END OF THE PRINT INTERVAL
C  XO- THE VALUE OF X AT THE BEGINNING OF THE PRINT INTERVAL
C  XMAX-  THE FINAL VALUE OF X
C  Y-  THE DEPENDENT VARIABLE
C
C**********************************************************************
C
      SUBROUTINE RKUTA1(F,DX,X,Y,XMAX,DXPNT,IPNT)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
C  CHECK DXPNT
      IF(DX.GT.DXPNT)WRITE(6,88)
   88 FORMAT( ' DX HAS BEEN SPECIFIED TO BE LESS THAN DXPNT')
      IF(DX.GT.DXPNT)STOP
C  PRINT INITIAL CONDITIONS IF IPNT = 1
      IF(IPNT.EQ.1)WRITE(6,10)X,Y
C  CALCULATE THE NUMBER OF STEPS PER DXPNT
      DXP=DXPNT
      NS=IFIX(XMAX/DXPNT)+1
      XO=X
      PRO=XMAX-DXPNT*DFLOAT(NS-1)
C
C  CALL SUBROUTINE RK1 TIL X=XMAX
C
      DO 1000 I=1,NS
      IF(I.EQ.NS)DXP=PRO
      IF((I.EQ.NS).AND.(PRO.EQ.1.E-6))GO TO 1000
      XF=XO+DXP
C  CALL SUBROUTINE RK1
      CALL RK1(F,DX,XO,XF,Y)
C  PRINT INTERMEDIATE RESULTS IF IPNT=1
      IF(IPNT.EQ.1)WRITE(6,10)XF,Y
   10 FORMAT( 10X,3H X=,D14.7,3X,3H Y=,D14.7)
C  INCREMENT X
 1000 XO=XO+DXP
      RETURN
      END
C
C*************************  ABSTRACT  ********************************
C
C     THIS APPLIES A FOURTH ORDER RUNGE-KUTTA METHOD FOR THE INTEGRATION
C  OF A SINGLE ODE FROM X=XO TO X=XF.
C
C****************************  NOMENCLATURE  *************************
C
C  CAY1-  K1 IN EQUATION 4.5
C  CAY2-  K2 IN EQUATION 4.5
C  CAY3-  K3 IN EQUATION 4.5
C  CAY4-  K4 IN EQUATION 4.5
C  DX-  THE STEP SIZE
C  DXX- INTERNALLY USED STEP SIZE
C  DYDX-  THE DERIVATIVE OF Y WITH RESPECT TO X
C  F- THE NAME OF THE SUBROUTINE THAT CALCULATES F(X,Y)
C  NS- THE NUMBER OF DX INTERVAL STEPS
C  PRO- THE VALUE OF DX NEAR XF
C  X- THE INDEPENDENT VARIABLE
C  XF- THE VALUE OF X AT THE END OF THE PRINT INTERVAL
C  XO- THE VALUE OF X AT THE BEGINNING OF THE PRINT INTERVAL
C  XMAX-  THE FINAL VALUE OF X
C  Y-  THE DEPENDENT VARIABLE
C
C************************************************************************
C
C
      SUBROUTINE RK1(F,DX,XO,XF,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      X=XO
C  CALCULATE THE NUMBER OF STEPS TO REACH XF
      XD=XF-XO
      NS=IFIX(XD/DX)+1
      DXX=DX
      PRO=XD-DX*DFLOAT(NS-1)
      DO 100 I=1,NS
C  CALCULATE DX FOR END OF PRINT INTERVAL
      IF(I.EQ.NS)DXX=PRO
C  CALCULATE K1
      CALL F(X,Y,DYDX)
      CAY1=DYDX
      Y1=Y+DYDX*DXX/2.
      X1=X+DXX/2.
C  CALCULATE K2
      CALL F(X1,Y1,DYDX)
      CAY2=DYDX
      Y1=Y+DYDX*DXX/2.
C  CALCULATE K3
      CALL F(X1,Y1,DYDX)
      CAY3=DYDX
      Y1=Y+DYDX*DXX
      X1=X+DXX
C  CALCULATE K4
      CALL F(X1,Y1,DYDX)
      CAY4=DYDX
C  CALCULATE DX AT END OF INTEGRATION INTERVAL
      IF(XF-X.LT.DXX)DXX=XF-X
C  CALCULATE NEW VALUE OF Y USING RK ALGORITHM
      Y=Y+DXX*(CAY1+2.*CAY2+2.*CAY3+CAY4)/6.
  100 X=X+DXX
      RETURN
      END