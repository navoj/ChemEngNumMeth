C
C***********************   ABSTRACT   ********************************
C
C      THIS SOLVES ONE-DIMENSIONAL BVP'S USING THE SOR METHOD.  THE USER
C   MUST SUPPLY SUBROUTINES RECUR, BCL, AND BCR WHICH APPLY THE RECURSION
C   RELATIONS FOR THE INTERIOR NODES, THE LEFT BOUNDARY, AND THE RIGHT
C   BOUNDARY.  ADDITIONAL INPUT PARAMETERS ARE SUPPLIED AS ARGUEMENTS
C   TO THE CALL OF SUBROUTINE SOR.
C
C**********************   NOMENCLATURE   ******************************
C
C   DX-  SEPARATION BETWEEN GRID POINTS (M)
C   ERLIM-  CONVERGENCE CRITERIA (RELATIVE ERROR)
C   ITER-  ITERATION NUMBER
C   N-  NUMBER OF GRID POINTS
C   MAXITR-  MAXIMUM NUMBER OF ITERATIONS USED IN SEEKING A CONVERGED
C            SOLUTION
C   NUMERS-  NUMBER OF POINTS THAT HAVE NOT MET CONVERGENCE CRITERIA FOR
C            A PARTIULAR ITERATION
C   SRF-  SUCCESSIVE RELAXATION FACTOR
C   TEST-  RELATIVE ERROR FOR A GRID POINT
C   Y(I)-  THE DEPENDENT VARIABLE VALUE AT THE ITH NODE
C   YC-  THE DEPENDENT VARIABLE VALUE CALCULATED BY THE RECURSION RELATION
C   YG-  INITIAL GUESS FOR ALL THE GRID POINTS (DEG C)
C   TOTAL-  TOTAL ERROR
C
C
      SUBROUTINE SOR(RECUR,BCL,BCR,N,SRF,ERLIM,MAXITR,XL,XR,Y,ITER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(1)
      EXTERNAL RECUR,BCL,BCR
C
C*************	 MAKE INITIAL GUESS FOR THE DEPENDENT VARIABLE VALUES
C
      DX=(XR-XL)/FLOAT(N-1)
      ITER=0
C                   *******************************
C*******************   BEGIN ITERATIVE PROCEDURE   *********************
C                   *******************************
  100 ITER=ITER+1
      NUMERS=0
      TOTAL=0.0
C
C*************   UPDATE BOUNDARY POINTS
C
      Y1=Y(1)
      Y2=Y(2)
      Y3=Y(3)
C  UPDATE LEFT BOUNDARY CONDITION
      CALL BCL(Y1,Y2,Y3,DX)
      Y(1)=Y1
      YN=Y(N)
      YNN=Y(N-1)
      YNNN=Y(N-2)
C  UPDATE THE RIGHT BOUNDARY CONDITION
      CALL BCR(YN,YNN,YNNN,DX)
      Y(N)=YN
      NM=N-1
C
C**************   UPDATE INTERIOR POINTS
C
      DO 20 I=2,NM
      YP=Y(I+1)
      YM=Y(I-1)
      YC=Y(I)
      CALL RECUR(YM,YC,YP,DX)
      TEST=ABS((YC-Y(I))/YC)
      TOTAL=TOTAL+DABS(YC-Y(I))
      Y(I)=Y(I)+SRF*(YC-Y(I))
   20 IF(TEST.GT.ERLIM)NUMERS=NUMERS+1
      IF(ITER.GT.MAXITR)GO TO 200
      IF(NUMERS.NE.0)GO TO 100
C
C  CONVERGED SOLUTION; RETURN TO MAIN PROGRAM
C
      RETURN
  200 CONTINUE
C		***************************************
C***************   PRINT OUT RESULTS IF NOT CONVERGED  ******************
C		***************************************
  202 FORMAT( ///)
      WRITE(6,202)
      WRITE(6,201)ITER,NUMERS,TOTAL
  201 FORMAT( 7H AFTER ,I3,22H ITERATIONS THERE ARE ,I3,61H POINTS THAT
     1HAVE NOT CONVERGED.  THE TOTAL RELATIVE ERROR IS,E10.3)
C***********   PRINT OUT DEPENDENT VARIABLE VALUES
      WRITE(6,202)
      WRITE(6,203)
  203 FORMAT( '       X          Y')
      WRITE(6,209)
  209 FORMAT( /)
      X=0.0
      DO 210 I=1,N
      WRITE(6,204)X,Y(I)
  210 X=X+DX
  204 FORMAT( 5X,E11.4,5X,E12.5)
      RETURN
      END
