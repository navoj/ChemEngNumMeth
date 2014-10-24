C**********************   ABSTRACT   **********************************
C
C      THIS PROGRAM CALCULATES THE CURRENT DISTRIBUTION ON THE BASE OF
C   A CYLINDER WITH A COUNTER ELECTRODE ON THE UPPER HALF OF THE OUTSIDE
C   SURFACE.   BY SYMMETRY  CONSIDERATIONS, THE PROBLEM IS REDUCED TO
C   A TWO-DIMENSIONAL BOUNDARY VALUE PROBLEM WHICH IS SOLVED USING
C   FINITE DIFFERENCE APPROXIMATIONS WITH A SUCCESSIVE RELAXATION CON-
C   VERGENCE SCHEME.
C
C**********************   NOMENCLATURE   ******************************
C
C   COND-  ELECTROLYTE CONDUCTIVITY (OHM-M)
C   CUR(I)-  CURRENT DENSITY ON THE BASE ELECTRODE (AMPS/M**2) WHERE
C            I DENOTES THE LOCATION ON THE ELECTRODE
C   DX-  GRID POINT SEPARATION (M)
C   ERLIM-  CONVERGENCE CRITERIA
C   H-  HEIGHT OF THE CYLINDER (M)
C   I-  GRID POINT COLUMN NUMBER
C   ITER- ITERATION NUMBER
C   J-  GRID POINT ROW NUMBER
C   MAXITR-  MAXIMUM NUMBER OF ITERATIONS ALLOW
C   MID-  VALUE OF J FOR WHICH THE ELECTRODE BEGINS ON THE OUTER SURFACE
C         OF THE CYLINDER
C   N-  NUMBER OF GRID POINTS USED IN ANY ROW OR ANY COLUMN
C   NUMERS-  NUMBER OF GRID POINTS THAT HAVE NOT MET THE CONVERGENCE
C            CRITERIA FOR A PARTICULAR ITERATION
C   P(I,J)-  THE ELECTRIC POTENTIAL FOR THE GRID POINT WHERE I AND J
C            INDICATE THE LOCATION(VOLTS)
C   PCAL-  THE VALUE OF THE POTENTIAL AS CALCULATED BY THE RECURSION
C          RELATION(VOLTS)
C   PGUESS-  THE INITIAL GUESS FOR THE POTENTIAL OF ALL THE GRID POINTS
C            (VOLTS)
C   R-  THE RADIAL LOCATION OF A GRID POINT (M)
C   SRF-  SUCCESSIVE RELAXATION FACTOR
C   TEST-  THE RELATIVE CHANGE IN THE POTENTIAL BETWEEN ITERATIONS
C   TOTAL-  TOTAL RELATIVE CHANGE FOR A PARTICULAR ITERATION
C   V-  VOLTAGE DIFFERENCE BETWEEN THE TWO ELECTRODES
C
C*********************   INPUT DESCRIPTION   **************************
C
C      INPUT IS PROVIDED TO THIS PROGRAM BY STATEMENTS AT BEGINING OF
C   THE PROGRAM.  PHYSICAL PARAMETERS AND NUMERICAL PARAMETERS ARE
C   PROVIDED IN SEPARATE SECTIONS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(25,25),CUR(25)
C*******************************************************************************
C****                             INPUT DATA                               *****
C*******************************************************************************
C******************     SET PHYSICAL PARAMETERS
      COND=.002
      H=1.
      R=1.
      V=10.0
C******************     SET NUMERICAL PARAMETERS
      N=12
      DX=R/FLOAT(N-1)
      SRF=1.5
      ERLIM=1.E-5
      PGUESS=V/2.
      MAXITR=1000
C     INITIAL GUESS FOR EACH GRID POINT
      NM=N-1
      DO 10 I=1,N
      DO 10 J=1,N
   10 P(I,J)=PGUESS
C*******************************************************************************
C****                     BEGIN ITERATIVE PROCEDURE                        *****
C*******************************************************************************
C
      ITER=0
  100 ITER=ITER+1
      NUMERS=0
      TOTAL=0.0
C******************     UPDATE BOUNDARY POINTS
C   CENTERLINE
      DO 11 J=2,NM
   11 P(1,J)=(4.*P(2,J)+P(1,J+1)+P(1,J-1))/6.
C   TOP
      DO 12 I=2,NM
   12 P(I,N)=(4.*P(I,N-1)-P(I,N-2))/3.
C   BOTTOM
      DO 13 I=1,N
   13 P(I,1)=0.0
C   OUTSIDE SURFACE
C         INSULATED PART
      MID=N/2-1
      DO 14 I=2,MID
   14 P(N,I)=(4.*P(N-1,I)-P(N-2,I))/3.
C          POSITIVE ELECTRODE
      MID=MID+1
      DO 15 I=MID,N
   15 P(N,I)=V
C******************     UPDATE INTERIOR POINTS
      DO 200 I=2,NM
      DO 200 J=2,NM
      R=DX*FLOAT(I-1)
C     APPLY RECURRSION RELATION
      PCAL=.25*(P(I+1,J)+P(I-1,J)+P(I,J+1)+P(I,J-1))+DX/(8.*R)*(P(I+1,J)
     1-P(I-1,J))
      TEST=DABS((PCAL-P(I,J))/PCAL)
      TOTAL=TOTAL+TEST
      P(I,J)=P(I,J)+SRF*(PCAL-P(I,J))
C     APPLY RELAXATION FACTOR
C     CONVERGENCE CHECK FOR POINT I
  200 IF(TEST.GT.ERLIM)NUMERS=NUMERS+1
C     CHECK TO SEE IF MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED
      IF(ITER.GT.MAXITR)GO TO 250
C     CHECK OVERALL CONVERGENCE
      IF(NUMERS.NE.0)GO TO 100
  250 CONTINUE
C
C*******************************************************************************
C****                         PRINT OUT RESULTS                            *****
C*******************************************************************************
C
  302 FORMAT( ////)
      WRITE(6,302)
      WRITE(6,201)ITER,NUMERS,TOTAL
  201 FORMAT( 7H AFTER ,I3,22H ITERATIONS THERE ARE ,I3,61H POINTS THAT
     1HAVE NOT CONVERGED.  THE TOTAL RELATIVE ERROR IS,D10.3)
C*********************   CALCULATE THE CURRENT DISTRIBUTION
      DO 300 I=1,N
  300 CUR(I)=(3.*P(I,1)-4.*P(I,2)+P(I,3))/(2.*DX)*COND
      WRITE(6,302)
      WRITE(6,308)
  308 FORMAT( 20H RADIAL POSITION (M),5X,27H CURRENT DENSITY (AMP/M**2))
      DO 305 I=1,N
      R=DX*FLOAT(I-1)+DX/2.
  305 WRITE(6,304)R,CUR(I)
  304 FORMAT( 7X,F6.4,19X,D14.7)
      WRITE(6,302)
      STOP
      END

