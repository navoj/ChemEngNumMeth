C
C****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM APPLIES THE EXPLICIT METHOD TO THE INTEGRATION OF
C  OF AN INITIAL VALUED PARTIAL DIFFERENTIAL EQUATION.  THE IV-PDE
C  FOR THIS PROBLEM REPRESENTS UNSTEADY STATE HEAT CONDUCTION IN A
C  SLAB.
C
C***********************************************************************
C
C
C************************  NOMENCLATURE  *******************************
C
C   ALFA-  THERMAL DIFFUSIVITY (M**2/SEC)
C   DX-  SPATIAL DESCRETIZATION (M)
C   DT-  TIME STEP (SEC)
C   NTMAX-  THE NUMBER  OF TIME STEPS
C   NX- THE NUMBER OF SPATIAL NODES
C   T(I,J)-  TEMPERATUREOF THE ITH SPATIAL NODE AND THE JTH TIME NODE
C   TI- INSIDE TEMPERATURE (DEG K)
C   TO- OUTSIDE TEMPERATURE (DEG K)
C   XL-  THE THICKNESS OF THE WALL (M)
C
C**********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(25,25)
C  SPECIFY PROBLEMS PARAMETERS
      DT=1.
      DX=.125
      ALFA=.005
      TO=0.0
      TI=25.
      TMAX=10.
      XL=.50
C  CALCULATE GRID DIMENSIONS
      NX=IFIX(XL/DX)+1
      NXM=NX-1
      NTMAX=IFIX(TMAX)+1
C  SET INITIAL CONDITIONS
      DO 10 I=1,NX
   10 T(I,1)=TO
C  SET BOUNDARY CONDITIONS
      DO 20 J=1,NTMAX
      T(1,J)=TO
   20 T(NX,J)=TI
      DO 100 J=2,NTMAX
C  PERFORM INTEGRATION
      DO 100 I=2,NXM
  100 T(I,J)=T(I,J-1)+ALFA*DT/DX/DX*(T(I+1,J-1)-2.*T(I,J-1)+T(I-1,J-1))
C  PRINT OUT RESULTS
      WRITE(6,16)
   16 FORMAT( '                     TEMPERATURE  (DEG C) ' )
      WRITE(6,15)
   15 FORMAT( '            ___________________________________________')
   14 FORMAT( '     TIME     X=0   X=.125  X=.25  X=.375   X=.5     ')
      WRITE(6,14)
      DO 200 J=1,NTMAX
  200 WRITE(6,17)J,(T(I,J),I=1,NX)
   17 FORMAT( 3H T=,I3,5X,6(F6.3,2X))
      END
