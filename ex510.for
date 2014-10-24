
C
C****************************  ABSTRACT  ********************************
C
C     THIS PROGRAM SOLVES FOR THE CURRENT DISTRIBUTION ON THE LOWER
C  ELECTRODE SHOWN IN FIGURE 5.15.  THE LIBRARY ROUTINE FEC IS USED TO
C  SOLVE FOR THE POTENTIAL DISTRIBUTION AND THE POTENTIAL DISTRIBUTION
C  IS USED TO CALCULATE THE CURRENT DISTRIBUTION.
C
C***************************  NOMENCLATURE  *****************************
C
C  ALAM- LAMBA IN THE GOVERNING DIFERENTIAL EQUATION THAT IS SOLVED IN
C        THE BVP (ALAM .EQ. 0)
C  AMP(I)- THE CURRENT DENSITY AT THE ITH NODE (X-DIRECTION) (AMP/CM**2)
C  CAY- K IN THE GOVERNING DIFFERENTIAL EQUATION THAT IS SOLVED IN THE
C       BVP (CAY .EQ. 1)
C  COORD(I,J)- THE COORDINATES OF THE NODE POINTS USED TO SOLVE THE BVP
C  DVOLT- THE EFFECTIVE VOLTAGE (VOLTS)
C  DY- THE STEP SIZE IN THE Y-DIRECTION (CM)
C  DX- THE STEP SIZE IN THE X-DIRECTION (CM)
C  H1- THE HEIGHT OF THE GEOMETRY ON THE LEFT SIDE (CM)
C  H2- THE HEIGHT OF THE GEOMETRY ON THE RIGHT SIDE (CM)
C  IE(I,J)- THIS MATRIX CONTAINS THE NODE POINT NUMBERS FOR THE ITH
C           ELEMENT
C  IP(I,2)- INDICATES THE NODE TYPE FOR THE ITH NODE
C  ML- THE LOWER BAND WIDTH
C  MU- THE UPPER BAND WIDTH
C  MW- THE MOLECULAR WEIGHT OF THE ANODE
C  NERS- THE NUMBER OF ERRORS IN THE CURRENT AND TEMPERATURE CALCULATIONS
C  NTEST- THE NUMBER OF NODES THAT NOT CONVERGED ON THE ANODE SURFACE
C  NVCT(I)- INPUT VECTOR FOR FEC
C  NX- THE NUMBER OF NODES IN THE X-DIRECTION
C  NY- THE NUMBER OF NODES IN THE Y-DIRECTION
C  TC(I)- CONTAINS THE CONSTANT VALUED NODES FOR FEC
C  Y(I)- THE SEPARATION THE ELECTRODES FOR THE ITH NODE
C  W-  THE WIDTH OF THE GEOMETRY (CM)
C
C************************************************************************
C
      DIMENSION COORD(128,2),IE(32,9),IP(128,2),U(80),NVCT(8),TC(48)
      DIMENSION Y(10),AMP(10)
C
C  INPUT THE PROBLEM SPECIFICATIONS
C
      DVOLT=10.
      CON=1.
      W=1.
      H1=1.
      H2=1.5
C
C  INPUT THE NUMERICAL PARAMETERS
C
      NX=9
      NY=5
      DX=W/FLOAT(NX-1)
C
C  INPUT THE DATA FOR FEC
C
      NVCT(1)=(NX-1)*(NY-1)
      NVCT(2)=NX*NY
      NVCT(3)=2*NX
      NVCT(4)=2
      NVCT(5)=0
      NVCT(6)=0
      NVCT(7)=0
      ALAM=0.0
      CAY=1.0
      F=0.0
      MU=NX+1
      ML=NX+1
C  SET UP IP(K,1)
      K=0
      DO 11 J=1,NY
      DO 11 I=1,NX
      K=K+1
      IP(K,1)=1
      IF((I.EQ.1).OR.(I.EQ.NX))IP(K,1)=2
   11 IF((J.EQ.1).OR.(J.EQ.NY))IP(K,1)=0
C  SET UP IE(K,L)
      K=0
      DO 12 J=2,NX
      DO 12 I=2,NY
      K=K+1
      IE(K,2)=NX*FLOAT(I-1)+J
      IE(K,1)=IE(K,2)-1
      IE(K,3)=NX*FLOAT(I-2)+J
   12 IE(K,4)=IE(K,3)-1
C
C  SET UP TC(I), THE CONSTANT VALUED NODES
C
      DO 15 I=1,NX
      TC(I)=DVOLT
   15 TC(I+NX)=0.0
C
C  SETUP THE UPPER BOUNDARY Y(I)
C
      DO 119 I=1,NX
  119 Y(I)=H1+FLOAT(I-1)/FLOAT(NX-1)*(H2-H1)
C
C  SET UP COORD(K,L)
C
      K=0
      DO 19 I=1,NY
      DO 19 J=1,NX
      DY=Y(J)/FLOAT(NY-1)
      K=K+1
      COORD(K,1)=DX*FLOAT(J-1)
   19 COORD(K,2)=DY*FLOAT(I-1)
C
C  CALL FEC TO SOLVE THE POTENTIAL DISTRIBUTION PROBLEM
C
      CALL FEC(NVCT,COORD,IE,IP,TC,MU,ML,ALAM,CAY,F,U)
C  CALCULATE THE CURRENT DISTRIBUTION
      DO 20 I=1,NX
      DY=Y(I)/FLOAT(NY-1)
      I1=I
      I2=NX+I
      I3=2.*NX+I
      X=-CON*(-3.*U(I1)+4.*U(I2)-U(I3))/(2.*DY)
   20 AMP(I)=X
C  PRINT OUT RESULTS
C
  805 FORMAT( //)
      WRITE(6,805)
  458 FORMAT( 30H THE ELECTRODE SEPARATION (CM))
      WRITE(6,458)
      WRITE(6,457)(Y(I),I=1,NX)
  457 FORMAT(5X,5(E12.5,2X))
      WRITE(6,805)
  452 FORMAT( 37H THE CURRENT DISTRIBUTION (AMP/CM**2))
      WRITE(6,452)
      WRITE(6,457)(AMP(I),I=1,NX)
      STOP
      END
