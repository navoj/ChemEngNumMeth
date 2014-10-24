      SUBROUTINE SVD (M,N,A,W,U,V,RV1,CN)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL A(M,N),W(N),U(M,N),V(M,N),RV1(N)
      REAL C,F,G,H,X,Y,Z,EPS,SCALE,MACHEP
      REAL DSQRT,AMAX1,DABS,DSIGN

C  ******* MACHEPIS A MACHINE DEPENDENT PARAMETER SPECIFYING 
C          THE RELATIVE PRECISION OF FLOATING POINT ARITHMATIC.
C                     ***********

      MACHEP = 1.E-5

      IERR = 0

      DO 100 I = 1, M

        DO 100 J = 1, N
         U(I,J) = A(I,J)

  100 CONTINUE

C  ********* HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ***********

      G = 0 
      SCALE = 0
      X = 0

      DO 300 I = 1, N
         
         L=I+1
         RV1(I) = SCALE * G
         G = 0
         S= 0
         SCALE = 0
         IF(I.GT.M) GO TO 210

         DO 120 K = I, M
  120	 SCALE = SCALE + DABS(U(K,I))

         IF(SCALE.EQ.0) GO TO 210 

         DO 130 K = I, M
            U(K,I) = U(K,I) / SCALE
            S = S + U(K,I)**2
  130    CONTINUE

         F = U(I,I)
	 G = -DSIGN(DSQRT(S),F)
         H = F * G - S
         U(I,I) = F - G
         IF(I.EQ.N) GO TO 190

         DO 150 J = L, N
            S = 0
         DO 140 K = I, M
  140    S = S + U(K,I) * U(K,J)

         F = S / H

         DO 150 K = I, M
            U(K,J) = U(K,J) + F * U(K,I)
  
  150 CONTINUE
      
  190 DO 200 K = I, M  
  200 U(K,I) = SCALE * U(K,I)

  210 W(I) = SCALE * G
      G = 0
      S = 0
      SCALE = 0
      IF (I .GT. M .OR. I .EQ. N) GO TO 290

      DO 220 K = L, N
  220 SCALE = SCALE + DABS(U(I,K))

      IF(SCALE .EQ. 0) GO TO 290

      DO 230 K = L, N
        U(I,K) = U(I,K) / SCALE
        S = S + U(I,K)**2
  230 CONTINUE

      F = U(I,L)
      G = -DSIGN(DSQRT(S),F)
      H = F * G - S
      U(I,L) = F - G

      DO 240 K = L, N
  240 RV1(K) = U(I,K) / H

      IF (I .EQ. M) GO TO 270 
      
      DO 260 J = L, M
         S = 0

         DO 250 K = L, N
  250    S = S + U(J,K) * U(I,K) 
  
         DO 260 K = L, N
           U(J,K) = U(J,K) + S * RV1(K)
  260 CONTINUE

  270 DO 280 K = L, N
  280 U(I,K) = SCALE * U(I,K)

  290 X = AMAX1(X,DABS(W(I))+DABS(RV1(I)))
  300 CONTINUE


C ******** ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS **********

      DO 400 II = 1, N
        I = N + 1 - II
        IF(I .EQ. N) GO TO 390
        IF(G .EQ. 0) GO TO 360

        DO 320 J = L, N

C **** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ******
  320 V(J,I) = (U(I,J) / U(I,L)) / G

      DO 350 J = L, N
         S = 0

         DO 340 K = L, N
  340    S = S + U(I,K) * V(K,J)

         DO 350 K = L, N 
         V(K,J) = V(K,J) + S * V(K,I)
  350 CONTINUE   

  360 DO 380 J = L, N
         V(I,J) = 0
         V(J,I) = 0
  380 CONTINUE

  390 V(I,I) = 1
      G = RV1(I)
      L = I
  400 CONTINUE

C ******* ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ********
  
C **** FOR I=MIN(M,N) STEP -1 UNTIL 1 DO ********
      MN = N
      IF(M .LT. N) MN = M

      DO 500 II = 1, MN
          I = MN + 1 - II
          L = I + 1
          G = W(I)
          IF (I .EQ. N) GO TO 430

          DO 420 J = L, N
  420     U(I,J) = 0

  430     IF (G .EQ. 0) GO TO 475
          IF (I .EQ. MN) GO TO 460 
          
          DO 450 J = L, N
             S = 0

             DO 440 K = L, M
  440        S = S + U(K,I) * U(K,J)

C ****** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW *******
             F = (S / U(I,I)) / G

             DO 450 K = I, M
             U(K,J) = U(K,J) + F * U(K,I)
  
  450     CONTINUE

  460     DO 470 J = I, M
  470     U(J,I) = U(J,I) / G      
  
          GO TO 490

  475     DO 480 J = I, M
  480     U(J,I) = 0

  490     U(I,I) = U(I,I) + 1

  500 CONTINUE

C ******** DIAGONALISATION OF THE BIDIAGONAL FORM ********

  510 EPS = MACHEP * X

C ******** FOR K = N STEP -1 UNTIL 1 DO -- *******

      DO 700 KK = 1, N
          K1 = N - KK      
          K = K1 + 1
          ITS = 0

C ******** TEST FOR SPLITTING
C                 FOR L=K STEP -1 UNTIL 1 DO-- **********

  520     DO 530 LL = 1, K
             L1 = K - LL
             L = L1 + 1 
	     IF(DABS(RV1(L)) .LE. EPS) GO TO 565
C ***** RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT THROUGH THE BOTTOM 
C              OF THE LOOP *******
	     IF(DABS(W(L1)) .LE. EPS) GO TO 540
  530     CONTINUE

C ******  CANCELLATION OF RV1(L) IF L GREATER THAN 1 **********    

  540     C = 0
          S = 1.0
          
          DO 560 I = L, K
             F = S * RV1(I)
             RV1(I) = C * RV1(I)
	     IF (DABS(F) .LE. EPS) GO TO 565
             G = W(I)
	     H = DSQRT(F*F+G*G)
             W(I) = H
             C = G / H
             S = -F / H
             DO 550 J = 1, M
                Y = U(J,L1)
                Z = U(J,I)
                U(J,L1) = Y * C + Z * S
                U(J,I) = -Y * S + Z * C
  550        CONTINUE    
  
  560     CONTINUE

C ******** TEST FOR CONVERGENCE **********
  565     Z = W(K)
          IF (L .EQ. K) GO TO 650

C ***** SHIFT FROM BOTTOM 2 BY 2 MINOR *******

          IF (ITS .EQ. 30) GO TO 1000
          ITS = ITS + 1
          X = W(L)
          Y = W(K1)
          G = RV1(K1)
          H = RV1(K)
          F = ((Y - Z) * (Y + Z) + (G - H) * (G + H)) / (2. * H * Y)
	  G = DSQRT(F*F+1.)
	  F = ((X - Z) * (X + Z) + H * (Y / (F + DSIGN(G,F)) - H)) / X

C ***** NEXT QR TRANSFORMATION ******* 

          C = 1.
          S = 1.

          DO 600 I1 = L, K1
             I = I1 + 1
             G = RV1(I)
             Y = W(I)
             H = S * G
             G = C * G
	     Z = DSQRT(F*F + H*H)
             RV1(I1) = Z
             C = F / Z
             S = H / Z
             F = X * C + G * S
             G = -X * S + G * C
             H = Y * S
             Y = Y * C
             DO 570 J = 1, N
                X = V(J,I1)
                Z = V(J,I)
                V(J,I1) = X * C + Z * S
                V(J,I) = -X * S + Z * C
  570        CONTINUE
	     Z = DSQRT(F*F+H*H)
             W(I1) = Z

C ******* ROTATION CAN BE ARBITRARY IF Z IS ZERO *******

             IF (Z .EQ. 0) GO TO 580
             C = F / Z
             S = H / Z
  580        F = C * G + S * Y
             X = -S * G + C * Y     
             DO 590 J = 1, M
                Y = U(J,I1)
                Z = U(J,I)
                U(J,I1) = Y * C + Z * S
                U(J,I) = -Y * S + Z * C
  590        CONTINUE

  600     CONTINUE

          RV1(L) = 0
          RV1(K) = F
          W(K) = X
          GO TO 520

C ***** CONVERGENCE ******

  650   IF (Z .GE. 0) GO TO 700     
  
C ****  W(K) IS MADE NONNEGATIVE ********

        W(K) = -Z
        DO 690 J = 1, N
  690   V(J,K) = -V(J,K)
  700	CONTINUE
      GO TO 1001

C ****** SET ERROR -- NO CONVERGENCE TO A SINGULAR VALUE AFTER 
C                     30 ITERATIONS *****                

 1000 WRITE (6,707)
  707 FORMAT( ' PROGRAM DID NOT CONVERGE IN 30 ITERATIONS')
      STOP
 1001 CONTINUE
      WMAX=W(1)
      WMIN=W(1)
      DO 832 I=2,N
      IF(DABS(W(I)).GT.WMAX)WMAX=W(I)
  832 IF(DABS(W(I)).LT.WMIN)WMIN=W(I)
      IF(DABS(WMIN).LT.1.E-20)WMIN=1.E-20
      CN=DABS(WMAX/WMIN)
      RETURN
      END


    
