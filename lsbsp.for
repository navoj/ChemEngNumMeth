      SUBROUTINE LSBSP(X,Y,BP,YFIT,JORDER,M,MDG,NBP,G,IERR)
      DIMENSION X(1),Y(1),BP(1),YFIT(1),G(MDG,1)
C     REDIMENSION THE FOLLOWING IF NECESSARY
      DIMENSION T(32),C(32),Q(5)
C------------------------------------------------------------------------
      IERR=0
      IF(JORDER.GE.2.AND.JORDER.LE.5) GO TO 5
      IERR=-1
      RETURN
C------------------------------------------------------------------------
    5 JRM1=JORDER-1
C---------- COMPLETE KNOT SEQUENCE T() ----------------------------------
C---------- TAKE CARE OF ENDS; USE SYMMETRIC END INTERVALS --------------
      DO 21 I=1,JRM1
C   LEFT END
      T(JORDER-I)=BP(1)+BP(1)-BP(I+1)
C   RIGHT END
   21 T(NBP+JRM1+I)=BP(NBP)+BP(NBP)-BP(NBP-I)
C-------- FILL MIDDLE OF T() FROM JORDER THROUGH JORDER-1+NBP -----------
      DO 27 I=JORDER,JRM1+NBP
      T(I)=BP(I-JRM1)
C  ORDER BREAK POINTS
      IF(T(I)-T(I-1)) 37,37,27
   27 CONTINUE
C------------------------------------------------------------------------
      IR=1
      IP=1
      JT=1
      J1=JORDER+1
      NBAND=JORDER
C----------------- FIND ONSET OF INTERVAL FOR FIT IN X() ----------------
      I=1
   30 IF(X(I)-T(JORDER)) 35,39,39
   35 I=I+1
      IF(I.LT.M) GO TO 30
   37 IERR=-2
      RETURN
C------------------------------------------------------------------------
   39 I1ST=I
      JB=JORDER
C------------------------------------------------------------------------
   40 MT=0
   50 CONTINUE
      IF(X(I).GT.T(JB+1)) GO TO 60
      CALL BSPL(T,JORDER,X(I),JB,Q)
      IG=IR+MT
      DO 55 K=1,JORDER
   55 G(IG,K)=Q(K)
      G(IG,J1)=Y(I)
      MT=MT+1
      IF(I.EQ.M) GO TO 60
      I=I+1
      GO TO 50
C----------------- SEND BLODK OF DATA TO PROCESSOR ----------------------
   60 CONTINUE
C     WRITE(6,63)MT,I
C  63 FORMAT(1X,'DATA BLOCK SIZE',I4,'; TOTAL POINT=',I4)
      CALL BNDACC(G,MDG,NBAND,IP,IR,MT,JT)
      JB=JB+1
      JT=JT+1
      IF(I.EQ.M.OR.JT.EQ.NBP) GO TO 70
      GO TO 40
C------------------- COMPUTE SOLUTION -----------------------------------
   70 NC=NBAND+NBP-2
      ILAST=I
      CALL BNDSOL(1,G,MDG,NBAND,IP,IR,C,NC,RNORM)
C------------------- COMPUTE FIT ----------------------------------------
      JB=JORDER
      IF(I1ST-1) 79,79,75
   75 DO 77 IG=1,I1ST-1
   77 YFIT(IG)=Y(IG)
   79 DO 110 IG=I1ST,ILAST
   80 IF(X(IG).LE.T(JB+1))GO TO 90
      JB=JB+1
      GO TO 80
   90 CALL BSPL(T,JORDER,X(IG),JB,Q)
      SUM=0.0
      DO 100 L=1,JORDER
      IC=JB-JORDER+L
  100 SUM=SUM+C(IC)*Q(L)
      YFIT(IG)=SUM
  110 CONTINUE
      IF(ILAST-M) 120,199,199
  120 DO 127 IG=ILAST+1,M
  127 YFIT(IG)=Y(IG)
  199 RETURN
      END
      SUBROUTINE BSPL(T,JORDER,X,LEFT,BX)
      DIMENSION T(1),BX(1),DL(7),DR(7)
      J=1
      BX(1)=1.
   20 JP1=J+1
      BX(JP1)=0.0
      DR(J)=T(LEFT+J)-X
      DL(J)=X-T(LEFT+1-J)
      SAVE=0.0
      DO 27 I=1,J
      TERM=BX(I)/(DR(I)+DL(JP1-I))
      BX(I)=SAVE+DR(I)*TERM
   27 SAVE=TERM*DL(JP1-I)
      BX(JP1)=SAVE
      J=JP1
      IF(J-JORDER) 20,99,99
   99 RETURN
      END
      SUBROUTINE BNDSOL(MODE,G,MDG,NB,IP,IR,X,N,RNORM)
      DIMENSION G(MDG,1),X(1)
      ZERO=0.0
      RNORM=0.0
      GO TO (10,90,50),MODE
   10 NBP1=NB+1
      DO 20 J=1,N
   20 X(J)=G(J,NBP1)
      RSQ=ZERO
      NP1=N+1
      IRM1=IR-1
      IF(NP1.GT.IRM1)GO TO 40
      DO 30 J=NP1,IRM1
   30 RSQ=RSQ+G(J,NBP1)**2
      RNORM=SQRT(RSQ)
   40 CONTINUE
C--------- MODE=3 -------------------------------------------------------
   50 DO 80 II=1,N
      I=N+1-II
      S=ZERO
      L=MAX0(0,I-IP)
      IF(I.EQ.N) GO TO 70
      IE=MIN0(N+1-I,NB)
      DO 60 J=2,IE
      JG=J+L
      IX=I-1+J
   60 S=S+G(I,JG)*X(IX)
   70 IF(G(I,L+1)) 80,130,80
   80 X(I)=(X(I)-S)/G(I,L+1)
      RETURN
C------------- MODE=2 ---------------------------------------------------
   90 DO 120 J=1,N
      S=ZERO
      IF(J.EQ.1) GO TO 110
      I1=MAX0(1,J-NB+1)
      I2=J-1
      DO 100 I=I1,I2
      L=J-I+1+MAX0(0,I-IP)
  100 S=S+X(I)*G(I,L)
  110 L=MAX0(0,J-IP)
      IF(G(J,L+1)) 120,130,120
  120 X(J)=(X(J)-S)/G(J,L+1)
      RETURN
  130 WRITE(6,140) MODE,I,J,L
  140 FORMAT(1X,'ERROR: ZERO DIAGONAL TERM IN BNDSOL()'/
     11X,'MODE,I,J,L=',4(1X,I4))
      STOP
      END
      SUBROUTINE BNDACC(G,MDG,NB,IP,IR,MT,JT)
      DIMENSION G(MDG,1)
      ZERO=0.0
C------- ALGOR STEPS 1-4 PERFORMED EXTERNAL TO THIS ROUTINE--------------
      NBP1=NB+1
      IF(MT.LE.0)RETURN
      IF(JT.EQ.IP)GO TO 70
      IF(JT.LE.IR)GO TO 30
      DO 10 I=1,MT
      IG1=JT+MT-I
      IG2=IR+MT-I
      DO 10 J=1,NBP1
      G(IG1,J)=G(IG2,J)
   10 CONTINUE
      IE=JT-IR
      DO 20 I=1,IE
      IG=IR+I-1
      DO 20 J=1,NBP1
      G(IG,J)=ZERO
   20 CONTINUE
      IR=JT
   30 MU=MIN0(NB-1,IR-IP-1)
      IF(MU.EQ.0)GO TO 60
      DO 50 L=1,MU
      K=MIN0(L,JT-IP)
      LP1=L+1
      IG=IP+L
      DO 40 I=LP1,NB
      JG=I-K
   40 G(IG,JG)=G(IG,I)
      DO 50 I=1,K
      JG=NBP1-I
   50 G(IG,JG)=ZERO
   60 IP=JT
   70 MH=IR+MT-IP
      KH=MIN0(NBP1,MH)
      DO 80 I=1,KH
   80 CALL H12(1,I,MAX0(I+1,IR-IP+1),MH,G(IP,I),1,RHO,
     1G(IP,I+1),1,MDG,NBP1-I)
      IR=IP+KH
      IF(KH.LT.NBP1) GO TO 100
      DO 90 I=1,NB
   90 G(IR-1,I)=ZERO
  100 CONTINUE
      RETURN
      END
      SUBROUTINE H12(MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
      DIMENSION U(IUE,M),C(1)
      DOUBLE PRECISION SM,B
      ONE=1.
      IF(0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M)RETURN
      CL=ABS(U(1,LPIVOT))
      IF(MODE.EQ.2)GO TO 60
      DO 10 J=L1,M
   10 CL=AMAX1(ABS(U(1,J)),CL)
      IF(CL) 130,130,20
   20 CLINV=ONE/CL
      SM=(DBLE(U(1,LPIVOT))*CLINV)**2
      DO 30 J=L1,M
   30 SM=SM+(DBLE(U(1,J))*CLINV)**2
      SM1=SM
      CL=CL*SQRT(SM1)
      IF(U(1,LPIVOT)) 50,50,40
   40 CL=-CL
   50 UP=U(1,LPIVOT)-CL
      U(1,LPIVOT)=CL
      GO TO 70
   60 IF(CL) 130,130,70
   70 IF(NCV.LE.0)RETURN
      B=DBLE(UP)*U(1,LPIVOT)
      IF(B) 80,130,130
   80 B=ONE/B
      I2=1-ICV+ICE*(LPIVOT-1)
      INCR=ICE*(L1-LPIVOT)
      DO 120 J=1,NCV
      I2=I2+ICV
      I3=I2+INCR
      I4=I3
      SM=C(I2)*DBLE(UP)
      DO 90 I=L1,M
      SM=SM+C(I3)*DBLE(U(1,I))
   90 I3=I3+ICE
      IF(SM) 100,120,100
  100 SM=SM*B
      C(I2)=C(I2)+SM*DBLE(UP)
      DO 110 I=L1,M
      C(I4)=C(I4)+SM*DBLE(U(1,I))
  110 I4=I4+ICE
  120 CONTINUE
  130 RETURN
      END
