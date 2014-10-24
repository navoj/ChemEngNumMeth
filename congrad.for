      SUBROUTINE CONGRAD(FUN,DER,N,W,GX,NFUNCT,NDRV,IPRINT,FY)
C
C Fletcher-Reeves conjugate gradient search method.
C
	
	Implicit Double Precision (A-H,O-Z)

	DOUBLE PRECISION W(1),yfr(99),s(99),GX(1),DELX(99)
	EXTERNAL FUN,DER
C
      ITER = 0
      IRESET = N+1
      INDEX = IRESET
C
C Evaluate starting point.
C
      CALL FUN(N,W,FX,NFUNCT)
      CALL DER(N,W,GX,NDRV)
      IF(IPNINT.EQ.1)PRINT 2000, ITER,NFUNCT,NDRV,FX,(W(I),I=1,N)
C
C Calculate squared norm of gradient.
C
   10 SQNOR1 = 0.
      DO 20 I=1,N
   20 SQNOR1 = SQNOR1 + GX(I)*GX(I)
      IF (INDEX.NE.IRESET) GO TO 50
C
C  Set search direction to negative gradient.
C
   30 IF (IPRINT.EQ.1) PRINT 2100
 2100 FORMAT( ' Gradient Step')
      INDEX = 0
      DO 40 I=1,N
   40 S(I) = -GX(I)
      GO TO 70
C
C Set search direction using ratio of squared norms.
C
   50 DO 60 I=1,N
   60 S(I) = -GX(I) + S(I)*SQNOR1/SQNOR2
C
C Find next point.
C
   70 CALL SEARCH(W,N,INDIC,IPRINT,NFUNCT,FX,S,YFR,FY,DELX,GX,ITER,NDRV)
C
C Check whether search was a success.  If not, take a gradient step.
C
      IF (FY.GE.FX) GO TO 30
      CALL DER(N,Yfr,GX,NDRV)
      INDEX = INDEX + 1
      ITER = ITER + 1
      CALL CONVRG(GX,IPASS,W,N,FX,FY,YFR)
      IF(IPASS.EQ.1) GO TO 90
C
C Convergence criteria not satisfied.  Continue search.
C
      IF (IPRINT.EQ.1) PRINT 2000, ITER,NFUNCT,NDRV,FY,(YFR(I),I=1,N)

C
C Save information for next stage.
C
      DO 80 I=1,N
	DELX(I) = YFR(I) - W(I)
   80   W(I) = YFR(I)
      FX = FY
      SQNOR2 = SQNOR1
      GO TO 10
C
C Convergence criteria satisfied.
C
 2000 FORMAT (1X,3I7,E16.8,(17E16.8))
   90 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
C                           Subroutine SEARCH
C
C-----------------------------------------------------------------------------
C
      SUBROUTINE SEARCH(W,N,INDIC,IPRINT,NFUNCT,FX,S,YFR,FY,DELX,GX,
     $ITER,NDRV)
C
C Coggin method of unidimensional search.
C
	
	Implicit Double Precision (A-H,O-Z)
      DOUBLE PRECISION W(10),yfr(10),s(10),DELX(10),GX(10)

C
C The initial vaiables are in X, and the corresponding function value is
C FX.  The search direction vector is S, and the initial step size is STEP.
C
      IEXIT = 0
      NTOL = 0
      FTOL = 0.001
      FTOL2 = FTOL/100.
      FA=FX
      FB=FX
      FC=FX
      DA=0.
      DB=0.
      DC=0.
      K = -2
      M = 0
      STEP = 1.0
      D = STEP
C
C Use the parameter INDIC to indicate how hte search vector length should
C be scaled.
C       If INDIC = 2, then do not scale.  Take length given by mini
C     calculation.
C       If INDIC = 1, then scale only if the length of the last step was
C     shorter than the length of the search vector.  Scale to length of
C     last step.
C       If INDIC = anything but 1 or 2, scale to length of last step.
C
      IF(INDIC.EQ.2.OR.ITER.EQ.0) GO TO 1
C
C  Find norm of S and norm of DELX
C
      DXNORM = 0.
      SNORM = 0.
      DO 102 I=1,N
	DXNORM = DXNORM + DELX(I)*DELX(I)
  102   SNORM = SNORM + S(I)*S(I)
      IF (INDIC.EQ.1.AND.DXNORM.GE.SNORM) GO TO 1
      RATIO = DXNORM/SNORM
      STEP = SQRT(RATIO)
      D = STEP
C
C Start the search the bound the minimum.
C
    1 DO 2 I=1,N
    2   YFR(I) = W(I) + D*S(I)
      CALL FUN(N,Yfr,F,NFUNCT)
      K = K+1
      IF (F-FA) 5,3,6
C
C No change in function value.  Return with vector corresponding to
C function value of FA, because if the function value is independent
C of this search direction,then changes in the variable values may
C upset the main program convergence testing.
C
    3 DO 4 I=1,N
    4  YFR(I) = W(I) + DA*S(I)
      FY = FA
      IF (IPRINT.EQ.1) PRINT 2100
 2100 FORMAT ( ' Search failed.  Function value independent of search
     + direction.')
      GO TO 326
C
C The function is still decreasing.  Increase the step size by double the
C previous increase in step size.
C
    5 FC=FB
      FB=FA
      FA=F
C
      DC=DB
      DB=DA
      DA=D
C
      D = 2.0*D + STEP
      GO TO 1
C
C Minimum is bounded in at least one direction.
C
    6 IF(K) 7,8,9
C
C Minimum is bounded in one direction only.  Reverse the search and direction
C and recycle.
C
    7 FB = F
      DB=D
      D=-D
      STEP=-STEP
      GO TO 1
C
C Minimum is bounded in both directions after only two function
C evaluations, one on either side of the origin.  Proceed to the parabolic
C interpolation.
C
    8 FC=FB
      FB=FA
      FA=F
C
      DC=DB
      DB=DA
      DA=D
C
      GO TO 21
C
C The minimum is bounded after at least two function evaluations in the
C same direction.  Evaluate the function at step size = (DA+DB)/2.
C This will yield 4 equally spaced points bounding the minimum.
C
    9 DC=DB
      DB=DA
      DA=D
C
      FC=FB
      FB=FA
      FA=F
C
   10 D = 0.5*(DA + DB)
      DO 11 I=1,N
   11 YFR(I) = W(I) + D*S(I)
      CALL FUN(N,Yfr,F,NFUNCT)
C
C Now have that FA*FB = FC and that FA*F = FC assuming that the function is
C unimodal.  Remove either point A or point B in such a way that the
C function is bounded and FA*FB = FC.  The corresponding step sizes are
C DA*DB*DC or DA = DB = DC = Z.
C
   12 IF ((DC-D)*(D-DB)) 15,13,18
C
C Location of minimum is limited by rounding errors.  Return with B.
C
   13 DO 14 I=1,N
   14   YFR(I) = W(I)+DB*S(I)
      FY = FB
      IF (IEXIT.EQ.1) GO TO 32
      IF (IPRINT.EQ.1) PRINT 2200
 2200 FORMAT( ' Search failed.  Loc. of minimum limited by rounding.')
      GO TO 325
C
C The point D is in the range DA to DB.
C
   15 IF(F-FB) 16,13,17
   16 FC=FB
      FB=F
      DC=DB
      DB=D
C
      GO TO 21
C
   17 FA=F
      DA=D
C
      GO TO 21
C
C The point D is in the range DB to DC.
C
   18 IF (F-FB) 19,13,20
   19 FA=FB
      FB=F
      DA=DB
      DB=D
C
      GO TO 21
C
   20 FC=F
      DC=D
C
C Now perform the parabolic interpolation.
C
   21 A = FA*(DB - DC) + FB*(DC - DA) + FC*(DA - DB)
      IF (A) 22,30,22
   22 D = 0.5*((DB*DB-DC*DC)*FA+(DC*DC-DA*DA)*FB+(DA*DA-DB*DB)*FC)/A
C
C Check that the point is good.  If so, evaluate the function.
C
      IF((DA-D)*(D-DC)) 13,13,23
   23 DO 24 I=1,N
   24   YFR(I) = W(I) + D*S(I)
      CALL FUN(N,Yfr,F,NFUNCT)
C
C Check for convergence.  If not achieved, recycle.
C
      IF (ABS(FB) - FTOL2) 25,25,26
   25 A = 1.0
      GO TO 27
   26 A = 1.0/FB
   27 IF ((ABS(FB-F)*A) - FTOL) 28,28,12
C
C Convergence achieved.  Return with the smaller F and FB.
C
   28 IEXIT = 1
      IF (F-FB) 29,13,13
   29 FY=F
      GO TO 32
C
C The parabolic interpolation was prevented by the divisor being zero.
C If this is the first time that it has happened, try an intermediate step
C size and recycle.  Otherwise, give up as it looks like a lost cause.
C
   30 IF (M) 31,31,13
   31 M=M+1
      GO TO 10
   32 DO 99 I=1,N
	IF(YFR(1).NE.W(1)) GO TO 325
   99 CONTINUE
      GO TO 33
  325 IF(NTOL.NE.0.AND.IPRINT.EQ.1) PRINT 3000,NTOL
 3000 FORMAT (1X,' Tolerance reduced ',I1,' time(s).')
  326 IF (FY.LT.FX) RETURN
      IF ((S(1).NE.-GX(1)).OR.(FY.LT.FX)) RETURN
      PRINT 5000
 5000 FORMAT ( ' Search failed on a gradient step.  Job terminated.')
      PRINT 5100, ITER,NFUNCT,NDRV,FY,(YFR(I),I=1,N)
 5100 FORMAT (1X,3I7,E16.8,(5E16.8))
      STOP
   33 IF(NTOL.EQ.5) GO TO 34
      IEXIT = 0
      NTOL = NTOL+1
      FTOL = FTOL/10.
      GO TO 12
   34 IF (IPRINT.EQ.1) PRINT 2000
 2000 FORMAT ( ' A pt. better than the entering pt. cannot be found.')
C
      RETURN
      END

c**************************************************************************
c***********************************************************************
C-----------------------------------------------------------------------------
C
C                        Subroutine CONVRG
C
C-----------------------------------------------------------------------------
C
      SUBROUTINE CONVRG(GY,IPASS,W,N,FX,FY,YFR)
	
	Implicit Double Precision (A-H,O-Z)
C
      DOUBLE PRECISION gy(10),W(10),yfr(10)
C
      
      XTOL = 0.00001
      FTOL = 1.E-10
      GTOL = 1.E-6

C        xtol=1
C        ftol=1
C        gtol=1

C
C Check function values.
C
      IF(ABS(FX).LE.FTOL) GO TO 10
      IF(ABS((FX - FY)/FX).GT.FTOL) GO TO 60
      GO TO 20
   10   IF (ABS(FX - FY).GT.FTOL) GO TO 60
C
C Check test point.
C
   20 DO 40 I=1,N
	IF (ABS(W(I)).LE.XTOL) GO TO 30
	IF (ABS((W(I) - YFR(I))/W(I)).GT.XTOL) GO TO 60
	GO TO 40
   30   IF (ABS(W(I) - YFR(I)).GT.XTOL) GO TO 60
   40 CONTINUE
C
C Check gradient.
C
      DO 50 I=1,N
   50 IF (ABS(GY(I)).GT.GTOL) GO TO 60
C
C All convergence criteria satisfied.
C
      IPASS = 1
      RETURN
C
C Convergence not achieved.
C
   60 IPASS = 2
C
      RETURN
      END
