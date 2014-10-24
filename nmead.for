      SUBROUTINE NMEAD(A,NC,HX,IPRINT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(10)
C
C       PROGRAM SOLVES USING NELDER-MEAD SIMPLEX METHOD.
C
       COMMON /NM1/X(14,10),F(14),N,ICON
      COMMON /NM2/H,ALPHA,BETA,GAMA
       COMMON /ONE/IPRNT
      COMMON /PRINT/IPNT
C
C      H,ALPHA AND BETA ARE NELDER-MEAD PARAMETERS
C
      IPRNT=IPRINT
      IPNT=IPRINT
      H=HX
      ALPHA=1.
      BETA=.5
      GAMA=2.
       N=NC
      DO 100 J=1,N
  100 X(1,J)=A(J)
       ICON=0
       INDEX=1
       NP1=N+1
       NY1=N+2
       NY2=N+3
       NY3=N+4
       DO 20 J=2,NP1
       JM=J-1
       DO 20 K=1,N
       X(J,K)=X(1,K)
       IF(JM.NE.K) GO TO 20
       IF( X(J,K) .EQ. 0.) GO TO 19
       X(J,K)=X(J,K)*(1.+H)
       GO TO 20
   19  X(J,K)=H
   20  CONTINUE
       CALL FUNC(1,NP1)
       CALL ORDER
  201  DO 22 J=1,N
       X(NY1,J)=0.0
       DO 21 I=2,NP1
   21  X(NY1,J)=X(NY1,J)+X(I,J)/N
   22  X(NY2,J)=X(NY1,J)*(1.+ALPHA)-ALPHA*X(1,J)
       CALL FUNC(NY2,NY2)
       IF( F(NY2) .LT. F(NP1) ) GO TO 50
       DO 30 IC1=2,N
       IC = N+2-IC1
       IF( F(NY2) .LT. F(IC) ) GO TO 40
   30  CONTINUE
       IF( F(NY2) .LT. F(1) ) CALL INSERT(NY2,1)
       DO 31 J=1,N
   31  X(NY3,J)=BETA*X(1,J)+(1.-BETA)*X(NY1,J)
       CALL FUNC(NY3,NY3)
       IF ( F(NY3) .LT. F(1) ) GO TO 35
       DO 32 K=1,N
       DO 32 J=1,N
   32  X(K,J)=0.5*(X(K,J)+X(NP1,J))
       CALL FUNC(1,N)
       CALL ORDER
       KIL=5
       GO TO 80
   35  DO 36 IC1=2,NP1
       IC=NP1+2-IC1
       IF ( F(NY3) .LE. F(IC) ) GO TO 37
   36  CONTINUE
       CALL INSERT(NY3,IC)
       GO TO 38
   37  CALL INSERT(NY3,IC)
   38  KIL=4
       GO TO 80
   40  CALL INSERT(NY2,IC)
       KIL=3
       GO TO 80
   50  DO 51 J=1,N
   51  X(NY3,J)=GAMA*X(NY2,J)+(1.-GAMA)*X(NY1,J)
       CALL FUNC(NY3,NY3)
       IF( F(NY3) .LT. F(NY2) ) GO TO 70
       CALL INSERT(NY2,NP1)
       KIL=2
       GO TO 80
   70  CALL INSERT(NY3,NP1)
       KIL=1
   80  CALL ERR(KIL,INDEX)
       IF( INDEX .EQ. 1 ) GO TO 201
       IF(IPRINT.EQ.0)GO TO 357
       WRITE(6,*)(X(NP1,IREC),IREC=1,NP1),F(NP1)
  357  CONTINUE
       CALL FUNC(NP1,NP1)
       DO 555 J=1,N
  555 A(J)=X(NP1,J)
       CALL NSOLV(A,FXC)
       WRITE(6,337)
  337  FORMAT( ///)
       WRITE(6,338)
  338  FORMAT( 5X,37H THE FINAL VALUES AT THE OPTIMUM ARE:)
       DO 339 I=1,NC
  339  WRITE(6,340)I,A(I)
  340  FORMAT( 10X,2HX(,I1,3H)= ,D14.7)
       WRITE(6,341)FXC
  341 FORMAT( 5X,45H THE FINAL VALUE OF THE OBJECTIVE FUNCTION = ,D14.7)
       WRITE(6,337)
       RETURN
       END
C
C
       SUBROUTINE FUNC(INIT,IEND)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /NM1/X(14,10),F(14),N,ICON
       DIMENSION ANS(10)
       DO 10 I=INIT,IEND
C       HERE IS THE PLACE TO COMPUTE YOUR FORMULA
C       EXAMPLE :
C               -  F(I)=SUM(YPRED-YEXP)**2
C                  YEXP CAN BE IN A LOOP
C                  YPRED IS COMPUTED FROM X(I , ......)
C                                              YOUR UNKNOWNS
      DO 100 J=1,N
  100 ANS(J)=X(I,J)
      CALL NSOLV(ANS,F(I))
   10 CONTINUE
       RETURN
       END
C
C
       SUBROUTINE INSERT(NHY,NLOW)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /NM1/X(14,10),F(14),N,ICON
       IF( NLOW .EQ. 1) GO TO 100
       DO 10 I=2,NLOW
       J=I-1
       F(J)=F(I)
       DO 10 K=1,N
   10  X(J,K)=X(I,K)
  100  F(NLOW)=F(NHY)
       DO 11 K=1,N
   11  X(NLOW,K)=X(NHY,K)
       RETURN
       END
C
C
       SUBROUTINE ORDER
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /NM1/X(14,10),F(14),N,ICON
       DIMENSION XA(10)
       NP1=N+1
       DO 20 J=2,NP1
       JK=J-1
       DO 18 K=J,NP1
       IF( F(JK) .GT. F(K) ) GO TO 18
       JK=K
  18   CONTINUE
       JM1=J-1
       DO 19 L=1,N
       XA(L)=X(JM1,L)
       X(JM1,L)=X(JK,L)
  19   X(JK,L)=XA(L)
       FD=F(JM1)
       F(JM1)=F(JK)
   20  F(JK)=FD
       RETURN
       END
C
C
       SUBROUTINE ERR(I,K)
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /NM1/X(14,10),F(14),N,ICON
       COMMON /ONE/IPRINT
       NP1=N+1
       IF( K .EQ. 2 ) GO TO 30
       IF(IPRINT.EQ.0)GO TO 350
       GO TO (12,13,14,15,16) ,I
   12  WRITE(6,121)
  121  FORMAT(' THE SIMPLEX HAS BEEN EXPANDED')
       GO TO 35
   13  WRITE(6,131)
  131  FORMAT(' THE EXPANSION FAILED,EXHANGE OF POINTS')
       GO TO 35
   14  WRITE(6,141)
  141  FORMAT(' THE SIMPLEX REFLECTED SUCCESFULLY')
       GO TO 35
   15  WRITE(6,151)
  151  FORMAT(' THE SIMPLEX HAS BEEN REFLECTED, AND CONTRACTED')
       GO TO 35
   16  WRITE(6,161)
  161  FORMAT(' THE SIMPLEX HAS BEEN CONTRACTED')
   35  WRITE(6,36)
   36  FORMAT('  ')
  350 CONTINUE
      DO 50 JJ=1,N
      DO 60 II=1,N
      IF(DABS(1-X(II+1,JJ)/X(II,JJ)).LT.1D-6) GOTO 60
      GOTO 70
   60 CONTINUE
   50 CONTINUE
      K=2
   70 CONTINUE
   30 RETURN
      END
