C
C*****************************  ABSTRACT  *******************************
C
C     THIS PROGRAM USES THE VARIABLE ORDER LEAST SQUARES SPLINE TO SMOOTH
C  SOME RADIO ACTIVE TRACER DATA FOR A KIDNEY.
C
C***************************  NOMENCLATURE  *****************************
C
C  BP(I)-  THE ITH NODE
C  DT-  THE TIME DIFFERENCE BETWEEN EACH KIDNEY READING
C  G(I,J)-  A WORK VECTOR OF LSBSP
C  JORDER-  THE ORDER OF THE SPLINE USED TO SMOOTH THE DATA
C           IF =2 LINEAR SPLINES; IF =3 PARABOLIC SPLINES;
C           IF =4 CUBIC SPLINES; IF=5 4TH ORDER SPLINES.
C  M-  THE NUMBER OF DATA POINTS USE FOR THE SMOOTHING PROCESS
C  MDG-  THE FIRST DIMENSION OF G(I,J)
C  NBP-  THE NUMBER OF NODES USED
C  ND-  THE NUMBER OF KIDNEY READINGS AVAILABLE
C  NXX-  THE NUMBER OF KIDNEY READINGS TO BE USED
C  RRK(I)-  THE ITH READING FOR THE RIGHT KIDNEY
C  X(I)-  THE VALUE OF THE INDEPENDENT VARIABLE FOR THE ITH DATA POINT
C  Y(I)-  THE VALUE OF THE DEPENDENT VARIABLE FOR THE ITH  DATA POINT
C  YFIT(I)-  THE SMOOTHED VALUE OF THE DEPENDENT VARIABLE AT X(I)
C
C************************************************************************
C
      DIMENSION RRK(200)
      DIMENSION X(200),Y(200),YFIT(200),G(35,6),BP(30)
C
C  SET NUMBER OF DATA POINTS AND TIME STEP
C
      ND=116
      DT=6
C
C  INPUT THE KIDNEY DATA
C
      RRK(1)=0.0
      RRK(2)=0.0
      RRK(3)=0.0
      RRK(4)=451.
      RRK(5)=577.
      RRK(6)=237.
      RRK(7)=220.
      RRK(8)=223.
      RRK(9)=188.
      RRK(10)=170.
      RRK(11)=298.
      RRK(12)=278.
      RRK(13)=301.
      RRK(14)=299.
      RRK(15)=302.
      RRK(16)=340.
      RRK(17)=356.
      RRK(18)=366.
      RRK(19)=364.
      RRK(20)=398.
      RRK(21)=392.
      RRK(22)=397.
      RRK(23)=406.
      RRK(24)=418.
      RRK(25)=407.
      RRK(26)=453.
      RRK(27)=449.
      RRK(28)=456.
      RRK(29)=490.
      RRK(30)=520.
      RRK(31)=479.
      RRK(32)=501.
      RRK(33)=514.
      RRK(34)=535.
      RRK(35)=557.
      RRK(36)=539.
      RRK(37)=555.
      RRK(38)=566.
      RRK(39)=550.
      RRK(40)=594.
      RRK(41)=580.
      RRK(42)=562.
      RRK(43)=585.
      RRK(44)=601.
      RRK(45)=579.
      RRK(46)=602.
      RRK(47)=599.
      RRK(48)=598.
      RRK(49)=627.
      RRK(50)=556.
      RRK(51)=604.
      RRK(52)=589.
      RRK(53)=601.
      RRK(54)=608.
      RRK(55)=601.
      RRK(56)=615.
      RRK(57)=568.
      RRK(58)=585.
      RRK(59)=577.
      RRK(60)=607.
      RRK(61)=583.
      RRK(62)=593.
      RRK(63)=584.
      RRK(64)=558.
      RRK(65)=542.
      RRK(66)=539.
      RRK(67)=544.
      RRK(68)=533.
      RRK(69)=519.
      RRK(70)=527.
      RRK(71)=504.
      RRK(72)=504.
      RRK(73)=500.
      RRK(74)=496.
      RRK(75)=484.
      RRK(76)=472.
      RRK(77)=477.
      RRK(78)=459.
      RRK(79)=458.
      RRK(80)=468.
      RRK(81)=494.
      RRK(82)=437.
      RRK(83)=440.
      RRK(84)=453.
      RRK(85)=441.
      RRK(86)=432.
      RRK(87)=427.
      RRK(88)=457.
      RRK(89)=437.
      RRK(90)=405.
      RRK(91)=403.
      RRK(92)=421.
      RRK(93)=427.
      RRK(94)=391.
      RRK(95)=404.
      RRK(96)=427.
      RRK(97)=401.
      RRK(98)=393.
      RRK(99)=406.
      RRK(100)=398.
      RRK(101)=360.
      RRK(102)=388.
      RRK(103)=385.
      RRK(104)=391.
      RRK(105)=376.
      RRK(106)=397.
      RRK(107)=366.
      RRK(108)=385.
      RRK(109)=366.
      RRK(110)=377.
      RRK(111)=372.
      RRK(112)=399.
      RRK(113)=360.
      RRK(114)=393.
      RRK(115)=396.
      RRK(116)=394.
C
C  SET THE NUMBER OF DATA POINTS TO BE USED AND SETUP X(I) AND Y(I)
C
      NXX=ND-9
      DO 10 I=1,NXX
      X(I)=DT*FLOAT(I-1)
   10 Y(I)=RRK(I+9)
C
C  INPUT THE NODE LOCATIONS
C
      BP(1)=0.0
      BP(2)=10.*DT
      BP(3)=20.*DT
      BP(4)=30.*DT
      BP(5)=40.*DT
      BP(6)=50.*DT
      BP(7)=70.*DT
      BP(8)=90.*DT
      BP(9)=DT*FLOAT(ND-1-9)
C
C  SET JORDER, M, MDG, AND NBP
C
      JORDER=5
      M=ND-9
      MDG=35
      NBP=9
C
C  CALL LSBSP
C
      CALL LSBSP(X,Y,BP,YFIT,JORDER,M,MDG,NBP,G,IERR)
C
C  PRINT OUT RESULTS
C
      DO 20 I=1,NXX
   20 WRITE(6,21)X(I),YFIT(I),Y(I)
   21 FORMAT( 3H X=,F6.1,2X,6H YFIT=,F6.1,2X,3H Y=,F6.1)
      STOP
      END
