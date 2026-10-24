C******PROGRAM WIDESEP.FOR*********************************************
C     To Determine the NxN Transmission-Line Inductance and Capacitance
C     Matrices of an (N+1) Conductor Transmission Line Consisting of
C     N Wires and a Reference Conductor Using the Wide-Separation
C     Approximations.  There are Three Choices for the Reference
C     Conductor:
C               (1) Another Wire
C               (2) An Infinite, Perfectly-Conducting Ground Plane
C               (3) An Overall, Cylindrical Shield.
C     The Wires are Assumed to be Immersed in a Homogeneous Medium
C     with Relative Permittivity ER and Relative Permeability MUR.
C
C     The Input Data Describing the Cross-sectional Dimensions and
C     the Reference Conductor are contained in the File WIDESEP.IN.
C     The Upper Triangles of the Per-Unit-Length Inductance and
C     Capacitance Matrices are Written to the File PUL.DAT with the
C     Inductance Matrix First Followed by the Capacitance Matrix.
C     The Array Dimensions are Contained in the Parameter MSIZE=N.
C
C         Written by: Clayton R. Paul
C                     Department of Electrical Engineering
C                     University of Kentucky
C                     Lexington, KY 40506
C
C**********************************************************************
      PARAMETER(MSIZE=2)
      INTEGER REFTYPE
      REAL L(MSIZE,MSIZE),LINV(MSIZE,MSIZE),C(MSIZE,MSIZE)
      REAL R(MSIZE),X(MSIZE),Y(MSIZE),MUR,MUO2PI,MUO4PI
      REAL IDENT(MSIZE,MSIZE)
      INTEGER IPIV(MSIZE),INDXR(MSIZE),INDXC(MSIZE)
      DO 1 I=1,MSIZE
      DO 2 J=1,MSIZE
2     IDENT(I,J)=0.0
1     IDENT(I,I)=1.0
      OPEN(UNIT=5,FILE='WIDESEP.IN',STATUS='OLD')
      OPEN(UNIT=6,FILE='PUL.DAT',STATUS='UNKNOWN')
      PI=4.0*ATAN(1.0)
      V0=2.997925E8
      V02=V0*V0
      RADEG=180./PI
      CMTM=2.54E-5
      MUO2PI=2.E-7
      MUO4PI=1.E-7
      READ(5,*) N
      READ(5,*) REFTYPE
      READ(5,*) ER
      READ(5,*) MUR
      IF (REFTYPE.LT.1.OR.REFTYPE.GT.3) GO TO 999
      IF(N.GT.MSIZE) GO TO 999
      READ(5,*) R0
      DO 4 I=1,N
      READ(5,*) NUM
      READ(5,*) R(NUM)
      READ(5,*) X(NUM)
      READ(5,*) Y(NUM)
4     CONTINUE
      GO TO (100,200,300) REFTYPE
C********REFERENCE CONDUCTOR IS A WIRE**********************************
100   CONTINUE
      R0=R0*CMTM
      DO 5 I=1,N
      R(I)=R(I)*CMTM
      DI02=X(I)**2+Y(I)**2
      L(I,I)=MUO2PI*MUR*ALOG(DI02/(R(I)*R0))
5     LINV(I,I)=L(I,I)
      DO 6 I=1,N
      IP=I+1
      DO 7 J=IP,N
      DI02=X(I)**2+Y(I)**2
      DJ02=X(J)**2+Y(J)**2
      DIJ2=(X(I)-X(J))**2+(Y(I)-Y(J))**2
      DI0=SQRT(DI02)
      DJ0=SQRT(DJ02)
      DIJ=SQRT(DIJ2)
      L(I,J)=MUO2PI*MUR*ALOG((DI0*DJ0)/(DIJ*R0))
      L(J,I)=L(I,J)
      LINV(I,J)=L(I,J)
7     LINV(J,I)=L(J,I)
6     CONTINUE
      GO TO 2000
C********REFERENCE CONDUCTOR IS A GROUND PLANE**************************
200   CONTINUE
      DO 15 I=1,N
      R(I)=R(I)*CMTM
      L(I,I)=MUO2PI*MUR*ALOG(2.*X(I)/R(I))
15    LINV(I,I)=L(I,I)
      DO 16 I=1,N
      IP=I+1
      DO 17 J=IP,N
      SIJ2=(X(I)-X(J))**2+(Y(I)-Y(J))**2
      L(I,J)=MUO4PI*MUR*ALOG(1+4.*X(I)*X(J)/SIJ2)
      L(J,I)=L(I,J)
      LINV(I,J)=L(I,J)
17    LINV(J,I)=L(J,I)
16    CONTINUE
      GO TO 2000
C********REFERENCE CONDUCTOR IS AN OVERALL CYLINDRICAL SHIELD***********
300   CONTINUE
      R0=R0*CMTM
      DO 25 I=1,N
      R(I)=R(I)*CMTM
      X(I)=X(I)*CMTM
      Y(I)=Y(I)/RADEG
      L(I,I)=MUO2PI*MUR*ALOG((R0**2-X(I)**2)/(R0*R(I)))
25    LINV(I,I)=L(I,I)
      DO 26 I=1,N
      IP=I+1
      DO 27 J=IP,N
      TOP=(X(I)*X(J))**2+R0**4-2.*X(I)*X(J)*R0*R0*COS(Y(I)-Y(J))
      BOT=(X(I)*X(J))**2+X(J)**4-2.*X(I)*(X(J)**3)*COS(Y(I)-Y(J))
      L(I,J)=MUO2PI*MUR*ALOG((X(J)/R0)*SQRT(TOP/BOT))
      L(J,I)=L(I,J)
      LINV(I,J)=L(I,J)
27    LINV(J,I)=L(J,I)
26    CONTINUE
C********DETERMINE CAPACITANCE MATRIX AS C=MU*EPS*L-1*******************
2000  CONTINUE
      CALL GAUSSJ(LINV,N,MSIZE,IDENT,N,MSIZE,IPIV,INDXR,INDXC)
      DO 30 I=1,N
      DO 30 J=1,N
30    C(I,J)=(MUR*ER/V02)*LINV(I,J)
C********PRINT INDUCTANCE AND CAPACITANCE MATRICES**********************
      DO 40 I=1,N
      DO 40 J=I,N
      WRITE(6,1100) I,J,L(I,J),I,J
1100  FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=L(',I3,',',I3,')')
40    CONTINUE
      DO 50 I=1,N
      DO 50 J=I,N
      WRITE(6,1200) I,J,C(I,J),I,J
1200  FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=C(',I3,',',I3,')')
50    CONTINUE
C********PRINT INPUT DATA***********************************************
      GO TO (400,500,600) REFTYPE
C********REFERENCE CONDUCTOR IS A WIRE**********************************
400   CONTINUE
      R0=R0/CMTM
      WRITE(6,1000) N,ER,MUR,R0
1000  FORMAT(//,'NUMBER OF WIRES=',I3/
     1'RELATIVE PERMITTIVITY OF HOMOGENEOUS MEDIUM=',F5.1/
     2'RELATIVE PERMEABILITY OF HOMOGENEOUS MEDIUM=',F5.1/
     3'REFERENCE CONDUCTOR IS A WIRE WITH RADIUS (mils)=',1PE10.3/)
      WRITE(6,1001)
1001  FORMAT('WIRE #',2X,'WIRE RADIUS (mils)',2X,
     1'X COORDINATE (meters)',2X,'Y COORDINATE (meters)')
      DO 60 I=1,N
      R(I)=R(I)/CMTM
      WRITE (6,1002) I,R(I),X(I),Y(I)
1002  FORMAT(I3,7X,1PE10.3,11X,1PE10.3,13X,1PE10.3)
60    CONTINUE
      GO TO 3000
C********REFERENCE CONDUCTOR IS A GROUND PLANE**************************
500   CONTINUE
      WRITE(6,1010) N,ER,MUR
1010  FORMAT(//,'NUMBER OF WIRES=',I3/
     1'RELATIVE PERMITTIVITY OF HOMOGENEOUS MEDIUM=',F5.1/
     2'RELATIVE PERMEABILITY OF HOMOGENEOUS MEDIUM=',F5.1/
     3'REFERENCE CONDUCTOR IS A GROUND PLANE')
      WRITE(6,1011)
1011  FORMAT('WIRE #',2X,'WIRE RADIUS (mils)',2X,
     1'Y COORDINATE (meters)',2X,'HEIGHT ABOVE GND (meters)')
      DO 70 I=1,N
      R(I)=R(I)/CMTM
      WRITE (6,1012) I,R(I),Y(I),X(I)
1012  FORMAT(I3,7X,1PE10.3,11X,1PE10.3,14X,1PE10.3)
70    CONTINUE
      GO TO 3000
C********REFERENCE CONDUCTOR IS AN OVERALL CYLINDRICAL SHIELD***********
600   CONTINUE
      R0=R0/CMTM
      WRITE(6,1020) N,ER,MUR,R0
1020  FORMAT(//,'NUMBER OF WIRES=',I3/
     1'RELATIVE PERMITTIVITY OF HOMOGENEOUS MEDIUM=',F5.1/
     2'RELATIVE PERMEABILITY OF HOMOGENEOUS MEDIUM=',F5.1/
     3'REFERENCE CONDUCTOR IS A SHIELD WITH RADIUS (mils)=',1PE10.3)
      R0=R0/CMTM
      WRITE(6,1021)
1021  FORMAT('WIRE #',2X,'WIRE RADIUS (mils)',2X,
     1'POSITION RADIUS (mils)',2X,'POSITION ANGLE (degrees)')
      DO 80 I=1,N
      R(I)=R(I)/CMTM
      X(I)=X(I)/CMTM
      Y(I)=Y(I)*RADEG
      WRITE (6,1022) I,R(I),X(I),Y(I)
1022  FORMAT(I3,7X,1PE10.3,11X,1PE10.3,14X,1PE10.3)
80    CONTINUE
      GO TO 3000
C*********************CLOSE FILES***************************************
999   CONTINUE
      WRITE(6,998)
998   FORMAT('INPUT DATA ERROR')
3000  CONTINUE
      CLOSE(5)
      CLOSE(6)
      STOP
      END
C******SUBROUTINE GAUSSJ FOR INVERTING REAL MATRICES*******************
C     Copyright (C) 1987-1992 Numerical Recipes Software
C     Permission is granted for use only with the programs contained
C     within this book.
C
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP,IPIV,INDXR,INDXC)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NP),INDXR(NP),INDXC(NP)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
