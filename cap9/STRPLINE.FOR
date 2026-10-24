C******PROGRAM STRPLINE.FOR*********************************************
C     To Determine the NxN Capacitance and Inductance Matrices
C     of a N-Conductor Coupled Stripline  Consisting of N Strips of
C     Equal-Width Spaced Midway Between Two Infinite Ground Planes.
C     The Number of Divisions per Strip for the Charge Distribution
C     on each Land is Denoted by NCDIV. The Maximum Problem Size is
C     MSIZE=N*NCDIV.
C     Output Data are in PUL.DAT and Input Data are in STRPLINE.IN
C
C         Written by: Clayton R. Paul
C                     Department of Electrical Engineering
C                     University of Kentucky
C                     Lexington, KY 40506
C
C**********************************************************************
      PARAMETER(MSIZE=100)
      REAL A(MSIZE,MSIZE),CAP(MSIZE,MSIZE),INDUCT,IDENT(MSIZE,MSIZE)
      REAL CAP0(MSIZE,MSIZE)
      INTEGER IPIV(MSIZE),INDXR(MSIZE),INDXC(MSIZE)
      COMMON FAC,DELTA,ER,T,PI
      DO 50 I=1,MSIZE
      DO 51 J=1,MSIZE
51    IDENT(I,J)=0.0
50    IDENT(I,I)=1.0
      OPEN(UNIT=5,FILE='STRPLINE.IN',STATUS='OLD')
      OPEN(UNIT=6,FILE='PUL.DAT',STATUS='UNKNOWN')
      PI=4.0*ATAN(1.0)
      V0=2.997925E8
      V02=V0*V0
      EPS=1.0/(V02*4.E-7*PI)
      FAC=1./(2.*PI*EPS)
      CMTM=2.54E-5
      READ(5,*) N
      READ(5,*) NCDIV
      READ(5,*) W
      W=W*CMTM
      READ(5,*) SEP
      SEP=SEP*CMTM
      READ(5,*) T
      T=T*CMTM
      READ(5,*) ER
C*******************ERROR CHECKING AND INITIAL COMPUTATIONS*************
      IF(N.LT.1) GO TO 999
      IF(N*NCDIV.GT.MSIZE) GO TO 999
      DELTA=W/FLOAT(NCDIV)
C******FILL A1**********************************************************
      DIST=0.
      DO 22 I=1,NCDIV
      A(1,I)=PHIB(DIST)
      A(I,1)=A(1,I)
22    DIST=DIST+DELTA
      DO 23 I=2,NCDIV
      DO 23 J=I,NCDIV
      A(I,J)=A(1,J-I+1)
23    A(J,I)=A(I,J)
C******FILL A2--AN******************************************************
      DO 24 I=2,N
      DIST=FLOAT(I-1)*(W+SEP)
      DO 25 J=1,NCDIV
      XJM1=FLOAT(J-1)
25    A(1,J+(I-1)*NCDIV)=PHIB(DIST+XJM1*DELTA)
      DO 26 J=2,NCDIV
      XJM1=FLOAT(J-1)
26    A(J,(I-1)*NCDIV+1)=PHIB(DIST-XJM1*DELTA)
      DO 27 K=2,NCDIV
      DO 27 M=K,NCDIV
      A(K,M+(I-1)*NCDIV)=A((K-1),(M-1)+(I-1)*NCDIV)
27    A(M,K+(I-1)*NCDIV)=A((M-1),(K-1)+(I-1)*NCDIV)
24    CONTINUE
C******FILL AT2-ATN*****************************************************
      NTOT=N*NCDIV
      DO 28 I=1,NCDIV
      DO 28 J=I,NTOT
28    A(J,I)=A(I,J)
C*******FILL A MATRIX****************************************************
      NM1=N-1
      DO 29 I=1,NM1
      NMI=N-I
      DO 29 J=1,NMI
      DO 29 K=1,NCDIV
      DO 29 L=1,NCDIV
29    A(K+J*NCDIV,L+J*NCDIV+(I-1)*NCDIV)=A(K,L+(I-1)*NCDIV)
      NM2=N-2
      DO 30 I=1,NM2
      NMI=N-(I+1)
      DO 30 J=1,NMI
      DO 30 K=1,NCDIV
      DO 30 L=1,NCDIV
30    A(K+(J+I)*NCDIV,L+J*NCDIV)=A(K+I*NCDIV,L)
C******DETERMINE INDUCTANCE AND CAPACITANCE MATRICES*******************
      CALL GAUSSJ(A,NTOT,MSIZE,IDENT,NTOT,MSIZE,IPIV,INDXR,INDXC)
      DO 32 I=1,N
      DO 32 J=1,N
      SUM=0.0
      DO 33 K=1,NCDIV
      DO 33 M=1,NCDIV
33    SUM=SUM+A(K+(I-1)*NCDIV,M+(J-1)*NCDIV)
32    CAP0(I,J)=SUM
      DO 13 I=1,N
      DO 13 J=1,N
13    CAP(I,J)=CAP0(I,J)
      CALL GAUSSJ(CAP,N,MSIZE,IDENT,N,MSIZE,IPIV,INDXR,INDXC)
      DO 20 I=1,N
      DO 20 J=I,N
      INDUCT=CAP(I,J)/V02
      WRITE(6,21) I,J,INDUCT,I,J
21    FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=L(',I3,',',I3,')')
20    CONTINUE
      DO 41 I=1,N
      DO 41 J=I,N
      CAP(I,J)=ER*CAP0(I,J)
      WRITE(6,42) I,J,CAP(I,J),I,J
42    FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=C(',I3,',',I3,')')
41    CONTINUE
      DO 43 I=1,N
      DO 43 J=I,N
      WRITE(6,44) I,J,CAP0(I,J),I,J
44    FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=C0(',I3,',',I3,')')
43    CONTINUE
      GO TO 997
999   CONTINUE
      WRITE(6,998)
998   FORMAT('INPUT DATA ERROR')
997   CONTINUE
      W=W/CMTM
      SEP=SEP/CMTM
      T=T/CMTM
      WRITE(6,60) N,NCDIV,W,SEP,T,ER
60    FORMAT(///,'NUMBER OF LANDS=',I3/'NUMBER OF DIVISIONS PER LAND='
     1,I3/'LAND WIDTH (mils)=',1PE10.3/
     2'EDGE-TO-EDGE SEPARATION (mils)=',1PE10.3/
     3'SEPARATION BETWEEN GROUND PLANES (mils)=',1PE10.3/
     4'RELATIVE DIELECTRIC CONSTANT=',1PE10.3)
      CLOSE(5)
      CLOSE(6)
      STOP
      END
C******POTENTIAL CALCULATION*******************************************
      REAL FUNCTION PHIB(D)
      COMMON FAC,DELTA,ER,T,PI
C*****4 POINT NUMERICAL INTEGRATION OF THE POTENTIAL FUNCTION
C*****USING THE LEGENDRE-GAUSS GAUSSIAN QUADRATURE METHOD
      X1=DELTA*0.339981/2.
      X2=DELTA*0.861136/2.
      U1=EXP(PI*(D+X1)/T)
      U2=EXP(PI*(D+X2)/T)
      U3=EXP(PI*(D-X2)/T)
      U4=EXP(PI*(D-X1)/T)
      F1=0.652145
      F2=0.347855
      F3=0.347855
      F4=0.652145
      PHIB=F1*ALOG((U1+1./U1+2.)/(U1+1./U1-2.))
     1+F2*ALOG((U2+1./U2+2.)/(U2+1./U2-2.))
     2+F3*ALOG((U3+1./U3+2.)/(U3+1./U3-2.))
     3+F4*ALOG((U4+1./U4+2.)/(U4+1./U4-2.))
      PHIB=PHIB*FAC/4.
      RETURN
      END
C******SUBROUTINE GAUSSJ FOR INVERTING REAL MATRICES*******************
C     Copyright (C) 1987-1992 Numerical Recipes Software
C     Permission is granted for use only with the programs contained
C     within this book
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
