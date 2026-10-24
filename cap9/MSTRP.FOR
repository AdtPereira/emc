C******PROGRAM MSTRP.FOR************************************************
C     To Determine the NxN Capacitance and Inductance Matrices
C     of a N-Conductor Coupled Microstrip Consisting of N Strips of
C     Equal-Width above an Infinite Ground Plane.  The Number
C     of Divisions per Strip for the Charge Distribution
C     on each Land is Denoted by NCDIV. The Maximum Problem Size is
C     MSIZE=N*NCDIV.
C     Output Data are in PUL.DAT and Input Data are in MSTRP.IN
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
      COMMON FAC,DELTA,ER,T
      DO 50 I=1,MSIZE
      DO 51 J=1,MSIZE
51    IDENT(I,J)=0.0
50    IDENT(I,I)=1.0
      OPEN(UNIT=5,FILE='MSTRP.IN',STATUS='OLD')
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
C******************WITHOUT DIELECTRIC***********************************
C******FILL A1**********************************************************
      DIST=0.
      DO 1 I=1,NCDIV
      A(1,I)=PHI(DIST)
      A(I,1)=A(1,I)
1     DIST=DIST+DELTA
      DO 2 I=2,NCDIV
      DO 2 J=I,NCDIV
      A(I,J)=A(1,J-I+1)
2     A(J,I)=A(I,J)
C******FILL A2--AN******************************************************
      DO 3 I=2,N
      DIST=FLOAT(I-1)*(W+SEP)
      DO 4 J=1,NCDIV
      XJM1=FLOAT(J-1)
4     A(1,J+(I-1)*NCDIV)=PHI(DIST+XJM1*DELTA)
      DO 5 J=2,NCDIV
      XJM1=FLOAT(J-1)
5     A(J,(I-1)*NCDIV+1)=PHI(DIST-XJM1*DELTA)
      DO 6 K=2,NCDIV
      DO 6 M=K,NCDIV
      A(K,M+(I-1)*NCDIV)=A((K-1),(M-1)+(I-1)*NCDIV)
6     A(M,K+(I-1)*NCDIV)=A((M-1),(K-1)+(I-1)*NCDIV)
3     CONTINUE
C******FILL AT2-ATN*****************************************************
      NTOT=N*NCDIV
      DO 7 I=1,NCDIV
      DO 7 J=I,NTOT
7     A(J,I)=A(I,J)
C*******FILL A MATRIX****************************************************
      NM1=N-1
      DO 8 I=1,NM1
      NMI=N-I
      DO 8 J=1,NMI
      DO 8 K=1,NCDIV
      DO 8 L=1,NCDIV
8     A(K+J*NCDIV,L+J*NCDIV+(I-1)*NCDIV)=A(K,L+(I-1)*NCDIV)
      NM2=N-2
      DO 9 I=1,NM2
      NMI=N-(I+1)
      DO 9 J=1,NMI
      DO 9 K=1,NCDIV
      DO 9 L=1,NCDIV
9     A(K+(J+I)*NCDIV,L+J*NCDIV)=A(K+I*NCDIV,L)
C******DETERMINE CAPACITANCE AND INDUCTANCE MATRICES WITHOUT DIELECTRIC*
      CALL GAUSSJ(A,NTOT,MSIZE,IDENT,NTOT,MSIZE,IPIV,INDXR,INDXC)
      DO 11 I=1,N
      DO 11 J=1,N
      SUM=0.0
      DO 12 K=1,NCDIV
      DO 12 M=1,NCDIV
12    SUM=SUM+A(K+(I-1)*NCDIV,M+(J-1)*NCDIV)
11    CAP0(I,J)=SUM
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
C******************WITH DIELECTRIC**************************************
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
C******DETERMINE CAPACITANCE MATRIX WITH DIELECTRIC**********************
      CALL GAUSSJ(A,NTOT,MSIZE,IDENT,NTOT,MSIZE,IPIV,INDXR,INDXC)
      DO 32 I=1,N
      DO 32 J=1,N
      SUM=0.0
      DO 33 K=1,NCDIV
      DO 33 M=1,NCDIV
33    SUM=SUM+A(K+(I-1)*NCDIV,M+(J-1)*NCDIV)
32    CAP(I,J)=SUM
      DO 41 I=1,N
      DO 41 J=I,N
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
     3'BOARD THICKNESS (mils)=',1PE10.3/
     4'RELATIVE DIELECTRIC CONSTANT=',1PE10.3)
      CLOSE(5)
      CLOSE(6)
      STOP
      END
C******POTENTIAL WITHOUT THE BOARD ************************************
      REAL FUNCTION PHI(D)
      COMMON FAC,DELTA,ER,T
      DOW=2.*D/DELTA
      TOW=4.*T/DELTA
      DOWM=DOW-1.
      DOWP=DOW+1.
      IF(D.EQ.0.) GO TO 1
      PHI=DOWM*ALOG(DOWM)-DOWP*ALOG(DOWP)
     1-0.5*DOWM*ALOG(DOWM**2+TOW**2)
     2+0.5*DOWP*ALOG(DOWP**2+TOW**2)
     3-TOW*(ATAN(DOWM/TOW)-ATAN(DOWP/TOW))
      PHI=PHI*FAC/2.
      GO TO 2
1     PHI=FAC*(0.5*ALOG(1.+TOW*TOW)+TOW*ATAN(1./TOW))
2     CONTINUE
      RETURN
      END
C******POTENTIAL WITH THE BOARD ***************************************
      REAL FUNCTION PHIB(D)
      COMMON FAC,DELTA,ER,T
      TOW=4.*T/DELTA
      DOW=2.*D/DELTA
      DOWM=DOW-1.
      DOWP=DOW+1.
      XK=(ER-1.)/(ER+1.)
      IF(D.EQ.0.) GO TO 1
      P0=2.*XK*ALOG(DELTA/2.)+DOWM*ALOG(DOWM)-DOWP*ALOG(DOWP)
      P0=P0+(XK+1.)*(0.5*DOWP*ALOG(DOWP**2+TOW**2)
     1-0.5*DOWM*ALOG(DOWM**2+TOW**2)
     2-TOW*(ATAN(DOWM/TOW)-ATAN(DOWP/TOW)))
      FIN=0.
      N=1
100   CONTINUE
      N=N+1
      TEMP=FIN
      XN=FLOAT(N)
      TNT=2.*XN*T
      TNT2=2.*TNT/DELTA
      F1=XK**(N-1)
      F1=F1*((-1.)**N)
      F2=-2.*ALOG(DELTA/2.)+0.5*DOWM*ALOG(DOWM**2+TNT2**2)
     1-0.5*DOWP*ALOG(DOWP**2+TNT2**2)
     2+TNT2*(ATAN(DOWM/TNT2)-ATAN(DOWP/TNT2))
      FIN=F1*F2
      FIN=FIN+TEMP
      TEST=ABS(FIN)-ABS(TEMP)
      IF(ABS(TEST).GE.1.E-5) GO TO 100
      PHIB=P0+FIN*(1.+XK)
      PHIB=PHIB*FAC/(ER+1.)
      GO TO 2
1     CONTINUE
      P0=2.*XK*ALOG(DELTA/2.)+(1.+XK)*(ALOG(1.+TOW**2)
     1+2.*TOW*ATAN(1./TOW))
      FIN=0.
      N=1
200   CONTINUE
      N=N+1
      TEMP=FIN
      XN=FLOAT(N)
      TNT=2.*T*XN
      TNT2=TNT**2
      WO2=DELTA/2.
      F1=XK**(N-1)
      F1=F1*((-1.)**N)
      F2=ALOG(WO2**2+TNT2)+(2.*TNT/WO2)*ATAN(WO2/TNT)
      FIN=-F1*F2
      FIN=FIN+TEMP
      TEST=ABS(FIN)-ABS(TEMP)
      IF(ABS(TEST).GE.1.E-5) GO TO 200
      PHIB=P0+FIN*(1.+XK)
      PHIB=PHIB*FAC/(ER+1.)
2     CONTINUE
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
