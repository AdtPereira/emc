C******PROGRAM PCB.FOR*************************************************
C     To Determine the NxN Generalized Capacitance Matrix
C     of a N-Conductor Printed Circuit Board.  The (N-1)x(N-1) Transmission
C     Line Capacitance and Inductance Matrices are also Computed from
C     these Results.  The Number of Divisions for the Charge Distribution
C     on each Land is Denoted by NCDIV. The Maximum Problem Size is
C     MSIZE=N*NCDIV.
C     Output Data are in PUL.DAT and Input Data are in PCB.IN.
C
C         Written by: Clayton R. Paul
C                     Department of Electrical Engineering
C                     University of Kentucky
C                     Lexington, KY 40506
C
C**********************************************************************
      PARAMETER(MSIZE=100)
      REAL A(MSIZE,MSIZE)
      REAL CGEN(MSIZE,MSIZE),CAP(MSIZE,MSIZE),INDUCT,IDENT(MSIZE,MSIZE)
      REAL CAP0(MSIZE,MSIZE)
      INTEGER IPIV(MSIZE),INDXR(MSIZE),INDXC(MSIZE)
      COMMON FAC,DELTA,ER,T
      DO 50 I=1,MSIZE
      DO 51 J=1,MSIZE
51    IDENT(I,J)=0.0
50    IDENT(I,I)=1.0
      OPEN(UNIT=5,FILE='PCB.IN',STATUS='OLD')
      OPEN(UNIT=6,FILE='PUL.DAT',STATUS='UNKNOWN')
      PI=4.0*ATAN(1.0)
      V0=2.997925E8
      V02=V0*V0
      EPS=1.0/(V02*4.E-7*PI)
      FAC=1./(2.*PI*EPS)
      CMTM=2.54E-5
      READ(5,*) N
      READ(5,*) NCDIV
      READ(5,*) IREF
      READ(5,*) W
      W=W*CMTM
      READ(5,*) SEP
      SEP=SEP*CMTM
      READ(5,*) T
      T=T*CMTM
      READ(5,*) ER
C*******************ERROR CHECKING AND INITIAL COMPUTATIONS*************
      IF(N.EQ.1.OR.N.EQ.0) GO TO 999
      IF(IREF.GT.N) GO TO 999
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
C******DETERMINE GENERALIZED CAPACITANCE MATRIX WITHOUT DIELECTRIC******
      DO 10 I=1,NTOT
      DO 10 J=1,NTOT
10    CAP(I,J)=A(I,J)
      CALL GAUSSJ(CAP,NTOT,MSIZE,IDENT,NTOT,MSIZE,IPIV,INDXR,INDXC)
      DO 11 I=1,N
      DO 11 J=1,N
      SUM=0.0
      DO 12 K=1,NCDIV
      DO 12 M=1,NCDIV
12    SUM=SUM+CAP(K+(I-1)*NCDIV,M+(J-1)*NCDIV)
11    CGEN(I,J)=SUM*DELTA
C******DETERMINE TRANSMISSION-LINE CAPACITANCE AND INDUCTANCE***********
C******MATRICES WITHOUT DIELECTRIC**************************************
      SUM=0.
      DO 13 I=1,N
      DO 14 J=1,N
14    SUM=SUM+CGEN(I,J)
13    SUMCAP=SUM
      INDEXR=0
      DO 15 I=1,N
      IF(I.EQ.IREF) GO TO 16
      INDEXR=INDEXR+1
      INDEXC=0
      DO 17 J=1,N
      IF(J.EQ.IREF) GO TO 18
      INDEXC=INDEXC+1
      SUMC=0.
      SUMR=0.
      DO 19 II=1,N
      SUMC=SUMC+CGEN(II,J)
19    SUMR=SUMR+CGEN(I,II)
      CAP(INDEXR,INDEXC)=CGEN(I,J)-SUMC*SUMR/SUMCAP
      CAP0(INDEXR,INDEXC)=CAP(INDEXR,INDEXC)
18    CONTINUE
17    CONTINUE
16    CONTINUE
15    CONTINUE
      NM1=N-1
      CALL GAUSSJ(CAP,NM1,MSIZE,IDENT,NM1,MSIZE,IPIV,INDXR,INDXC)
      DO 20 I=1,NM1
      DO 20 J=I,NM1
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
28     A(J,I)=A(I,J)
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
C******DETERMINE GENERALIZED CAPACITANCE MATRIX WITH DIELECTRIC*********
      DO 31 I=1,NTOT
      DO 31 J=1,NTOT
31    CAP(I,J)=A(I,J)
      CALL GAUSSJ(CAP,NTOT,MSIZE,IDENT,NTOT,MSIZE,IPIV,INDXR,INDXC)
      DO 32 I=1,N
      DO 32 J=1,N
      SUM=0.0
      DO 33 K=1,NCDIV
      DO 33 M=1,NCDIV
33    SUM=SUM+CAP(K+(I-1)*NCDIV,M+(J-1)*NCDIV)
32    CGEN(I,J)=SUM*DELTA
C******DETERMINE TRANSMISSION-LINE CAPACITANCE MATRIX WITH DIELECTRIC***
      SUM=0.
      DO 34 I=1,N
      DO 35 J=1,N
35    SUM=SUM+CGEN(I,J)
34    SUMCAP=SUM
      INDEXR=0
      DO 36 I=1,N
      IF(I.EQ.IREF) GO TO 37
      INDEXR=INDEXR+1
      INDEXC=0
      DO 38 J=1,N
      IF(J.EQ.IREF) GO TO 39
      INDEXC=INDEXC+1
      SUMC=0.
      SUMR=0.
      DO 40 II=1,N
      SUMC=SUMC+CGEN(II,J)
40    SUMR=SUMR+CGEN(I,II)
      CAP(INDEXR,INDEXC)=CGEN(I,J)-SUMC*SUMR/SUMCAP
39    CONTINUE
38    CONTINUE
37    CONTINUE
36    CONTINUE
      NM1=N-1
      DO 41 I=1,NM1
      DO 41 J=I,NM1
      WRITE(6,42) I,J,CAP(I,J),I,J
42    FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=C(',I3,',',I3,')')
41    CONTINUE
      DO 43 I=1,NM1
      DO 43 J=I,NM1
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
      WRITE(6,60) N,NCDIV,IREF,W,SEP,T,ER
60    FORMAT(///,'NUMBER OF LANDS=',I3/'NUMBER OF DIVISIONS PER LAND='
     1,I3/'REFERENCE LAND=',I3/'LAND WIDTH (mils)=',1PE10.3/
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
      DP=D+DELTA/2.
      DM=D-DELTA/2.
      IF(D.EQ.0.) GO TO 1
      PHI=FAC*(DM*ALOG(DM)-DP*ALOG(DP)+DELTA)
      GO TO 2
1     PHI=-DELTA*FAC*(ALOG(DELTA/2.)-1.)
2     CONTINUE
      RETURN
      END
C******POTENTIAL WITH THE BOARD ***************************************
      REAL FUNCTION PHIB(D)
      COMMON FAC,DELTA,ER,T
      DP=D+DELTA/2.
      DM=D-DELTA/2.
      DM2=DM*DM
      DP2=DP*DP
      XK=(ER-1.)/(ER+1.)
      IF(D.EQ.0.) GO TO 1
      P0=FAC*(DM*ALOG(DM2)-DP*ALOG(DP2)+2.*DELTA)/(ER+1.)
      FIN=0.
      N=0
100   CONTINUE
      N=N+1
      TEMP=FIN
      XN=FLOAT(N)
      TNT=2.*XN*T
      TNT2=TNT*TNT
      F0=(1.+XK)*FAC/(ER+1.)
      F1=XK**(2*N-1)
      DMOA=DM2/TNT2
      DPOA=DP2/TNT2
      F2=D*ALOG((DMOA+1.)/(DPOA+1.))-(DELTA/2.)*(ALOG(DM2+TNT2)
     1+ALOG(DP2+TNT2))
      F3=2.*DELTA+2.*TNT*(ATAN(DM/TNT)-ATAN(DP/TNT))
      FIN=F0*F1*(F2+F3)
      FIN=FIN+TEMP
      IF(ABS(FIN-TEMP).GE.1.E-5) GO TO 100
      PHIB=P0+FIN
      GO TO 2
1     CONTINUE
      P0=DELTA*FAC*(2.-ALOG((DELTA/2.)**2))/(ER+1.)
      FIN=0.
      N=0
200   CONTINUE
      N=N+1
      TEMP=FIN
      XN=FLOAT(N)
      TNT=2.*XN*T
      TNT2=TNT*TNT
      F1=XK**(2*N-1)
      F2=2.-ALOG((DELTA/2.)**2+TNT2)-(4./DELTA)*TNT*ATAN(DELTA/
     1(2*TNT))
      FIN=FAC*DELTA*(1.+XK)*F1*F2/(ER+1.)
      FIN=FIN+TEMP
      IF(ABS(FIN-TEMP).GE.1.E-5) GO TO 200
      PHIB=P0+FIN
2     CONTINUE
      RETURN
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
