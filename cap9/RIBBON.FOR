C******PROGRAM RIBBON.FOR**********************************************
C     To Determine the NxN Generalized Capacitance Matrix
C     of a N-Wire Ribbon Cable.  The (N-1)x(N-1) Transmission-Line
C     Capacitance and Inductance Matrices are also Computed from these
C     Results.  The Number of Fourier Coefficients for the Charge
C     Distribution Around the Wire Peripheries is Denoted by NF.  The
C     Maximum Problem Size is MSIZE=N*NF.  Output Data are
C     PUL.DAT and Input Data are RIBBON.IN.
C
C         Written by: Clayton R. Paul
C                     Department of Electrical Engineering
C                     University of Kentucky
C                     Lexington, KY 40506
C
C**********************************************************************
      PARAMETER(MSIZE=100)
      REAL A(MSIZE,MSIZE),B(MSIZE,MSIZE),C(MSIZE,MSIZE),D(MSIZE,MSIZE)
      REAL CGEN(MSIZE,MSIZE),CAP(MSIZE,MSIZE),INDUCT,IDENT(MSIZE,MSIZE)
      REAL CAP0(MSIZE,MSIZE)
      INTEGER IPIV(MSIZE),INDXR(MSIZE),INDXC(MSIZE)
      DO 50 I=1,MSIZE
      DO 51 J=1,MSIZE
51    IDENT(I,J)=0.0
50    IDENT(I,I)=1.0
      OPEN(UNIT=5,FILE='RIBBON.IN',STATUS='OLD')
      OPEN(UNIT=6,FILE='PUL.DAT',STATUS='UNKNOWN')
      PI=4.0*ATAN(1.0)
      V0=2.997925E8
      V02=V0*V0
      EPS=1.0/(V02*4.E-7*PI)
      CMTM=2.54E-5
      READ(5,*) N
      READ(5,*) NF
      READ(5,*) IREF
      READ(5,*) RW
      RW=RW*CMTM
      READ(5,*) TD
      TD=TD*CMTM
      READ(5,*) ER
      READ(5,*) S
      S=S*CMTM
C*******************ERROR CHECKING AND INITIAL COMPUTATIONS*************
      IF(N.EQ.1.OR.N.EQ.0) GO TO 999
      IF(IREF.GT.N) GO TO 999
      IF((2.*(RW+TD)).GT.S) GO TO 999
      IF(N*NF.GT.MSIZE) GO TO 999
      RD=RW+TD
      ANGLE=2.0*PI/FLOAT(NF)
      DELTA=ANGLE/4.0
C******FILL A1 AND B1***************************************************
      DO 1 I=1,NF
      ANG=FLOAT(I-1)*ANGLE+DELTA
      DO 1 J=1,NF
      IF(J.NE.1) GO TO 2
      A(I,J)=-RW*ALOG(RW)/EPS
      B(I,J)=-RD*ALOG(RD)/EPS
      GO TO 1
2     XM1=FLOAT(J-1)
      A(I,J)=(RW**J)*COS(XM1*ANG)/(2.*EPS*XM1*(RW**(J-1)))
      B(I,J)=(RW**(J-1))*COS(XM1*ANG)/(2.*EPS*XM1*(RD**(J-2)))
1     CONTINUE
C******FILL A2--AN AND B2--BN*******************************************
      DO 3 I=2,N
      SEP=FLOAT(I-1)*S
      DO 3 J=1,NF
      ANG=FLOAT(J-1)*ANGLE+DELTA
      HORZ=SEP-RW*COS(ANG)
      VERT=RW*SIN(ANG)
      TH=PI-ATAN2(VERT,HORZ)
      RP=SQRT(HORZ**2+VERT**2)
      DO 3 K=1,NF
      IF(K.NE.1) GO TO 4
      A(J,K+(I-1)*NF)=-RW*ALOG(RP)/EPS
      B(J,K+(I-1)*NF)=-RD*ALOG(RP)/EPS
      GO TO 3
4     XKM1=FLOAT(K-1)
      A(J,K+(I-1)*NF)=(RW**K)*COS(XKM1*TH)/(2.*EPS*XKM1*(RP**(K-1)))
      B(J,K+(I-1)*NF)=(RD**K)*COS(XKM1*TH)/(2.*EPS*XKM1*(RP**(K-1)))
3     CONTINUE
C******FILL AM2--AMN AND BM2--BMN***************************************
      DO 5 I=2,N
      SEP=FLOAT(I-1)*S
      DO 5 J=1,NF
      ANG=FLOAT(J-1)*ANGLE+DELTA
      HORZ=SEP+RW*COS(ANG)
      VERT=RW*SIN(ANG)
      TH=ATAN2(VERT,HORZ)
      RP=SQRT(HORZ**2+VERT**2)
      DO 5 K=1,NF
      IF(K.NE.1) GO TO 6
      A(J+(I-1)*NF,K)=-RW*ALOG(RP)/EPS
      B(J+(I-1)*NF,K)=-RD*ALOG(RP)/EPS
      GO TO 5
6     XKM1=FLOAT(K-1)
      A(J+(I-1)*NF,K)=(RW**K)*COS(XKM1*TH)/(2.*EPS*XKM1*(RP**(K-1)))
      B(J+(I-1)*NF,K)=(RD**K)*COS(XKM1*TH)/(2.*EPS*XKM1*(RP**(K-1)))
5     CONTINUE
C******FILL C1 AND D1***************************************************
      DO 7 I=1,NF
      ANG=FLOAT(I-1)*ANGLE+DELTA
      DO 7 J=1,NF
      IF(J.NE.1) GO TO 8
      C(I,J)=(ER-1.)*RW/RD
      D(I,J)=-1.0
      GO TO 7
8     XM1=FLOAT(J-1)
      C(I,J)=(ER-1.)*((RW/RD)**J)*COS(XM1*ANG)/2.0
      D(I,J)=-(ER+1.0)*COS(XM1*ANG)/2.0
7     CONTINUE
C******FILL C2--CN AND D2--DN*******************************************
      DO 9 I=2,N
      SEP=FLOAT(I-1)*S
      DO 9 J=1,NF
      ANG=FLOAT(J-1)*ANGLE+DELTA
      HORZ=SEP-RD*COS(ANG)
      VERT=RD*SIN(ANG)
      TH=PI-ATAN2(VERT,HORZ)
      RP=SQRT(HORZ**2+VERT**2)
      RDN=COS(TH-ANG)
      TDN=-SIN(TH-ANG)
      DO 9 K=1,NF
      IF(K.NE.1) GO TO 10
      C(J,K+(I-1)*NF)=(ER-1.)*(RW/RP)*RDN
      D(J,K+(I-1)*NF)=(ER-1.)*(RD/RP)*RDN
      GO TO 9
10    XKM1=FLOAT(K-1)
      C(J,K+(I-1)*NF)=(ER-1.)*(((RW/RP)**K)/2.)*(COS(XKM1*TH)*RDN
     1+SIN(XKM1*TH)*TDN)
      D(J,K+(I-1)*NF)=(ER-1.)*(((RD/RP)**K)/2.)*(COS(XKM1*TH)*RDN
     1+SIN(XKM1*TH)*TDN)
9     CONTINUE
C******FILL CM2--CMN AND DM2--DMN***************************************
      DO 11 I=2,N
      SEP=FLOAT(I-1)*S
      DO 11 J=1,NF
      ANG=FLOAT(J-1)*ANGLE+DELTA
      HORZ=SEP+RD*COS(ANG)
      VERT=RD*SIN(ANG)
      TH=ATAN2(VERT,HORZ)
      RP=SQRT(HORZ**2+VERT**2)
      RDN=COS(TH-ANG)
      TDN=-SIN(TH-ANG)
      DO 11 K=1,NF
      IF(K.NE.1) GO TO 12
      C(J+(I-1)*NF,K)=(ER-1.)*(RW/RP)*RDN
      D(J+(I-1)*NF,K)=(ER-1.)*(RD/RP)*RDN
      GO TO 11
12    XKM1=FLOAT(K-1)
      C(J+(I-1)*NF,K)=(ER-1.)*(((RW/RP)**K)/2.)*(COS(XKM1*TH)*RDN
     1+SIN(XKM1*TH)*TDN)
      D(J+(I-1)*NF,K)=(ER-1.)*(((RD/RP)**K)/2.)*(COS(XKM1*TH)*RDN
     1+SIN(XKM1*TH)*TDN)
11    CONTINUE
C******FILL A,B,C,D MATRICES********************************************
      NM1=N-1
      DO 13 I=1,NM1
      NMI=N-I
      DO 13 J=1,NMI
      DO 13 K=1,NF
      DO 13 L=1,NF
      A(K+J*NF,L+J*NF+(I-1)*NF)=A(K,L+(I-1)*NF)
      B(K+J*NF,L+J*NF+(I-1)*NF)=B(K,L+(I-1)*NF)
      C(K+J*NF,L+J*NF+(I-1)*NF)=C(K,L+(I-1)*NF)
13    D(K+J*NF,L+J*NF+(I-1)*NF)=D(K,L+(I-1)*NF)
      NM2=N-2
      DO 14 I=1,NM2
      NMI=N-(I+1)
      DO 14 J=1,NMI
      DO 14 K=1,NF
      DO 14 L=1,NF
      A(K+(J+I)*NF,L+J*NF)=A(K+I*NF,L)
      B(K+(J+I)*NF,L+J*NF)=B(K+I*NF,L)
      C(K+(J+I)*NF,L+J*NF)=C(K+I*NF,L)
14    D(K+(J+I)*NF,L+J*NF)=D(K+I*NF,L)
C******DETERMINE GENERALIZED CAPACITANCE MATRIX WITHOUT DIELECTRIC******
      NNF=N*NF
      DO 101 I=1,NNF
      DO 101 J=1,NNF
101   CAP(I,J)=A(I,J)
      CALL GAUSSJ(CAP,NNF,MSIZE,IDENT,NNF,MSIZE,IPIV,INDXR,INDXC)
      DO 15 I=1,N
      DO 15 J=1,N
      SUM=0.0
      DO 16 K=1,NF
16    SUM=SUM+CAP((I-1)*NF+1,K+(J-1)*NF)
      CGEN(I,J)=2.*RW*PI*SUM
15    CONTINUE
C******DETERMINE TRANSMISSION-LINE CAPACITANCE AND INDUCTANCE***********
C******MATRICES WITHOUT DIELECTRIC**************************************
      SUM=0.
      DO 17 I=1,N
      DO 18 J=1,N
18    SUM=SUM+CGEN(I,J)
17    SUMCAP=SUM
      INDEXR=0
      DO 19 I=1,N
      IF(I.EQ.IREF) GO TO 20
      INDEXR=INDEXR+1
      INDEXC=0
      DO 21 J=1,N
      IF(J.EQ.IREF) GO TO 22
      INDEXC=INDEXC+1
      SUMC=0.
      SUMR=0.
      DO 23 II=1,N
      SUMC=SUMC+CGEN(II,J)
23    SUMR=SUMR+CGEN(I,II)
      CAP(INDEXR,INDEXC)=CGEN(I,J)-SUMC*SUMR/SUMCAP
      CAP0(INDEXR,INDEXC)=CAP(INDEXR,INDEXC)
22    CONTINUE
21    CONTINUE
20    CONTINUE
19    CONTINUE
      NM1=N-1
      CALL GAUSSJ(CAP,NM1,MSIZE,IDENT,NM1,MSIZE,IPIV,INDXR,INDXC)
      DO 24 I=1,NM1
      DO 24 J=1,NM1
      INDUCT=CAP(I,J)/V02
      IF(J.LT.I) GO TO 333
      WRITE(6,25) I,J,INDUCT,I,J
25    FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=L(',I3,',',I3,')')
333   CONTINUE
24    CONTINUE
C******DETERMINE GENERALIZED CAPACITANCE MATRIX WITH DIELECTRIC*********
      CALL GAUSSJ(D,NNF,MSIZE,IDENT,NNF,MSIZE,IPIV,INDXR,INDXC)
      DO 26 I=1,NNF
      DO 26 J=1,NNF
      SUM=0.
      DO 27 K=1,NNF
27    SUM=SUM+D(I,K)*C(K,J)
26    CAP(I,J)=SUM
      DO 28 I=1,NNF
      DO 28 J=1,NNF
      SUM=0.
      DO 29 K=1,NNF
29    SUM=SUM+B(I,K)*CAP(K,J)
28    A(I,J)=A(I,J)-SUM
      CALL GAUSSJ(A,NNF,MSIZE,IDENT,NNF,MSIZE,IPIV,INDXR,INDXC)
      DO 30 I=1,NNF
      DO 30 J=1,NNF
      SUM=0.
      DO 31 K=1,NNF
31    SUM=SUM-CAP(I,K)*A(K,J)
30    B(I,J)=SUM
      DO 32 I=1,N
      DO 32 J=1,N
      SUM1=0.
      SUM2=0.
      DO 33 K=1,NF
      SUM1=SUM1+A((I-1)*NF+1,K+(J-1)*NF)
33    SUM2=SUM2+B((I-1)*NF+1,K+(J-1)*NF)
      CGEN(I,J)=2.*RW*PI*SUM1+2.*RD*PI*SUM2
32    CONTINUE
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
44    FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=C0(',I3,','I3,')')
43    CONTINUE
      GO TO 997
999   CONTINUE
      WRITE(6,998)
998   FORMAT('INPUT DATA ERROR')
997   CONTINUE
      RW=RW/CMTM
      TD=TD/CMTM
      S=S/CMTM
      WRITE(6,60) N,NF,IREF,RW,TD,ER,S
60    FORMAT(///,'NUMBER OF WIRES=',I3/'NUMBER OF FOURIER COEFFICIENTS='
     1,I3/'REFERENCE WIRE=',I3/'WIRE RADIUS (mils)=',1PE10.3/
     2'DIELECTRIC INSULATION THICKNESS (mils)=',1PE10.3/
     3'DIELECTRIC CONSTANT OF INSULATION=',1PE10.3/
     4'CENTER-TO-CENTER SEPARATION (mils)=',1PE10.3)
      CLOSE(5)
      CLOSE(6)
      STOP
      END
C******SUBROUTINE GAUSSJ FOR INVERTING REAL MATRICES*******************
C     Copyright (C) 1987-1992 Numerical Recipes Software
C     Permission is granted for use only with the programs contained
C     in this book.
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
