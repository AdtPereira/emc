C******PROGRAM SPICEMTL.FOR********************************************
C
C  Generates a SPICE subcircuit model of a Multiconductor Transmission
C  Line (MTL) by diagonalizing the NxN per-unit-length inductance, L,
C  and capacitance, C, matrices of the (N+1)-conductor MTL.  The total
C  number of conductors, (N+1), and the total line length are input data
C  SPICEMTL.IN.  The per-unit-length parameters in the upper triangle
C  of L and C are input data from PUL.DAT.
C
C  The SPICE subcircuit is stored in SPICEMTL.OUT as output with the
C  following port node numbering convention (node 0 is the global
C  ground or reference node and refers to the (N+1)st reference
C  conductor (lossless) of the MTL):
C
C   Node Numbering: 101-10N are input nodes at left end of MTL.
C                   201-20N are output nodes at the right end of MTL.
C
C  NOTE: Parameter MSIZE=9 dimensions arrays to handle problems of
C        N<10 conductors!  Array space can be minimized by setting
C        MSIZE=N.
C
C          Written by: Clayton R. Paul
C                      Department of Electrical Engineering
C                      University of Kentucky
C                      Lexington, KY 40506
C
C**********************************************************************
      PARAMETER(MSIZE=2)
      INTEGER LP(MSIZE),RP(MSIZE),LPI(MSIZE),RPI(MSIZE)
      INTEGER LPII(MSIZE),RPII(MSIZE)
      REAL LEN,C(MSIZE,MSIZE),L(MSIZE,MSIZE),TV(MSIZE,MSIZE)
      REAL TEMP(MSIZE,MSIZE),U(MSIZE,MSIZE),S(MSIZE,MSIZE)
      REAL TI(MSIZE,MSIZE)
      OPEN(UNIT=5,FILE='SPICEMTL.IN',STATUS='OLD')
      OPEN(UNIT=6,FILE='SPICEMTL.OUT',STATUS='UNKNOWN')
      OPEN(UNIT=7,FILE='PUL.DAT',STATUS='OLD')
C
C*****READ INPUT DATA FROM SPICEMTL.IN
C
      READ(5,*) NP1
      N=NP1-1
      READ(5,*) LEN
      IF(N.GT.MSIZE) GO TO 16
      DO 1 I=1,N
      DO 1 J=I,N
      READ(7,*) M,K,L(M,K)
      L(J,I)=L(I,J)
1     CONTINUE
      DO 2 I=1,N
      DO 2 J=I,N
      READ(7,*) M,K,C(M,K)
      C(J,I)=C(I,J)
2     CONTINUE
      WRITE(6,130) NP1,LEN
130   FORMAT('*SUBCIRCUIT MODEL OF A MULTICONDUCTOR TRANSMISSION LINE*'/
     1'*      NUMBER OF CONDUCTORS=',I3/
     2'*      TOTAL LINE LENGTH (METERS)=',1PE12.5)
      DO 20 I=1,N
      DO 20 J=I,N
      WRITE(6,131) I,J,L(I,J)
131   FORMAT('*      L(',I3,',',I3,')=',1PE12.5)
20    CONTINUE
      DO 21 I=1,N
      DO 21 J=I,N
      WRITE(6,132) I,J,C(I,J)
132   FORMAT('*      C(',I3,',',I3,')=',1PE12.5)
21    CONTINUE
C
C*****GENERATE ARRAYS HOLDING MODEL NODE NUMBERS
C     LP and RP HOLD EXTERNAL NODE NUMBERS AT LEFT AND RIGHT
C     LPI AND RPI HOLD INTERNAL NODE NUMBERS AFTER ZERO-VOLT SOURCES
C     LPII AND RPII HOLD INTERNAL NODE NUMBERS AT ENDS OF MODAL LINES
C
      DO 3 I=1,N
      LP(I)=100+I
      RP(I)=200+I
      LPI(I)=300+I
      RPI(I)=400+I
      LPII(I)=500+I
3     RPII(I)=600+I
C
C*****DIAGONALIZE L AND C TO GET MODE TRANSFORMATIONS TV, TI-1=TVtran
C     SO THAT TV-1xLxTI=lm  AND  TI-1xCxTV=cm
C
      CALL DIAG(C,L,TV,TI,N,MSIZE,TEMP,U,S)
C
C*****PRINT SUBCKT CARDS TO SPICEMTL.OUT
C     DO LEFT END FIRST
C
      WRITE(6,103)
103   FORMAT('.SUBCKT MTL')
      DO 4 I=1,N
      WRITE(6,104) LP(I)
104   FORMAT('+',1X,I3)
4     CONTINUE
      DO 5 I=1,N
      WRITE(6,105) RP(I)
105   FORMAT('+',1X,I3)
5     CONTINUE
      DO 6 I=1,N
      WRITE(6,106) LP(I),LP(I),LPI(I)
106   FORMAT('V',I3,2X,I3,2X,I3)
      WRITE(6,107) LP(I),LPI(I),N
107   FORMAT('EC',I3,2X,I3,2X,'0',2X,'POLY(',I2,')')
      DO 7 J=1,N
      WRITE(6,108) LPII(J)
108   FORMAT('+ (',I3,',0)')
7     CONTINUE
      WRITE(6,109)
109   FORMAT('+0')
      DO 8 J=1,N
      WRITE(6,110) TV(I,J)
110   FORMAT('+',2X,1PE13.6)
8     CONTINUE
      WRITE(6,111) LP(I),LPII(I),N
111   FORMAT('FC',I3,2X,'0',2X,I3,2X,'POLY(',I2,')')
      DO 9 J=1,N
      WRITE(6,112) LP(J)
112   FORMAT('+ V',I3)
9     CONTINUE
      WRITE(6,113)
113   FORMAT('+0')
      DO 10 J=1,N
      WRITE(6,114) TV(J,I)
114   FORMAT('+',2X,1PE13.6)
10    CONTINUE
C
C*****NOW DO RIGHT END
C
      WRITE(6,115) RP(I),RP(I),RPI(I)
115   FORMAT('V',I3,2X,I3,2X,I3)
      WRITE(6,116) RP(I),RPI(I),N
116   FORMAT('EC',I3,2X,I3,2X,'0',2X,'POLY(',I2,')')
      DO 11 J=1,N
      WRITE(6,117) RPII(J)
117   FORMAT('+ (',I3,',0)')
11    CONTINUE
      WRITE(6,118)
118   FORMAT('+0')
      DO 12 J=1,N
      WRITE(6,119) TV(I,J)
119   FORMAT('+',2X,1PE13.6)
12    CONTINUE
      WRITE(6,120) RP(I),RPII(I),N
120   FORMAT('FC',I3,2X,'0',2X,I3,2X,'POLY(',I2,')')
      DO 13 J=1,N
      WRITE(6,121) RP(J)
121   FORMAT('+ V',I3)
13    CONTINUE
      WRITE(6,122)
122   FORMAT('+0')
      DO 14 J=1,N
      WRITE(6,123) TV(J,I)
123   FORMAT('+',2X,1PE13.6)
14    CONTINUE
6     CONTINUE
C
C*****PRINT MODAL LINE DESCRIPTIONS
C
      DO 15 I=1,N
      WRITE(6,124)LP(I),LPII(I),RPII(I),C(I,I),LEN/L(I,I)
124   FORMAT('T',I3,2X,I3,2X,'0',2X,I3,2X,'0',2X,'Z0=',1PE13.6,2X,
     1'TD=',1PE13.6)
15    CONTINUE
      WRITE(6,125)
125   FORMAT('.ENDS',5X,'MTL')
16    CONTINUE
      CLOSE(5)
      CLOSE(6)
      CLOSE(7)
      STOP
      END
C**********************************************************************
      SUBROUTINE DIAG(C,L,TV,TI,N,MSIZE,TEMP,U,S)
C**********************************************************************
C        Written by: Clayton R. Paul
C                    Department of Electrical Engineering
C                    University of Kentucky
C                    Lexington, KY 40506
C     To Diagonalize Two NxN Real, Symmetric Matrices, C and L, with
C     C Positive Definite as TV-1xLxTI=lm(diag) and TI-1xCxTV=cm(diag)
C     and TI-1=TVtran.
C
C     C and L are DESTROYED.  On Return, C Contains the MODE
C     Characteristic Impedances, ZCmi=sqrt(lmi/cmi) and L Contains
C     the MODE Velocities of Propagation, vmi=1/sqrt(lmicmi).  These
C     are Stored on the Main Diagonal of these Arrays.
C
C     Arrays MUST Be Dimensioned MSIZE in Calling Program
C
C     Subroutines Used: JACOBI
C**********************************************************************
      REAL C(MSIZE,MSIZE),L(MSIZE,MSIZE),TEMP(MSIZE,MSIZE)
      REAL U(MSIZE,MSIZE),TV(MSIZE,MSIZE),S(MSIZE,MSIZE),TI(MSIZE,MSIZE)
C
C*****Diagonalize C as UtranxCxU=D
C
      CALL JACOBI(C,U,N,MSIZE)
      DO 1 I=1,N
1     C(I,I)=SQRT(C(I,I))
C
C*****Form D1/2xUtranxLxUxD1/2
C
      DO 2 I=1,N
      DO 2 J=1,N
      SUM=0.0
      DO 3 K=1,N
3     SUM=SUM+L(I,K)*U(K,J)*C(J,J)
2     TEMP(I,J)=SUM
C
      DO 4 I=1,N
      DO 4 J=1,N
      SUM=0.0
      DO 5 K=1,N
5     SUM=SUM+C(I,I)*U(K,I)*TEMP(K,J)
4     L(I,J)=SUM
C
C*****Diagonalize D1/2xUtranxLxUxD1/2 with Orthogonal S
C
      CALL JACOBI(L,S,N,MSIZE)
C
C*****Form T=UxD1/2xS
C
      DO 6 I=1,N
      DO 6 J=1,N
      SUM=0.
      DO 7 K=1,N
7     SUM=SUM+U(I,K)*C(K,K)*S(K,J)
6     TV(I,J)=SUM
C
C*****Normalize T as Tnorm=TxAlpha
C
      DO 8 I=1,N
      SUM=0.0
      DO 9 J=1,N
9     SUM=SUM+TV(J,I)*TV(J,I)
8     TEMP(I,I)=1.0/SQRT(SUM)
C
C*****Form TI=UxD1/2xSxAlpha
C
      DO 10 I=1,N
      DO 10 J=1,N
10    TI(I,J)=TV(I,J)*TEMP(J,J)
C
C*****Form TV=UxD-1/2xSxAlpha-1
      DO 11 I=1,N
      DO 11 J=1,N
      SUM=0.0
      DO 12 K=1,N
12    SUM=SUM+U(I,K)*S(K,J)/(C(K,K)*TEMP(J,J))
11    TV(I,J)=SUM
C
C*****Form ZCmi=sqrt(lmi/cmi) and vmi=1/sqrt(lmicmi)
C
      DO 13 I=1,N
      C(I,I)=TEMP(I,I)*TEMP(I,I)*SQRT(L(I,I))
13    L(I,I)=1.0/SQRT(L(I,I))
C
      RETURN
      END
C
C
C
      SUBROUTINE JACOBI(A,T,N,MSIZE)
C**********************************************************************
C        Written by: Clayton R. Paul
C                    Department of Electrical Engineering
C                    University of Kentucky
C                    Lexington, KY 40506
C
C     To Diagonalize a NxN Real, Symmetric Matrix, A, With an Orthogonal
C     Transformation, T, so that TtranxAxT=Diagonal.
C            A is NxN Real, Symmetric Input Matrix (DESTROYED)
C            Eigenvalues of A are returned on main diagonal of A
C            T is NxN Orthogonal Matrix (TtranxT=Identity)
C     A and T must be dimensioned as below in calling program as MSIZE.
C     N is actual dimension of A and T.
C     N must be less than or equal to MSIZE.
C
C     METHOD USED: JACOBI Method with Threshold Method for Convergence
C
C     REFERENCES:(1)A. Ralston,"A First Course in Numerical Analysis,"
C                   McGraw-Hill, 1965, Chapter 10.
C                (2)R.W. Hornbeck,"Numerical Methods,"
C                   Quantum, 1975, Chapter 10.
C**********************************************************************
      PARAMETER(EPS=1.E-6)
      REAL A(MSIZE,MSIZE),T(MSIZE,MSIZE),INORM
C
C*****Initialize T to Identity Matrix
C
      DO 1 I=1,N
      DO 2 J=1,N
2     T(I,J)=0.0
1     T(I,I)=1.0
C
C*****Compute Initial and Final Norms from Off-Diagonal Elements of A
C
      INORM=0.0
      DO 3 I=1,N-1
      DO 3 J=I+1,N
3     INORM=INORM+A(I,J)*A(I,J)
      THRS=SQRT(2.*INORM)/FLOAT(N)
      FNORM=THRS*EPS
C
C*****Annihilate Off-Diagonal Terms Greater Than Threshold
C*****IR is Row Location and IC is Column Location of Element that
C*****      Exceeds Threshold
C
310   CONTINUE
      DO 4 I=1,N-1
      DO 4 J=I+1,N
      IF(THRS-ABS(A(I,J))) 100,110,110
100   IR=I
      IC=J
      IF(A(IR,IR).EQ.A(IC,IC)) GO TO 7
      ALPH=0.5*(A(IR,IR)-A(IC,IC))
      ABALPH=ABS(ALPH)
      BETA=SQRT(A(IR,IC)*A(IR,IC)+ALPH*ALPH)
      C=SQRT(0.5+0.5*ABALPH/BETA)
      S=-0.5*ALPH*A(IR,IC)/(BETA*ABALPH*C)
      GO TO 8
7     CONTINUE
      C=1./SQRT(2.)
      S=-1./SQRT(2.)
8     CONTINUE
      XIR=A(IR,IR)
      XIC=A(IC,IC)
      XIRIC=A(IR,IC)
      DO 5 K=1,N
      TEMP=C*A(IR,K)-S*A(IC,K)
      A(IC,K)=S*A(IR,K)+C*A(IC,K)
      A(IR,K)=TEMP
      TEMP=C*A(K,IR)-S*A(K,IC)
      A(K,IC)=S*A(K,IR)+C*A(K,IC)
5     A(K,IR)=TEMP
      A(IR,IC)=0.
      A(IC,IR)=0.
      A(IR,IR)=C*C*XIR+S*S*XIC-2.*C*S*XIRIC
      A(IC,IC)=C*C*XIC+S*S*XIR+2.*C*S*XIRIC
C
C*****Compute T Matrix
C
      DO 6 K=1,N
      XIR=T(K,IR)
      XIC=T(K,IC)
      T(K,IR)=C*XIR-S*XIC
6     T(K,IC)=S*XIR+C*XIC
110   CONTINUE
4     CONTINUE
C
C*****Compute New Threshold
C
      THRS=THRS/FLOAT(N)
      IF(THRS-FNORM) 300,300,310
300   CONTINUE
      RETURN
      END
