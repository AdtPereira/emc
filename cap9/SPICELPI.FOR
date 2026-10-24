C******PROGRAM SPICELPI.FOR********************************************
C
C  Generates a SPICE subcircuit model of a Multiconductor Transmission
C  Line (MTL) as an approximate, lumped-PI equivalent circuit.  Model
C  is valid only for frequencies where the line is electrically short.
C  The number of conductors, (N+1), and the total line length are input data
C  SPICELPI.IN.  The per-unit-length parameters in the upper triangle
C  of L and C are input data from PUL.DAT.
C
C  The SPICE subcircuit is stored in SPICELPI.OUT as output with the
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
      INTEGER LP(MSIZE),RP(MSIZE)
      REAL LEN,LTOTAL,C(MSIZE,MSIZE),L(MSIZE,MSIZE)
      OPEN(UNIT=5,FILE='SPICELPI.IN',STATUS='OLD')
      OPEN(UNIT=6,FILE='SPICELPI.OUT',STATUS='UNKNOWN')
      OPEN(UNIT=7,FILE='PUL.DAT',STATUS='OLD')
C
C*****READ INPUT DATA FROM SPICELPI.IN AND PUL.DAT
C
      READ(5,*) NP1
      N=NP1-1
      READ(5,*) LEN
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
     1'*LUMPED-PI APPROXIMATE STRUCTURE*'/
     2'*      TOTAL NUMBER OF CONDUCTORS=',I3/
     3'*      TOTAL LINE LENGTH (METERS)=',1PE12.5)
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
C
      DO 3 I=1,N
      LP(I)=100+I
3     RP(I)=200+I
C
C*****PRINT SUBCKT CARDS TO SPICELPI.OUT
C
      WRITE(6,103)
103   FORMAT('.SUBCKT LPI')
      DO 4 I=1,N
      WRITE(6,104) LP(I)
104   FORMAT('+',1X,I3)
4     CONTINUE
      DO 5 I=1,N
      WRITE(6,105) RP(I)
105   FORMAT('+',1X,I3)
5     CONTINUE
      DO 8 I=1,N
      SUM=0.
      DO 9 J=1,N
9     SUM=SUM+C(I,J)
8     C(I,I)=SUM
      DO 6 I=1,N
      LTOTAL=L(I,I)*LEN
      WRITE(6,106) LP(I),LP(I),RP(I),LTOTAL
106   FORMAT('L',I3,2X,I3,2X,I3,3X,1PE12.5)
      CTOTAL=C(I,I)*LEN/2.
      WRITE(6,107) LP(I),LP(I),CTOTAL
107   FORMAT('C',I3,2X,I3,2X,'0',3X,1PE12.5)
      WRITE(6,108) RP(I),RP(I),CTOTAL
108   FORMAT('C',I3,2X,I3,2X,'0',3X,1PE12.5)
6     CONTINUE
      INDEX=100
      JNDEX=200
      DO 7 I=1,N
      DO 7 J=I+1,N
      INDEX=INDEX+1
      JNDEX=JNDEX+1
      XK=L(I,J)/SQRT(L(I,I)*L(J,J))
      INDLI=100+I
      INDLJ=100+J
      WRITE(6,109) INDEX,INDLI,INDLJ,XK
109   FORMAT('K',I3,2X,'L',I3,2X,'L',I3,3X,1PE12.5)
      CTOTAL=-C(I,J)*LEN/2.
      WRITE(6,110) INDEX,LP(I),LP(J),CTOTAL
110   FORMAT('CM',I3,2X,I3,2X,I3,3X,1PE12.5)
      WRITE(6,111) JNDEX,RP(I),RP(J),CTOTAL
111   FORMAT('CM',I3,2X,I3,2X,I3,3X,1PE12.5)
7     CONTINUE
      WRITE(6,112)
112   FORMAT('.ENDS     LPI')
      CLOSE(5)
      CLOSE(6)
      STOP
      END
