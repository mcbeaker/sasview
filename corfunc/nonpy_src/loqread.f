      SUBROUTINE LOQREAD(title1,title2,MaxDim,N,Q,C,E,UNITID,IRET)
C==== 16/1/98 R.K.Heenan, routine to read FISH/COLETTE 1d data file
C==== from i/o channel = UNITID
c
C      DIMENSION Q(MaxDim),C(MaxDim),E(MaxDim)
      INTEGER MAXDIM,N,UNITID,NREAD,IRET
      REAL*4 :: Q(MaxDim)
      REAL*4 :: C(MaxDim)
      REAL*4 :: E(MaxDim)
      character*80 TITLE1,TITLE2
C      DIMENSION NCH(7),FMT(19),RJUNK(2400)
      INTEGER :: NCH(7)
      CHARACTER :: DATAFMT(19)
      REAL*4 :: RJUNK(2400)
C==== DIMENSION IC(4)
C====
C==== DATA       1) TITLE ETC. (A80)  (first 20 are date and time)
C====            2) MORE TITLE and other input e.g. lambda, q calib const.
C====            3) NCH(7)    (7I5)   
C====                 Ntot, NL1,NL2,NMC=CENTRE*10, NR1,NR2,IDC= data type flag
C====                 e.g. file stores Ntot data but only NR1 to NR2 are good
C====            4) NSUM,IC(3)  4I10  monitor counts ( UNUSED )
C====            5) IFLAG, FORTRAN format for rest of data (I2,1X,19A4)
C====            6) DATA depending on IFLAG
C==== IFLAG=1      C( ) Y value only, ASSUME RAW COUNTS E=SQRT(C), 
C====       2      Q( ) C( ) 
C====       3      Q( ) C( ) ERROR( ) 
C
C===== define FORTRAN data streams for data, screen input, screen output
c===== (pass through call      UNITID=4 )
      IS=5
      JS=6
      WRITE(JS,10)
   10 FORMAT(' LOQ DATA input routine, file is read sequentially',/,
     >' ENTER 0 (return) to ignore a set, 1 to KEEP',/,
     >' 91 TO REWIND,  99 TO STOP')
C
C     Return flag; 1=ok, 0=error
      iret=1
C
    5 READ(UNITID,20,END=70)title1
      READ(UNITID,20,END=70)TITLE2
   20 FORMAT(A80)
C
      WRITE(JS,21)title1,TITLE2
   21 FORMAT(1x,a80,/,1X,A80)
      WRITE(JS,30)
   30 FORMAT(1X,'KEEP ? (ANSWER 1) ',$)
      READ(IS,31,ERR=90)I,IEND
   31 FORMAT(2I1)
      IF(IEND.EQ.1)GOTO 100
C
      IF(I.EQ.0 .OR. IEND.GT.1) GOTO 60
      READ(UNITID,51,ERR=80)(NCH(J),J=1,7)
   51 FORMAT(16I5)
C==== skip next record as valus not stored
      READ(UNITID,52)
C====      READ(UNITID,52,ERR=80)(IC(J),J=1,4)
   52 FORMAT(8I10)
      N=NCH(1)
   53 FORMAT(I2,1X,19A4)
      READ(UNITID,53)IFLAG,DATAFMT
      IF(IFLAG.EQ.1)THEN
      READ(UNITID,DATAFMT,ERR=80,END=70)(C(J),J=1,N)
      DO J=1,N
      Q(J)=FLOAT(J)
      END DO
      ELSE IF(IFLAG.EQ.2)THEN
      READ(UNITID,DATAFMT,ERR=80,END=70)(Q(J),C(J),J=1,N)
      ELSE IF(IFLAG.EQ.3)THEN
      READ(UNITID,DATAFMT,ERR=80,END=70)(Q(J),C(J),E(J),J=1,N)
      END IF
C
C==== for IFLAG=1 or 2 calculate errors here, assume COUNTS !
      IF(IFLAG.EQ.1 .OR. IFLAG.EQ.2)THEN
      DO J=1,N
      E(J)=SQRT(ABS(C(J)))
      END DO
      END IF
C

C==== COPY ONLY THE RHS DATA, NCH(5) TO NCH(6), THIS ENABLES ONLY PART
C==== OF INPUT DATA TO BE COPIED INTO GENIE
      N=NCH(6)-NCH(5)+1
      IF(NCH(5).GT.1)THEN
      DO J=1,N
      Q(J)=Q(NCH(5)+J-1)
      C(J)=C(NCH(5)+J-1)
      E(J)=E(NCH(5)+J-1)
      ENDDO
      ENDIF
      GOTO 999
C====
   60 READ(UNITID,51,ERR=80)N
C==== SKIP OVER DATA SET THAT IS NOT WANTED
      READ(UNITID,52)
      READ(UNITID,53)IFLAG,DATAFMT
      N=N*IFLAG
      READ(UNITID,DATAFMT,ERR=80,END=70)(RJUNK(J),J=1,N)
      IF(IEND.NE.0)GOTO 999
      GOTO 5
   70 WRITE(JS,71)
   71 FORMAT(1X,'END OF FILE IN DATIN')
      GOTO 100
   80 WRITE(JS,81)
   81 FORMAT(1X,'FORMAT ERROR IN DATIN')
   90 GOTO 999
  100 REWIND UNITID
      WRITE(JS,101)
  101 FORMAT(1X,'FILE REWOUND')
      GOTO 5
999   RETURN
      END
