      SUBROUTINE showprompt(prompts)
C     Displays a prompt
      CHARACTER*80 prompts
      CHARACTER*4 chunk

1000  FORMAT(A,$)
1010  FORMAT(/,1x,A,$)
1020  FORMAT(1x)

      DO i=1,76
        chunk=prompts(i:i+3)
        IF (chunk .EQ. '    ') THEN
          GOTO 10
        ENDIF
      END DO

      i=80
10    WRITE(6,1010)prompts(1:1)
      DO j=2,i-1
        WRITE(6,1000)prompts(j:j)
      END DO

      WRITE(6,1020)

      RETURN
      END



      SUBROUTINE defaultint(prompts,nvalue)
C     Reads in an integer from keyboard displaying prompt,
C     pressing return leaves value unaltered.
      CHARACTER*80 prompts
      CHARACTER*4 chunk
      INTEGER NVALUE

1000  FORMAT(A,$)
1010  FORMAT(/,1x,A,$)
1020  FORMAT(' [',I3,'] : ',/,$)
1030  FORMAT(I5)

      DO i=1,76
        chunk=prompts(i:i+3)
        IF (chunk .EQ. '    ') THEN
          GOTO 10
        ENDIF
      END DO

      i=80
10    WRITE(6,1010)prompts(1:1)
      DO j=2,i-1
        WRITE(6,1000)prompts(j:j)
      END DO

      WRITE(6,1020)nvalue
      READ(*,1030,ERR=20)nvalue2

      IF (nvalue2 .NE. 0) THEN
        nvalue=nvalue2
      ENDIF
      RETURN

20    nvalue=0
      RETURN

      END



      SUBROUTINE defint(prompts,nvalue)
C     Reads in an integer from keyboard displaying prompt,
C     pressing return leaves value unaltered.
      CHARACTER*80 prompts
      CHARACTER*4 chunk
      INTEGER NVALUE
	  
1000  FORMAT(A,$)
1010  FORMAT(1x,A,$)
1020  FORMAT(' [',I3,'] : ',/,$)
1030  FORMAT(I5)

      DO i=1,76
        chunk=prompts(i:i+3)
        IF (chunk .EQ. '    ') THEN
          GOTO 10
        ENDIF
      END DO

      i=80
10    WRITE(6,1010)prompts(1:1)
      DO j=2,i-1
        WRITE(6,1000)prompts(j:j)
      END DO

      WRITE(6,1020)nvalue
      READ(*,1030,ERR=20)nvalue2

      IF (nvalue2 .NE. 0) THEN
        nvalue=nvalue2
      ENDIF
      RETURN

20    nvalue=0
      RETURN

      END



      SUBROUTINE defaultletter(prompts,letter)
C     Reads in one character from keyboard displaying prompt,
C     pressing return leaves value unaltered.
      CHARACTER*80 prompts
      CHARACTER*4 chunk
      CHARACTER*1 letter
      CHARACTER*1 letter2

1000  FORMAT(A,$)
1010  FORMAT(/,1x,A,$)
1020  FORMAT(' [',A,'] : ',/,$)
1030  FORMAT(A1)

      DO i=1,76
        chunk=prompts(i:i+3)
        IF (chunk .EQ. '    ') THEN
          GOTO 10
        ENDIF
      END DO

      i=80
10    WRITE(6,1010)prompts(1:1)
      DO j=2,i-1
        WRITE(6,1000)prompts(j:j)
      END DO

      WRITE(6,1020)letter
      READ(*,1030,ERR=10)letter2

      IF (letter2 .NE. '') THEN
        letter=letter2
      ENDIF

      RETURN
      END



      SUBROUTINE defaultname(prompts,name1)
C     Reads in a filename etc. from keyboard displaying prompt,
C     pressing return leaves default
      CHARACTER*80 prompts
      CHARACTER*4 chunk
      CHARACTER*80 name1,name2

1000  FORMAT(A,$)
1010  FORMAT(/,1x,A,$)
1020  FORMAT('] : ',/,$)
1030  FORMAT(A80)
1040  FORMAT(' [',$)

      DO i=1,76
        chunk=prompts(i:i+3)
        IF (chunk .EQ. '    ') THEN
          GOTO 10
        ENDIF
      END DO

      i=80
10    WRITE(6,1010)prompts(1:1)
      DO j=2,i-1
        WRITE(6,1000)prompts(j:j)
      END DO

      WRITE(6,1040)

      DO i=1,77
        chunk=name1(i:i+3)
        IF (chunk .EQ. '    ') THEN
          GOTO 20
        ENDIF
      END DO

20    DO j=1,i-1
        WRITE(6,1000)name1(j:j)
      END DO

      WRITE(6,1020)
      READ(*,1030,ERR=10)name2

      IF (name2 .NE. '') THEN
        name1=name2
      ENDIF

      RETURN
      END



      SUBROUTINE showtitle(title)
C     Displays a title.
      CHARACTER*40 title
      CHARACTER*42 underline
      CHARACTER*4 chunk

1000  FORMAT(/,1x,A,$)
1010  FORMAT(A,$)
1020  FORMAT(/,A,$)
1030  FORMAT(1x)

      underline='------------------------------------------'

      DO i=1,36
        chunk=title(i:i+3)
        IF (chunk .EQ. '    ') THEN
          GOTO 10
        ENDIF
      END DO

      i=40
10    WRITE(6,1000)title(1:1)
      DO j=2,i-1
        WRITE(6,1010)title(j:j)
      END DO
      WRITE(6,1020)underline(1:1)
      DO j=2,i+1
        WRITE(6,1010)underline(j:j)
      END DO
      WRITE(6,1030)

      RETURN
      END



      SUBROUTINE defaultreal(prompts,realvalue)
C      SUBROUTINE defaultreal(prompts,value)
C     Reads in a real no. from keyboard displaying prompt,
C     pressing return leaves value unaltered.
C     Involves file output!
      CHARACTER*80 prompts,text
      CHARACTER*4 chunk
      REAL REALVALUE

1000  FORMAT(A,$)
1010  FORMAT(/,1x,A,$)
1020  FORMAT(' [',E12.6,'] : ',/,$)
1030  FORMAT(A80)

      DO i=1,76
        chunk=prompts(i:i+3)
        IF (chunk .EQ. '    ') THEN
          GOTO 10
        ENDIF
      END DO

      i=80
10    WRITE(6,1010)prompts(1:1)
      DO j=2,i-1
        WRITE(6,1000)prompts(j:j)
      END DO

C      WRITE(6,1020)value
      WRITE(6,1020)realvalue
      READ(*,1030,ERR=20)text

      IF(text.NE.' ') THEN
C        READ(text,*)value
        READ(text,*)realvalue
      ENDIF

C      IF (text .NE. '') THEN
C        OPEN(UNIT=8,FILE='dummy.num',ERR=20,STATUS='unknown')
C        WRITE(8,1030,ERR=20)text
C        CLOSE(8)
C        OPEN(UNIT=8,FILE='dummy.num',ERR=20,STATUS='old')
C        READ(8,*,ERR=20)value2
C        CLOSE(8)
C        value=value2
C      ENDIF

20    RETURN
      END



      SUBROUTINE changeotok(name1,name2,nerr)
C     Alters an otoko header name to a data filename
C     eg. X23000.ABS becomes X23001.ABS
C     nerr = 0 if everything goes smoothly.
C     nerr = 1 if change cannot be made.

C      CHARACTER*40 chunk*4
      CHARACTER*4 chunk
      CHARACTER*80 name1,name2
      INTEGER NERR

      nerr=1
      DO i=1,77
        chunk=name1(i:i+3)
        IF (chunk .EQ. '000.') THEN
          GOTO 10
        ENDIF
      END DO

      RETURN

10    nerr=0
      name2=name1(1:i+1)//'1.'//name1(i+4:80)

      RETURN
      END



      SUBROUTINE swapexten(filename,exten)
C     Changes a file extension eg. .XAX to .QAX

      CHARACTER*80 filename
      CHARACTER*3 exten
      INTEGER i

 1000  FORMAT(/,1x,'Error in subroutine swapexten: fatal...')

      DO i=1,80
        IF(filename(i:i).EQ.'.')THEN
          filename=filename(1:i)//exten
          RETURN
        ENDIF
      END DO

      WRITE(6,1000)
      STOP

      RETURN
      END
