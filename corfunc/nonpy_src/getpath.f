      SUBROUTINE getpath(fname,dirname)
C     get the path of fname and return it in dirname
      CHARACTER*80 fname,dirname
      CHARACTER letter

      DO 10 i=80,1,-1
          letter=fname(i:i)
      IF (letter.EQ.'/'.OR.letter.EQ.'\\') THEN
          dirname=fname(1:i)
          goto 20
      ENDIF
 10   CONTINUE
      dirname=""
 20   RETURN
      END

