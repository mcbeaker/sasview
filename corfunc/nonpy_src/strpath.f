      SUBROUTINE strippath(fname1,fname2)
C     get the path of fname and return it in dirname
      CHARACTER*80 fname1,fname2
      CHARACTER letter

      DO 10 i=80,1,-1
          letter=fname1(i:i)
      IF (letter.EQ.'/'.OR.letter.EQ.'\\') THEN
          fname2=fname1(i+1:i+11)
          goto 20
      ENDIF
 10   CONTINUE
      fname2=fname1
 20   RETURN
      END

