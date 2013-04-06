      SUBROUTINE addstrings(dirname,filename,outname)
C     adds a directory path to a filename
      CHARACTER*80 outname,filename,dirname

      DO 10 i=1,80
          IF(dirname(i:i).EQ.' ')THEN
              outname=dirname(1:i-1)//filename
              RETURN
          ENDIF
 10   CONTINUE
      outname=dirname//filename
      RETURN
      END

