C     LAST UPDATE 11/12/92
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      IMPLICIT NONE
C
C Purpose: Given the covariance matrix COVAR of a fit for MFIT of MA
C          total parameters and their ordering LISTA(I) repack the 
C          covariance matrix to the true order of the parameters.
C          Elements associated with fixed parameters will be zero.
C          NCVM is the physical dimension of COVAR. 
C
      INTEGER NCVM,MA,MFIT
C      INTEGER LISTA(MFIT)
C      REAL*8 COVAR(NCVM,NCVM) 
      INTEGER :: LISTA(MFIT)
      REAL*8 :: COVAR(NCVM,NCVM) 
C
C Calls   0:
C Called by: 
C
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C
      REAL*8 SWAP
      INTEGER I,J
C
C---------------------------------------------------------------------
      DO 20 J=1,MA-1
         DO 10 I=J+1,MA
            COVAR(I,J) = 0.0D0
 10      CONTINUE
 20   CONTINUE
      DO 40 I=1,MFIT-1
         DO 30 J=I+1,MFIT
            IF(LISTA(J).GT.LISTA(I))THEN
               COVAR(LISTA(J),LISTA(I)) = COVAR(I,J)
            ELSE
               COVAR(LISTA(I),LISTA(J)) = COVAR(I,J)
            ENDIF
 30      CONTINUE
 40   CONTINUE
      SWAP = COVAR(1,1)
      DO 50 J=1,MA
         COVAR(1,J) = COVAR(J,J)
         COVAR(J,J) = 0.0D0
 50   CONTINUE
      COVAR(LISTA(1),LISTA(1)) = SWAP
      DO 60 J=2,MFIT
         COVAR(LISTA(J),LISTA(J)) = COVAR(1,J)
 60   CONTINUE
      DO 80 J=2,MA
         DO 70 I=1,J-1
            COVAR(I,J) = COVAR(J,I)
 70      CONTINUE
 80   CONTINUE
      RETURN
      END
