C     LAST UPDATE 20/06/96
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SYMINV(A,W,V,N,NP)
      IMPLICIT NONE
C
C Purpose: Inverts a symmetric matrix specified by the 
C          arrays W,V as returned by TRED2,TQLI. N is the logical
C          dimension of A, the output matrix (square and symmetric). 
C          NP are the physical dimensions of A. 
C
      INTEGER N,NP
C      DOUBLE PRECISION A(NP,NP),V(NP,NP),W(NP)
      REAL*8 :: A(NP,NP)
      REAL*8 :: V(NP,NP)
      REAL*8 :: W(NP)
C
C Calls   0:
C Called by: 
C
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C
C      DOUBLE PRECISION S,WTMP,VTMP
      REAL*8 S,WTMP,VTMP
      INTEGER I,J,K
C
C---------------------------------------------------------------------
C
C========Loop over eigenvalues
C
      DO 50 K=1,N
         IF(W(K).NE.0.0D0)THEN
            WTMP = 1.0D0/W(K)
C
C========Accumulate diagonal elements
C
            DO 20 I=1,N
               A(I,I) = V(I,K)*WTMP*V(I,K)
 20         CONTINUE
C
C========Accumulate off-diagonal elements using symmetry
C
            DO 40 I=1,N-1
               VTMP = V(I,K)
               DO 30 J=I+1,N
                  S = VTMP*WTMP*V(J,K)
                  A(I,J) = S
                  A(J,I) = S
 30            CONTINUE
 40         CONTINUE
         ENDIF
 50   CONTINUE
      RETURN
      END
