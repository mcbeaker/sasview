C     LAST UPDATE 06/06/96
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE TBKSB(W,V,N,NP,B,X)
      IMPLICIT NONE
C
C Purpose: Solves A.X = B for a vector X, where A is specified by the 
C          arrays W,V as returned by TRED2,TQLI. N is the logical
C          dimension of A (originally square and symmetric). 
C          NP are the physical dimensions of A. B is the input right
C          hand side. X is the output solution vector. No input
C          quantities are destroyed so the routine may be called 
C          sequentially with different B's.
C
      INTEGER N,NP
C      DOUBLE PRECISION V(NP,NP),W(NP),B(NP),X(NP)
      REAL*8 :: V(NP,NP)
      REAL*8 :: W(NP)
      REAL*8 :: B(NP)
      REAL*8 :: X(NP)
C
C Calls   0:
C Called by: 
C
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C
C      DOUBLE PRECISION S,WTMP,VTMP,XTMP,BTMP
      REAL*8 S,WTMP,VTMP,XTMP,BTMP
      INTEGER I,J,K
C
C---------------------------------------------------------------------
C
C========Initialize X
C
      DO 10 I=1,N
         X(I) = 0.0D0
 10   CONTINUE
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
               X(I) = X(I) + V(I,K)*WTMP*V(I,K)*B(I)
 20         CONTINUE
C
C========Accumulate off-diagonal elements using symmetry
C
            DO 40 I=1,N-1
               XTMP = X(I)
               BTMP = B(I)
               VTMP = V(I,K)
               DO 30 J=I+1,N
                  S = VTMP*WTMP*V(J,K)
                  XTMP = XTMP + S*B(J)
                  X(J) = X(J) + S*BTMP
 30            CONTINUE
               X(I) = XTMP
 40         CONTINUE
         ENDIF
 50   CONTINUE
      RETURN
      END
