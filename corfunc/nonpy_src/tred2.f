C     LAST UPDATE 07/05/93
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE TRED2(A,N,NP,D,E)
      IMPLICIT NONE
C
C Purpose: Householder reduction of a real, symmetric, N by N matrix A,
C          stored in an NP by NP physical array. On output, A is
C          replaced by the orthogonal matrix Q effecting the transformation  
C          D returns the diagonal elements of the tridiagonal matrix and
C          E the off-diagonal elements, with E(1) = 0. 
C
      INTEGER N,NP 
C      REAL*8 A(NP,NP),D(NP),E(NP)
      REAL*8 :: A(NP,NP)
      REAL*8 :: D(NP)
      REAL*8 :: E(NP)
C
C Calls   0:
C Called by: 
C
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C
      REAL*8 H,SCALEFAC,F,G,HH  
      INTEGER I,J,K,L  
C
C---------------------------------------------------------------------
      DO 80 I=N,2,-1
         L = I - 1
         H = 0.0
         SCALEFAC = 0.0
         IF(L.GT.1)THEN 
            DO 10 K=1,L
               SCALEFAC = SCALEFAC + ABS(A(I,K))
 10         CONTINUE
            IF(SCALEFAC.EQ.0.0)THEN
               E(I) = A(I,L)
            ELSE
               DO 20 K=1,L
                  A(I,K) = A(I,K)/SCALEFAC
                  H = H + A(I,K)**2
 20            CONTINUE
               F = A(I,L)
               G = -SIGN(SQRT(H),F)
               E(I) = SCALEFAC*G
               H = H - F*G
               A(I,L) = F - G
               F = 0.0
               DO 50 J=1,L
                  A(J,I) = A(I,J)/H
                  G = 0.0
                  DO 30 K=1,J
                     G = G + A(J,K)*A(I,K)
 30               CONTINUE
                  DO 40 K=J+1,L
                     G = G + A(K,J)*A(I,K)
 40               CONTINUE
                  E(J) = G/H
                  F = F + E(J)*A(I,J)
 50            CONTINUE
               HH = F/(H+H)
               DO 70 J=1,L
                  F = A(I,J)
                  G = E(J) - HH*F
                  E(J) = G
                  DO 60 K=1,J
                     A(J,K) = A(J,K) - F*E(K) - G*A(I,K)
 60               CONTINUE
 70            CONTINUE
            ENDIF
         ELSE
            E(I) = A(I,L)
         ENDIF
         D(I) = H
 80   CONTINUE
      D(1) = 0.0
      E(1) = 0.0
      DO 130 I=1,N
         L = I - 1
         IF(D(I).NE.0.0)THEN
            DO 110 J=1,L
               G = 0.0
               DO 90 K=1,L
                  G = G + A(I,K)*A(K,J)
 90            CONTINUE
               DO 100 K=1,L
                  A(K,J) = A(K,J) - G*A(K,I)
 100           CONTINUE
 110        CONTINUE
         ENDIF
         D(I) = A(I,I)
         A(I,I) = 1.0
         DO 120 J=1,L
            A(I,J) = 0.0
            A(J,I) = 0.0
 120     CONTINUE
 130  CONTINUE 
      RETURN
C
      END 
