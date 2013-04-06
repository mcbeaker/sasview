C     LAST UPDATE 07/05/93
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE TQLI(D,E,N,NP,Z)
      IMPLICIT NONE
C
C Purpose: QL algorithm with implicit shifts, to determine the
C          eigenvalues and eigenvectors of a real symmetric tridiagonal
C          matrix, or of a real symmetric matrix previously produced by
C          TRED2. D is a vector of length NP. On input, its first N
C          elements are the diagonal elements of the tridiagonal matrix.
C          On output, it returns the eigenvalues. The vector E inputs
C          the subdiagonal elements of the tridiagonal matrix, with E(1)
C          arbitrary. On output, E is destroyed. If the eigenvectors of
C          a tridiagonal matrix are desired, the matrix Z is input as
C          the identity matrix. If the eigenvectors of a matrix that has
C          been reduced by TRED2 are required then Z is input as the
C          matric output by TRED2. In either case, the Kth column of Z
C          returns the normalized eigenvector corresponding to D(K)
C
      INTEGER N,NP 
C      REAL*8 D(NP),E(NP),Z(NP,NP)
      REAL*8 :: D(NP)
      REAL*8 :: E(NP)
      REAL*8 :: Z(NP,NP)
C
C Calls   0:
C Called by: 
C
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C
      REAL*8 C,R,S,G,DD,P,B,F  
      INTEGER I,K,L,M,ITER 
C
C---------------------------------------------------------------------
      DO 10 I=2,N
         E(I-1) = E(I)
 10   CONTINUE
      E(N) = 0.0
      DO 50 L=1,N
         ITER = 0
 1       DO 20 M=L,N-1
            DD = ABS(D(M)) + ABS(D(M+1))
            IF(ABS(E(M))+DD.EQ.DD)GOTO 2 
 20      CONTINUE
         M = N
 2       IF(M.NE.L)THEN
            IF(ITER.EQ.50)WRITE(6,1000)
            ITER = ITER + 1
            G = (D(L+1)-D(L))/(2.0*E(L))
            R = SQRT(G**2+1.0)
            G = D(M) - D(L) + E(L)/(G+SIGN(R,G))
            S = 1.0
            C = 1.0
            P = 0.0
            DO 40 I=M-1,L,-1
               F = S*E(I)
               B = C*E(I)
               IF(ABS(F).GE.ABS(G))THEN
                  C = G/F
                  R = SQRT(C**2+1.0)
                  E(I+1) = F*R
                  S = 1.0/R
                  C = C*S
               ELSE
                  S = F/G
                  R = SQRT(S**2+1.0)
                  E(I+1) = G*R
                  C = 1.0/R
                  S = S*C
               ENDIF
               G = D(I+1) - P
               R = (D(I)-G)*S + 2.0*C*B
               P = S*R
               D(I+1) = G + P
               G = C*R - B
               DO 30 K=1,N
                  F = Z(K,I+1)
                  Z(K,I+1) = S*Z(K,I) + C*F
                  Z(K,I) = C*Z(K,I) - S*F
 30            CONTINUE
 40         CONTINUE
            D(L) = D(L) - P
            E(L) = G
            E(M) = 0.0
            IF(ITER.LT.50)GOTO 1
         ENDIF
 50   CONTINUE
      RETURN
 1000 FORMAT(1X,'TQLI: Warning - no convergence in 50 iterations')
C
      END 
