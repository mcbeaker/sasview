C     LAST UPDATE 27/07/95
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE MRQMIN(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,ALPHA,NCA,
     &                  CHISQ,FUNCS,ALAMDA)
      IMPLICIT NONE
C
C Purpose: Non-linear least-squares by Levenberg-Marquardt method.
C          Attempts to reduce the value of CHISQ of a fit between a set
C          of NDATA points X(I),Y(I) with individual standard deviations
C          SIG(I) and a non-linear function dependent on MA coefficients
C          A. The array LISTA numbers the parameters A such that the
C          first MFIT elements correspond to values actually being
C          adjusted; the remaining MA-MFIT parameters are held fixed at
C          their input values. The program returns current best-fit
C          values for the  MA fit parameters A and CHISQ. The arrays
C          COVAR(NCA,NCA), ALPHA(NCA,NCA) with physical dimension NCA 
C          (>= MFIT) are used as working space during most iterations. 
C          Supply a subroutine FUNCS(X,A,YFIT,DYDA,MA) that evaluates
C          the fitting function YFIT and its derivatives DYDA with
C          respect to the fitting parameters A at X. On the first call
C          provide an initial guess for the parameters A and set ALAMDA
C          < 0 for initialization (which then sets ALAMDA=0.001). If a 
C          step succeeds CHISQ becomes smaller and ALAMDA decreases by a
C          factor of 10. You must call this routine repeatedly until
C          convergence is achieved. Then make one final call with
C          ALAMDA=0.0, so that COVAR(I,J) returns the covariance matrix.
C
      INTEGER NDATA,MA,NCA,MFIT
C      REAL X(NDATA),Y(NDATA),SIG(NDATA)
      REAL*4 :: X(NDATA)
      REAL*4 :: Y(NDATA)
      REAL*4 :: SIG(NDATA)
C      DOUBLE PRECISION A(NCA)
      REAL*8 :: A(NCA)
C      INTEGER LISTA(NCA)
      INTEGER :: LISTA(NCA)
C      DOUBLE PRECISION COVAR(NCA,NCA),ALPHA(NCA,NCA)
      REAL*8 :: COVAR(NCA,NCA)
      REAL*8 :: ALPHA(NCA,NCA)
      REAL CHISQ,ALAMDA
      EXTERNAL FUNCS
C
C Calls   6: TRED2 , TQLI , TBKSB , SYMINV , COVSRT , MRQCOF
C Called by: 
C
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C
      INTEGER NMAX 
C      DOUBLE PRECISION TINY,RMIN,FACTOR,ALAMIN,ALAMAX,DELMIN,DELMAX
      REAL*8 TEENSY,RMIN,FACTOR,ALAMIN,ALAMAX,DELMIN,DELMAX
      PARAMETER (NMAX=50,TEENSY=1.0D-30,RMIN=1.0D-06,FACTOR=10.0D0)
      PARAMETER (ALAMIN=1.0D-04,ALAMAX=1.0D+10,DELMIN=1.0D-30,
     &           DELMAX=1.0D+30)
C      DOUBLE PRECISION ATRY(NMAX),BETA(NMAX),DA(NMAX),D(NMAX),E(NMAX),
C     &                 PSCALE(NMAX)
      REAL*8 :: ATRY(NMAX)
      REAL*8 :: BETA(NMAX)
      REAL*8 :: DA(NMAX)
      REAL*8 :: D(NMAX)
      REAL*8 :: E(NMAX)
      REAL*8 :: PSCALE(NMAX)
      REAL DCHISQ,OCHISQ,RATIO,FAC,EPS
C      DOUBLE PRECISION DELTA,DELTST,THRESH,DMAX,TEMP
      REAL*8 DELTA,DELTST,THRESH,DMAX,TEMP
      INTEGER J,K,KK,IHIT,ICOUNT,JCOUNT
      LOGICAL START,GOOD1
      SAVE ATRY,DELTA,BETA,OCHISQ,PSCALE,GOOD1
C
C---------------------------------------------------------------------
C
C========Do initializations
C
      START = .FALSE.
      ICOUNT = 0
      JCOUNT = 0
      FAC = 1.0
      EPS = 10.0
      IF(ALAMDA.LT.0.0)THEN
         KK = MFIT + 1
         DO 20 J=1,MA
            IHIT = 0
            DO 10 K=1,MFIT
               IF(LISTA(K).EQ.J)IHIT = IHIT + 1
 10         CONTINUE
            IF(IHIT.EQ.0)THEN
               LISTA(KK) = J
               KK = KK + 1
            ELSEIF(IHIT.GT.1)THEN
               STOP' Improper permutation in LISTA'
            ENDIF
 20      CONTINUE 
         IF(KK.NE.(MA+1))STOP' Improper permutation in LISTA'
C
C========Fill out ALPHA and BETA arrays
C
         CALL MRQCOF(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NCA,
     &               CHISQ,FUNCS,PSCALE)
         OCHISQ = CHISQ
         DO 30 J=1,MA
            ATRY(J) = A(J)
 30      CONTINUE
         START = .TRUE.
         GOOD1 = .FALSE.
C
C========Initialize step limit scaling to the current gradient
C
         DELTA = 0.0D0
         DO 32 J=1,MFIT
            DELTA = DELTA + (A(J)*PSCALE(J))**2
 32      CONTINUE
         DELTA = MIN(FACTOR*SQRT(DELTA),DELMAX)
         IF(DELTA.LT.DELMIN)DELTA = DELMIN
         ALAMDA = ALAMIN
      ENDIF
C
C========Copy ALPHA into COVAR
C
      DO 50 J=1,MFIT
         DO 40 K=1,J-1
            TEMP = ALPHA(K,J)
            COVAR(K,J) = TEMP
            COVAR(J,K) = TEMP 
 40      CONTINUE
 50   CONTINUE
C
C========Apply L-M factor to diagonal elements
C
 51   DO 52 J=1,MFIT
         TEMP = DBLE(1.0+ALAMDA)
         COVAR(J,J) = ALPHA(J,J)*TEMP
 52   CONTINUE
C
C========Find eigenvalues and eigenvectors of COVAR
C
      CALL TRED2(COVAR,MFIT,NCA,D,E)
      CALL TQLI(D,E,MFIT,NCA,COVAR)
C
C========Filter eigenvalues
C
      DMAX = 0.0D0
      DO 53 J=1,MFIT
         IF(ABS(D(J)).GT.DMAX)DMAX = D(J)
 53   CONTINUE
      IF(ALAMDA.GT.0.0)THEN
         THRESH = RMIN*DMAX
      ELSE
         THRESH = TEENSY
      ENDIF
      DO 54 J=1,MFIT
         IF(ABS(D(J)).LT.THRESH)D(J) = 0.0D0
 54   CONTINUE
      IF(ALAMDA.EQ.0.0)THEN
C
C========User wishes to call a halt - copy eigenvectors into ALPHA
C
         DO 56 J=1,MFIT
            DO 55 K=1,MFIT
               ALPHA(K,J) = COVAR(K,J)/PSCALE(K)
 55         CONTINUE
 56      CONTINUE
C
C========Construct covariance matrix from D and ALPHA
C
         CALL SYMINV(COVAR,D,ALPHA,MFIT,NCA)
C
C========Calculate covariances
C
         CALL COVSRT(COVAR,NCA,MA,LISTA,MFIT)
C
C========Will need to initialize ALPHA if called again
C
         ALAMDA = -1.0
         RETURN
      ENDIF
C
C========Back substitute to solve for BETA - result in DA
C
      CALL TBKSB(D,COVAR,MFIT,NCA,BETA,DA)
C
C========Test step against step limit and calculate expected reduction 
C========in CHISQ
C
      DELTST = 0.0D0
      DCHISQ = 0.0D0
      DO 57 J=1,MFIT
         TEMP = DA(J)
         DELTST = DELTST + TEMP*TEMP
         DCHISQ = DCHISQ + TEMP*BETA(J)
 57   CONTINUE
      DELTST = SQRT(DELTST)
C
C========Reinitialize DELTA to first stepsize if necessary
C
      IF(START)THEN
         DELTA = MIN(DELTA,DELTST)
         DELTA = MAX(DELTA,DELMIN)
      ENDIF
C
C========Try to find a solution close to the boundary of trust region
C
      ICOUNT = ICOUNT + 1
      IF(DELTST.GT.DELTA)THEN
         IF(FAC.LT.1.0)THEN
            JCOUNT = JCOUNT + 1
            EPS = SQRT(EPS)
         ENDIF
         FAC = EPS
         ALAMDA = FAC*ALAMDA
         IF(ALAMDA.LT.ALAMAX)THEN
            GOTO 51
         ELSE
            ALAMDA = ALAMAX
         ENDIF
      ELSEIF(DELTST.LT.0.9*DELTA)THEN
         IF(FAC.GT.1.0)THEN
            JCOUNT = JCOUNT + 1
            EPS = SQRT(EPS)
         ENDIF
         FAC = 1.0/EPS
         ALAMDA = FAC*ALAMDA
         IF(ALAMDA.GT.ALAMIN)THEN
            IF(ICOUNT.LT.20.AND.JCOUNT.LT.10)GOTO 51
         ELSE
            ALAMDA = ALAMIN
         ENDIF
      ENDIF
C
C========Calculate actual reduction in CHISQ
C
      DO 60 J=1,MFIT
         ATRY(LISTA(J)) = A(LISTA(J)) + DA(J)/PSCALE(J)
 60   CONTINUE
      CALL MRQCOF(X,Y,SIG,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,CHISQ,
     &            FUNCS,PSCALE)
C
C========Find ratio and test 
C
      IF(DCHISQ.GT.0.0)THEN
         RATIO = (OCHISQ-CHISQ)/DCHISQ
      ELSE
         RATIO = 0.0
      ENDIF
      IF(RATIO.GT.0.0001)THEN
C
C========Good enough to accept this step
C
         GOOD1 = .TRUE.
         DO 80 J=1,MFIT
            DO 70 K=1,J-1
               TEMP = COVAR(K,J)
               ALPHA(K,J) = TEMP
               ALPHA(J,K) = TEMP
 70         CONTINUE
            ALPHA(J,J) = COVAR(J,J)
            BETA(J) = DA(J)
            A(LISTA(J)) = ATRY(LISTA(J))
 80      CONTINUE
         OCHISQ = CHISQ
         IF(RATIO.GT.0.75)THEN
C     
C========Good enough to reset ALAMDA and double step size
C     
            ALAMDA = MAX(ALAMDA/2.0,ALAMIN)
            DELTA = MIN(2.0D0*DELTA,2.0D0*DELTST,DELMAX)
         ELSEIF(RATIO.LT.0.25)THEN
C
C========Not so good - double ALAMDA and half step size
C
            ALAMDA = MIN(2.0*ALAMDA,ALAMAX)
            DELTA = MIN(DELTA,DELTST)/2.0D0
            DELTA = MAX(DELTA,DELMIN)
         ENDIF
      ELSE
C
C========This step is no good.
C========Increase ALAMDA by a factor of 10 and reduce DELTA likewise
C========Negate ALAMDA to signal possible convergence if we've had at 
C========least one good step
C
         IF(RATIO.GT.0.0)THEN
            DELTA = MIN(DELTA,DELTST)/2.0D0
            ALAMDA = MIN(2.0*ALAMDA,ALAMAX)
         ELSE
            DELTA = MIN(DELTA,DELTST)/10.0D0
            ALAMDA = MIN(10.0*ALAMDA,ALAMAX)
         ENDIF
         DELTA = MAX(DELTA,DELMIN)
         IF(GOOD1)ALAMDA = -ALAMDA
      ENDIF
      RETURN
      END

C     LAST UPDATE 03/12/92
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE MRQCOF(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NALP,
     &                  CHISQ,FUNCS,PSCALE)
      IMPLICIT NONE
C
C Purpose: Used by MRQMIN to evaluate the linearized fitting matrix
C          ALPHA and vector BETA. 
C
      INTEGER NDATA,NALP,MA,MFIT
C      REAL X(NDATA),Y(NDATA),SIG(NDATA)
      REAL*4 :: X(NDATA)
      REAL*4 :: Y(NDATA)
      REAL*4 :: SIG(NDATA)
C      REAL*8 ALPHA(NALP,NALP)
      REAL*8 :: ALPHA(NALP,NALP)
C      REAL*8 BETA(MA),A(NALP),PSCALE(MA)
      REAL*8 :: BETA(MA)
      REAL*8 :: A(NALP)
      REAL*8 :: PSCALE(MA)
C      INTEGER LISTA(NALP)
      INTEGER :: LISTA(NALP)
      REAL CHISQ
      EXTERNAL FUNCS
C
C Calls   1: FUNCS
C Called by: MRQMIN
C
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C Local variables:
C
      INTEGER NMAX
      PARAMETER(NMAX=50)
C      DOUBLE PRECISION TEENSY
      REAL*8 TEENSY
      PARAMETER(TEENSY=1.0D-10)
C      DOUBLE PRECISION DYDA(NMAX)   
      REAL*8 :: DYDA(NMAX)   
C      DOUBLE PRECISION DY,SIG2I,WT,YMOD,XI
      REAL*8 DY,SIG2I,WT,YMOD,XI
      INTEGER I,J,K 
C
C---------------------------------------------------------------------
      DO 20 J=1,MFIT
         DO 10 K=1,J
            ALPHA(J,K) = 0.0D0
 10      CONTINUE
         BETA(J) = 0.0D0
         PSCALE(J) = 0.0D0
 20   CONTINUE
      CHISQ = 0.0
      DO 50 I=1,NDATA
         IF(Y(I).GT.-0.99E+30)THEN
            XI = DBLE(X(I))
            CALL FUNCS(XI,A,YMOD,DYDA,MA)
            SIG2I = 1.0D0/(DBLE(SIG(I))*DBLE(SIG(I)))
            DY = DBLE(Y(I)) - YMOD
            DO 40 J=1,MFIT
               WT = DYDA(LISTA(J))*SIG2I
               PSCALE(J) = PSCALE(J) + WT*DYDA(LISTA(J))
               DO 30 K=1,J
                  ALPHA(J,K) = ALPHA(J,K) + WT*DYDA(LISTA(K))
 30            CONTINUE
               BETA(J) = BETA(J) + DY*WT
 40         CONTINUE
            CHISQ = CHISQ + SNGL(DY*DY*SIG2I)
         ENDIF 
 50   CONTINUE
      DO 70 J=1,MFIT
         IF(PSCALE(J).GT.TEENSY)THEN
            PSCALE(J) = SQRT(PSCALE(J))
         ELSE
            PSCALE(J) = 1.0D0
         ENDIF
         BETA(J) = BETA(J)/PSCALE(J)
         DO 60 K=1,J
            ALPHA(J,K) = ALPHA(J,K)/(PSCALE(J)*PSCALE(K))
            ALPHA(K,J) = ALPHA(J,K)
 60      CONTINUE
 70   CONTINUE 
      RETURN
      END
