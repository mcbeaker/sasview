C      TMWN JULY 94
C      Savitsky Golay smoothing subroutines.
C      Taken from "Numerical Recipes" by Flannery and Press.
C      Stored here in a separate file to assist debugging.

       SUBROUTINE savgol(c,np,nl,nr,ld,m)
C      Savitsky Golay smoothing.
C      Taken from Numerical Recipes, Flannery and Press.
C      Calculates coefficients for linear smoothing,
C      such that higher moments are preserved.
C      Re. correlation function analysis, smoothing used
C      to ease the join between calc. and expt. data.
       INTEGER ld,m,nl,np,nr,mmax
       REAL c(np)
       PARAMETER (mmax=6)
C       INTEGER imj,ipj,j,k,kk,mm,indx(mmax+1)
       INTEGER imj,ipj,j,k,kk,mm
       INTEGER :: indx(mmax+1)
       REAL d,fac,sums
C       REAL a(mmax+1,mmax+1),b(mmax+1)
       REAL*4 :: a(mmax+1,mmax+1)
       REAL*4 :: b(mmax+1)

      
       IF (np.LT.nl+nr+1 .OR.nl.LT.0 .OR.nr.LT.0 .OR.ld.GT.m
     & .OR.m.GT.mmax.OR.nl+nr.LT.m) then 
	      WRITE(6,*) 'Bad arguments passed to subroutine savgol'
       ENDIF
	  
       DO ipj=0,2*m
        sums=0.
        IF (ipj .EQ. 0) sums=1.

        DO k=1,nr
          sums=sums+FLOAT(k)**ipj
        END DO

        DO k=1,nl
          sums=sums+FLOAT(-k)**ipj
        END DO

        mm=MIN(ipj,2*m-ipj)
        DO imj=-mm,mm,2
          a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sums
        END DO
       END DO

       CALL ludcmp(a,m+1,mmax+1,indx,d)

       DO j=1,m+1
        b(j)=0.
       END DO
       b(ld+1)=1

       CALL lubksb(a,m+1,mmax+1,indx,b)

       DO kk=1,np
        c(kk)=0.
       END DO

       DO k=-nl,nr
        sums=b(1)
        fac=1.
        DO mm=1,m
          fac=fac*k
          sums=sums+b(mm+1)*fac
        END DO
        kk=MOD(np-k,np)+1
        c(kk)=sums
       END DO

      RETURN
      END



      SUBROUTINE ludcmp(a,n,np,indx,d)
C     LU decomposition.
C     Taken from Numerical Recipes, Flannery and Press.
C     Used by subroutine savgol.
      INTEGER n,np,nmax
      PARAMETER (nmax=500,teensy=1.E-8)
      INTEGER i,imax,j,k 
      REAL d,aamax,dum,sums,teensy
      INTEGER :: indx(n)
C      REAL a(np,np),vv(nmax)
      REAL*4 :: a(np,np)
      REAL*4 :: vv(nmax)

      d=1.
      DO i=1,n
        aamax=0.
        DO j=1,n
          IF (ABS(a(i,j)) .GT. aamax) aamax=ABS(a(i,j))
        END DO
        IF (aamax .EQ. 0.) write(6,*) 'Singular matrix in subroutine ludcmp'
        vv(i)=1./aamax
      END DO
      
      DO j=1,n
        DO i=1,j-1
          sums=a(i,j)
          DO k=1,i-1
            sums=sums-a(i,k)*a(k,j)
          END DO
          a(i,j)=sums
        END DO
        aamax=0.
        DO i=j,n
          sums=a(i,j)
          DO k=1,j-1
            sums=sums-a(i,k)*a(k,j)
          END DO
          a(i,j)=sums
          dum=vv(i)*ABS(sums)
          IF (dum .GE. aamax) THEN
            imax=i
            aamax=dum
          ENDIF
        END DO
        IF (j .NE. imax) THEN
          DO k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          END DO
          d=-d
          vv(imax)=vv(j)
        ENDIF
        indx(j)=imax
        IF (a(j,j) .EQ. 0.) a(j,j)=teensy
        IF (j .NE. n) THEN
          dum=1./a(j,j)
          DO i=j+1,n
            a(i,j)=a(i,j)*dum
          END DO
        ENDIF
      END DO
      
      RETURN
      END



      SUBROUTINE lubksb(a,n,np,indx,b)
C     Subroutine associated with LU decomposition.
C     Taken from "Numerical Recipes" Flannery and Press.
C     Used in Savitsky Golay smoothing.
      INTEGER n,np
      INTEGER i,ii,j,ll
      REAL sums
      INTEGER :: indx(n)
C      REAL a(np,np),b(n)
      REAL*4 :: a(np,np)
      REAL*4 :: b(n)

      ii=0
      DO i=1,n
        ll=indx(i)
        sums=b(ll)
        b(ll)=b(i)
        IF (ii .NE. 0) THEN
          DO j=ii,i-1
            sums=sums-a(i,j)*b(j)
          END DO
        ELSEIF (sums .NE. 0.) THEN
          ii=i
        ENDIF
        b(i)=sums
      END DO
      DO i=n,1,-1
        sums=b(i)
        DO j=i+1,n
          sums=sums-a(i,j)*b(j)
        END DO
        b(i)=sums/a(i,i)
      END DO

      RETURN
      END
