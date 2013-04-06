C     Polynomial Interpolation subroutine
C     Taken from "Numerical Recipes", 2nd Edition by Press, Teukolsky, Vetterling & Flannery.

      SUBROUTINE polint(xa,ya,n,x,y,dy)

      INTEGER n,NMAX
C      REAL dy,x,y,xa(n),ya(n)
      REAL dy,x,y
C	NMAX is the maximum anticipated polynomial order
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
C      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      REAL den,dif,dift,ho,hp,w
      REAL*4 :: xa(n)
      REAL*4 :: ya(n)
      REAL*4 :: c(NMAX)
      REAL*4 :: d(NMAX)

      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)write(6,*)'WARNING: polint failed'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
