      SUBROUTINE sigmoid(x,y,dyda,param)
C     Levenburg Marquart subroutine.
C     Returns y and grad y (wrt fitting parameters)
C     at a point x, according to the sigmoid function
C     Iobs=B + (K/Q^4)*EXP(-sigma^2*Q^2)
C
C     Often leads to "IEEE" errors.
C     Debug in progress.
C
C      DIMENSION param(3)
C      DIMENSION dyda(3)
      REAL*4 :: param(3)
      REAL*4 :: dyda(3)
      REAL x,y

      dyda(1)=1.

      temp1=x*x
      temp4=1./x
      temp4=temp4*temp4
      temp2=-temp1*param(3)*param(3)
      IF (temp2 .GT. 16.) THEN
        temp3=1.E+8
      ELSEIF (temp2 .LT. -16.) THEN
        temp3=0.
      ELSE
        temp3=EXP(temp2)
      ENDIF

C     Estimate size of results
      xlg=LOG10(ABS(x))
      a=-4.*xlg
      IF (param(2) .NE. 0.) THEN
        porlg=LOG10(ABS(param(2)))
        c=porlg-4.*xlg
      ELSE
        porlg=0.
        c=0.
      ENDIF
      IF (param(3) .NE. 0. .AND. param(2) .NE. 0.) THEN
        siglg=LOG10(ABS(param(3)))
        b=siglg+porlg-2.*xlg
      ELSE
        b=0.
      ENDIF
        

C     Take care of y
      IF (c .GT. 8.) THEN
        y=1.E+8
      ELSEIF (c .LT. -8.) THEN
        y=0.
      ELSE
        y=param(2)*temp3
        y=y*temp4
        y=y*temp4
        y=y+param(3)
      ENDIF

C     Take care of deriv. wrt K
      IF (a .GT. 8.) THEN
        dyda(2)=1.E+8
      ELSEIF (a .LT. -8.) THEN
        dyda(2)=0.
      ELSE
        dyda(2)=temp3*temp4
        dyda(2)=dyda(2)*temp4
      ENDIF
      
C     Take care of deriv. wrt sig
      IF (b .GT. 8.) THEN
        dyda(3)=1.E+8
      ELSEIF (b .LT. -8.) THEN
        dyda(3)=0.
      ELSE
        dyda(3)=-2.*param(2)*param(3)*temp3
        dyda(3)=dyda(3)*temp4
      ENDIF

C      IF (ABS(y) .GT. 1.0E+8) THEN
C        write(6,*)'Debug warning!'
C      ENDIF

      RETURN
      END
