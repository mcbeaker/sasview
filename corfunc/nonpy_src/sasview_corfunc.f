C     sasview_corfunc.f90
C     S King, April 2013

C     Data arrays are currently dimensioned for 4096 points (or multiples of)
C     Anything that was 512, 1024 or 2048 points is now 4096
C     DIMENSION statements have been changed to INTEGER ::, REAL*4 ::
C     or REAL*8 :: to satisfy f2py
C     Several variables had the names of intrinsic functions in later versions
C     of FORTRAN. These have been changed as a precaution.

C==============================================================================
C==============================================================================

      SUBROUTINE tailfit 
     &(filename, qaxname, sigmodel, limit1, limit2, realtime, ndata, 
     &numiteration, results, best, bestchannel, param, static_bn)
C
C     TMWN July 94
C     Performs tailfit
C
C     Background is Bonart, and tail sigmoidal.
C     See Koberstein and Stein, J.polymer.sci,phys.ed.,vol21,pp.2181
C
C     Non-linear fitting via Levenburg Marquart algorithm.
C     See Numerical Recipes by Flannery and Press, pp.679
C
C     Updated to handle realtime data.
C
C     30/7/94 Update.
C     Now choice of either sigmoid model or simple Porod.
C
C     5/8/94 Update
C     Singular matrix handling.
C
C     9/8/94 Update
C     Variable channel limits for tailfit.
C
C     Nov 05 Update
C     Redimensioned data arrays from 512 to MaxDim.
C     Added ASCII output of multi-frame data & array asciiresults

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

      EXTERNAL corfunc_io, wrdlen, getfiletype, readascii, getpath
      EXTERNAL addstr, byte_order, mrqmin, strpath, tred2, tqli
      EXTERNAL syminv, covsrt, tbksb, writeascii, loqread
C      EXTERNAL sigmoid

      INTEGER lim1,lim2,qzero,lower,upper,ij,i,npoints,ncof
      INTEGER plotflag,lorentzflag,lowerlim,upperlim
C      INTEGER best(MaxDim),datastart,lword
      INTEGER :: best(MaxDim)
      INTEGER datastart,lword,ndata,numiteration
      LOGICAL static_bn,good
      REAL gradient,intercept,lamda,oldchisqu,kp
      REAL chitst,estd,chisqu
      CHARACTER*80 filename,qaxname,othername,dirname
      CHARACTER*40 title,retrans,ascii,sigmodel,user,idftoggle
      CHARACTER*80 fname,axisname,fname2,arbabs,graphics,outname,stats,
C      CHARACTER*80 prompts(23),grlong(4),backex*40,stats,grname(4)
      CHARACTER*80 :: prompts(23)
      CHARACTER*80 :: grlong(4)
      CHARACTER*80 :: grname(4)
      CHARACTER*1 letter
      CHARACTER*40 graphdir,graphname,graphqax,graphtitle,backex
      CHARACTER*40 xlabel,ylabel
      CHARACTER*10 filetype
C      INTEGER bestchannel,realtime(4),realflag
C      INTEGER lista(3),OPTION
      INTEGER bestchannel,realflag,OPTION
      INTEGER :: lista(3)
      INTEGER :: realtime(4)
C      REAL xtemp(MaxDim),ytemp(MaxDim),stdev(MaxDim),temp
      REAL temp
      REAL*4 :: xtemp(MaxDim)
      REAL*4 :: ytemp(MaxDim)
      REAL*4 :: stdev(MaxDim)
      REAL*4 :: covar(3,3)
      REAL*4 :: alpha(3,3)
      REAL*4 :: param(3)
      REAL*4 :: newparam(3)
      REAL*4 :: asciidata(MaxDim*MaxDim)
C      DOUBLE PRECISION dcovar(3,3),dalpha(3,3),dparam(3)
      REAL*8 :: dcovar(3,3)
      REAL*8 :: dalpha(3,3)
      REAL*8 :: dparam(3)
      INTEGER :: notok(10)
C      DIMENSION xdata(MaxDim),ydata(MaxDim)
      REAL*4 :: xdata(MaxDim)
      REAL*4 :: ydata(MaxDim)
      REAL*4 :: deltaparam(3)
      REAL*4 :: beta(3)
C      DIMENSION data(MaxDim,MaxDim),results(MaxDim,5)
      REAL*4 :: dataz(MaxDim,MaxDim)
      REAL*4 :: results(MaxDim,5)
      REAL*4 :: asciiresults(MaxDim)
C      DIMENSION limit1(MaxDim),limit2(MaxDim)
      INTEGER :: limit1(MaxDim)
      INTEGER :: limit2(MaxDim)

C     Array "param"
C     Parameter 1: Bonart background
C     Parameter 2: K (interface per unit volume)
C     Parameter 3: sigma (diffuse boundary thickness

C     Array "data"
C     Stores intensities. dataz(frame no.,channel).
C     Added during realtime update.

C     "sigmodel" controls which model is applied.
C     If sigmodel='off' a straight forward Porod plot is done.

1000  FORMAT(A1)
1010  FORMAT(A80)
1020  FORMAT(2x,I3)
1030  FORMAT(10I8)
1040  FORMAT(A40)
1050  FORMAT(/,1x,'100: Initial guess at background: IQ^4 ',
     &'vs. Q^4 ',E12.6,'.',/,1x,
     &'100: Simple measurement ',E12.6,'.')
1060  FORMAT(/,1x,'100: Initial guess at sigmoid parameters:',
     &/,1x,'100: K ',E12.6,'.',/,1x,'100: Sigma ',E12.6,'.')
1070  FORMAT(10I8)
1080  FORMAT(A10)
1090  FORMAT(/,1x,'100: Parameter            Value',/,
     &     '----------------------------',/)
1100  FORMAT(1x,'100: ',A10,11x,E12.6,3x,E12.6)
1110  FORMAT(/,1x,'100: Chi squared          ',E12.6)
1120  FORMAT(/,1x,'Optimum channel for tail start = ',I3)
1130  FORMAT(E12.6)
1140  FORMAT(2x,I3,2x,I3,2x,I3,2x,I3)
1150  FORMAT(/,1x,'100: Frame  Stage    Background    ',
     &'K Porod       Sigma         Chisq ',
     &'        Std. Error    Reduced Chisq',/,
     &'----------------------------------------',
     &'----------------------------------------')
1160  FORMAT(1x,'100: ',I3,4x,A7,2x,E12.6,2x,E12.6,2x,E12.6,
     &2x,E12.6,2x,E12.6,2x,E12.6)
1170  FORMAT(1x,'100: ',I3,4x,A7,2x,E12.6,2x,E12.6,2x,E12.6)
1180  FORMAT(E12.6,2x,E12.6,2x,E12.6,2x,I3)
1190  FORMAT(/,1x,'WARNING: Fitting problem: Singular alpha matrix',/)
1200  FORMAT(2x,I3,2x,I3)
1210  FORMAT(1x,"100: Written thermal background v frame: ",A10)
1215  FORMAT(1x,"100: Written ASCII thermal background v frame: ",A10)
1220  FORMAT(1x,"100: Written Porod constant v frame: ",A10)
1225  FORMAT(1x,"100: Written ASCII Porod constant v frame: ",A10)
1230  FORMAT(1x,"100: Written sigma v frame: ",A10)
1235  FORMAT(1x,"100: Written ASCII sigma v frame: ",A10)
1240  FORMAT(1x,"100: Written x-axis: ",A10)
1241  FORMAT(1x,"100: Written ASCII x-axis: ",A10)
1245  FORMAT(1x,"100: Covergence achieved after ",I3," iterations")
1247  FORMAT(1x,"100: Standard error       ",E12.6)

      prompts(1)='ERROR: Error reading transfer parameter file: FATAL'
      prompts(2)='ERROR: Error reading data file: FATAL'
      prompts(3)='ERROR: Expecting static Q axis,'
     &//' received dynamic: FATAL'
      prompts(4)='All necessary files correctly loaded...'
      prompts(5)='ERROR: Error writing output otoko files: FATAL'
      prompts(6)='WARNING: Positive gradient during '
     &//'estimate of sigma...'
      prompts(7)='ERROR: Software error with subroutine'
     &//' changeotok: FATAL'
      prompts(8)='100: Using IQ^4 vs. Q^4 least squares background...'
      prompts(9)='WARNING: Ran out of iterations before convergence...'
      prompts(10)='NON-LINEAR FITTING OVER!'
      prompts(11)='100: Final parameter values:'
      prompts(12)='Updating parameter files...'
      prompts(13)='ERROR: Error writing file tailinfo.txt: FATAL'
      prompts(14)='Creating files for graphics output...'
      prompts(15)='100: Applying sigmoid tail model...'
      prompts(16)='100: Applying Porod tail model...'
      prompts(17)='ERROR: Error with status file: FATAL'
      prompts(18)='ERROR: Error reading limitinfo.txt: FATAL'
      prompts(19)='Started tailfitting'
      prompts(20)='Finished tailfitting'
      prompts(21)='100: TAIL-FITTING...'
      prompts(22)='100:'
C     Added by SMK, Nov 05
      prompts(23)='WARNING: Error writing ascii output'

      CALL WRDLEN(LWORD)
      limopt=0
      realflag=1
      if(static_bn)then
        realflag=0
        realtime(1)=1
        realtime(2)=1
        realtime(3)=1
        realtime(4)=1
      endif

      WRITE(6,*)
      WRITE(6,*)
C      title='TAIL-FITTING'
C      CALL showtitle(title)
      CALL showprompts(prompts(22))
      CALL showprompts(prompts(21))

      fname=filename
      axisname=qaxname
      CALL getpath(fname,dirname)

c     Read the data
      CALL getfiletype(filetype,fname)
      IF(filetype.EQ."ascii")THEN
          CALL readascii(fname,notok,asciidata,irc)
          IF(irc.EQ.0)THEN
              GOTO 5010
          ENDIF
          nndata=ndata*4
          ij=1
          DO nframe=realtime(1),realtime(2),realtime(3)
            DO i=1,ndata
              dataz(nframe,i)=asciidata(ij)
              ij=ij+1
            END DO
          END DO
      ELSE
          OPEN(UNIT=9,FILE=fname,STATUS='old',ERR=5010)
          READ(9,1000,ERR=5010)letter
          READ(9,1000,ERR=5010)letter
          READ(9,1030,ERR=5010)notok
          READ(9,1040,ERR=5010)othername
          CLOSE(9)

          nndata=4*ndata

C         READ intensities
          CALL addstrings(dirname,othername,fname2)
          OPEN(UNIT=9,FILE=fname2,STATUS='old',
     &    ACCESS='direct',RECL=nndata/lword,ERR=5010)
          DO nframe=realtime(1),realtime(2),realtime(3)
            READ(9,REC=nframe,ERR=5010)(dataz(nframe,i),i=1,ndata)
          END DO
          CLOSE(9)
      ENDIF

      IF(filetype.NE."ascii")THEN
C SMK   CALL swap(data,notok(4),512*512,4)
        CALL swap(data,notok(4),MaxDim*MaxDim,4)
      ENDIF

      nframe=1

C     Open Q axis
      CALL getfiletype(filetype,axisname)
      IF(filetype.EQ."ascii")THEN
          CALL readascii(axisname,notok,asciidata,irc)
          IF(irc.EQ.0)THEN
              GOTO 5010
          ENDIF
C         Check static
          IF (notok(2) .NE. 1) THEN
            CALL showprompts(prompts(3))
            STOP
          ENDIF
          DO i=1,ndata
            xdata(i)=asciidata(i)
          END DO
      ELSE
          OPEN(UNIT=9,FILE=axisname,STATUS='old',ERR=5010)
          READ(9,1000,ERR=5010)letter
          READ(9,1000,ERR=5010)letter
          READ(9,1030,ERR=5010)notok
          READ(9,1040,ERR=5010)othername
          CLOSE(9)

C         Check static
          IF (notok(2) .NE. 1) THEN
            CALL showprompts(prompts(3))
            STOP
          ENDIF

C         Read x axis data
          CALL addstrings(dirname,othername,fname2)
          OPEN(UNIT=9,FILE=fname2,STATUS='old',
     &    ACCESS='direct',RECL=nndata/lword,ERR=5010)
          READ(9,REC=nframe,ERR=5010)(xdata(i),i=1,ndata)
2020      CLOSE(9)
      ENDIF

      IF(filetype.NE."ascii")THEN
C SMK   CALL swap(xdata,notok(4),512,4)
        CALL swap(xdata,notok(4),MaxDim,4)
      ENDIF

C     Ok - data fine
      CALL showprompts(prompts(4))


C     START NUMBER CRUNCHING
C     ~~~~~~~~~~~~~~~~~~~~~

      IF (sigmodel .EQ. 'off' .OR. sigmodel .EQ. 'OFF') THEN
        CALL showprompts(prompts(16))
      ELSE
        CALL showprompts(prompts(15))
      ENDIF

C     Titles for text output
      IF (realflag .EQ. 1 .AND. realtime(4) .NE. 1) THEN
        WRITE(6,1150)
      ENDIF

C     Realtime update: loop through original static routine.
      DO nframe=realtime(1),realtime(2),realtime(3)

C     Load up data
      DO i=1,MaxDim
        ydata(i)=dataz(nframe,i)
      END DO

C     Set limits
      lim1=limit1(nframe)
      lim2=limit2(nframe)


C     Which model?
      IF (sigmodel .EQ. 'off' .OR. sigmodel .EQ. 'OFF') THEN

C       POROD MODEL!
C       Do an LSQ fit in I vs. 1/Q^4 world, with channel optimisation.
        DO i=1,3
          param(i)=0.
        END DO
        oldchisqu=1.E+20

C       optimise channel
        bestchannel=lim1-limopt
        DO newlim1=lim1-limopt,lim1+limopt

C         Load up data.
          npts=lim2-newlim1+1
          DO i=1,npts
            xtemp(i)=(1./(xdata(newlim1+i-1))**4)
            ytemp(i)=ydata(newlim1+i-1)
          END DO
C         LSQ fit
          CALL lsqfit(xtemp,ytemp,npts,intercept,gradient)

C         Get chisqu
          chisqu=0.
          DO i=newlim1,lim2
            xp4=xdata(i)**4
            chisqu=chisqu+(ydata(i)-intercept-(gradient/xp4))**2
          END DO
          chisqu=chisqu/npts

C         Improved?
          IF (chisqu .LT. oldchisqu) THEN
            oldchisqu=chisqu
            param(1)=intercept
            param(2)=gradient
            param(3)=0.
            bestchannel=newlim1
          ENDIF

C         Debug:
C          write(6,*)bestchannel,chisqu

C       End channel optimisation
        END DO

C       Output results.
        IF (realflag .EQ. 0 .OR. realtime(4) .EQ. 1) THEN
C         single frame
          WRITE(6,*)
          CALL showprompts(prompts(11))

          WRITE(6,1090)
          WRITE(6,1100)'Background',param(1)
          WRITE(6,1100)'K',param(2)
          WRITE(6,1110)chisqu
          WRITE(6,1120)bestchannel
          best(nframe)=bestchannel

        ELSE
C         More than one frame.
          best(nframe)=bestchannel
          DO i=1,3
            results(nframe,i)=param(i)
          END DO
          results(nframe,4)=chisqu
          WRITE(6,1160)nframe,'Final',param(1),param(2),
     &    param(3),chisqu

        ENDIF

      ELSE

C       SIGMOID MODEL!

C       Get initial guess for parameters
C       First: plot IQ^4 vs Q^4 and measure slope as guess for backgr.
C       This might not work very well: if not, change to a direct
C       measurement of intensity at high Q.
        npts=lim2-lim1+1
        DO i=1,npts
          xtemp(i)=(xdata(lim1+i-1))**4
          ytemp(i)=ydata(lim1+i-1)*xtemp(i)
        END DO
        CALL lsqfit(xtemp,ytemp,npts,intercept,gradient)
        guess1=gradient

C       Or use a simple measurement:
        guess2=0.
        DO i=lim2-15,lim2
          guess2=guess2+ydata(i)
        END DO
        guess2=guess2/16.0

C       More may be added here re. two methods of finding backgr.
        param(1)=guess1

C       Next guess at sigma and Kp.
C       Plot ln{[Iobs-Ib]q^4} vs q^2
C       Slope=-sigma**2, intercept=ln Kp
        nptstemp=0
        DO i=1,npts
          temp1=(ydata(i+lim1-1)-param(1))
          IF (temp1 .GT. 0.) THEN
            nptstemp=nptstemp+1
            temp2=xdata(i+lim1-1)**2
            xtemp(nptstemp)=temp2
            ytemp(nptstemp)=LOG(temp1*temp2*temp2)
C            debug:
C            write(6,*)(i+lim1-1),xtemp(nptstemp),ytemp(nptstemp)
          ENDIF
        END DO

        CALL lsqfit(xtemp,ytemp,nptstemp,intercept,gradient)

        param(2)=EXP(intercept)
        IF (gradient .GT. 0.) THEN
C         param(3)=0.   ???
C         Not really sure what to do: but a start of 0.0 always
C         leads to end of 0.0 so use:-
          param(3)=SQRT(ABS(gradient))
        ELSE
          param(3)=SQRT(-gradient)
        ENDIF

C       Debug
C        param(1)=0.2
C        param(2)=0.5E-3
C        param(3)=2.0


C       DEBUG
C       ~~~~~
C       Loop through values of sigma from 0 to 10.
C       For each value do an LSQ to determine backgr. and Porod.
C       Then calculate a chi squared for comparison with
C       final result.
C       Ignores channel optimisation for simplicity.
        ndebug=0
        IF (ndebug .EQ. 1) THEN
          WRITE(6,*)
          num=lim2-lim1+1
          DO k=0,10
            sig=FLOAT(k)
            sigsqu=sig*sig
            sum1=0.
            sum2=0.
            sum3=0.
            sum4=0.
            DO i=lim1,lim2
              x=xdata(i)
              y=ydata(i)
              temp=(EXP(-x*x*sigsqu))/(x**4)
              sum1=sum1+y
              sum2=sum2+temp
              sum3=sum3+y*temp
              sum4=sum4+temp*temp
            END DO
C           Evaluate B and Kp
            kp=(sum1*sum2-num*sum3)/(sum2*sum2-num*sum4)
            b=(sum1-kp*sum2)/num
C           Get chisqu
            chisqu=0.
            DO i=lim1,lim2
              x=xdata(i)
              y=ydata(i)
              chisqu=chisqu+(y-b-(kp/(x**4))*(EXP(-sigsqu*x*x)))**2
            END DO
            chisqu=chisqu/num
C           Output
7000        FORMAT(1x,'100: Sigma=',E12.6,'. K=',E12.6,
     &      '. B=',E12.6,'. Chisqu=',E12.6,'.')
            WRITE(6,7000)sig,kp,b,chisqu
          END DO
        ENDIF

C       END DEBUG
C       ~~~~~~~~~



C       Text output for single frame
        IF (realflag .EQ. 0 .OR. realtime(4) .EQ. 1) THEN
          WRITE(6,1050)guess1,guess2
          CALL showprompts(prompts(8))
          IF (gradient .GT. 0.) THEN
            CALL showprompts(prompts(6))
          ENDIF
          WRITE(6,1060)param(2),param(3)

        ELSE
C         text output for more than one frame
          WRITE(6,*)
          WRITE(6,1170)nframe,'Initial',param(1),param(2),param(3)

        ENDIF

      OPTION=2

      IF(OPTION.EQ.1)THEN

C       START LEVENBURG MARQUART FIT
C       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C       This is the equivalent of the Numerical Recipes subroutine
C       mrqmin, together with an external program.

C       Optimise position of tail start: search through all
C       possible starting points according to the user supplied
C       variable limopt.
        bestchannel=lim1
        DO newlim1=lim1-limopt,lim1+limopt

C         Initialise
C         Debug: may need changing
          lamda=0.001
C         The larger lamda is, the smaller the fitting steps.
C         This initial value may need changing later on.
          iteration=0

C         Do contents of mrqmin when lamda=-ve
C         Fill matrix alpha and vector beta; get chisqu for
C         initial parameter set.
          CALL mrqcof2(xdata,ydata,newlim1,lim2,param,
     &    alpha,beta,oldchisqu)


C         Start one iteration
2030      iteration=iteration+1

C         Show location.
C         Debug: replace display of information eventually.
C          CALL showpos(newlim1,iteration,param,alpha,
C     &    beta,lamda,oldchisqu)

C         Augment diagonal of linearized fitting matrix
          DO j=1,3
            DO k=1,3
              covar(j,k)=alpha(j,k)
            END DO
            covar(j,j)=alpha(j,j)*(1.+lamda)
C           Copy beta into array deltaparam for gaussj
            deltaparam(j)=beta(j)
          END DO

C         solve linear equns.
          CALL gaussj(covar,deltaparam,nerr)
          IF (nerr .EQ. 1) THEN
C           Singular matrix occurred!
C           Default to porod profile
C           Load up data.
            npts=lim2-newlim1+1
            DO i=1,npts
              xtemp(i)=(1./(xdata(newlim1+i-1))**4)
              ytemp(i)=ydata(newlim1+i-1)
            END DO
C           LSQ fit
            CALL lsqfit(xtemp,ytemp,npts,intercept,gradient)
            param(1)=intercept
            param(2)=gradient
            param(3)=0.
C           Ouput message and exit loop
            WRITE(6,1190)
            bestchannel=lim1
            GOTO 2050
          ENDIF

C         Increment params
          DO j=1,3
            newparam(j)=param(j)+deltaparam(j)
          END DO

C         Redo chisqu calculation.
C         Uncertainty: call mrqcov with covar or alpha?
C         I think covar is correct.
          CALL mrqcof2(xdata,ydata,newlim1,lim2,newparam,
     &    covar,deltaparam,chisqu)

C         Test the new solution
          IF (newparam(1) .LT. 0. .OR.
     &        newparam(2) .LT. 0. .OR.
     &        newparam(3) .LT. 0.) THEN
            chisqu=1.E+20
          ENDIF

          IF (chisqu .LT. oldchisqu) THEN
C           Reduce lamda for decreased sensitivity
C           Debug added:
C           cap lamda from getting too small.
            IF (lamda .GT. 0.001) THEN
              lamda=lamda*0.1
            ENDIF
            oldchisqu=chisqu
C           record which channel is best
            bestchannel=newlim1

C           update alpha, beta, and parameters
            DO j=1,3
              DO k=1,3
                alpha(j,k)=covar(j,k)
              END DO
              beta(j)=deltaparam(j)
              param(j)=newparam(j)
            END DO

          ELSE

C           chisqu has increased
            chisqu=oldchisqu
C           Increase lamda for increased sensitivity
            lamda=10.*lamda

          ENDIF

C         Now test for convergance or too many iterations
          IF (iteration .EQ. numiteration) THEN
            CALL showprompts(prompts(9))
            WRITE(6,*)
          ENDIF

C         lamda=1.E+8 is an arbitrary measure of convergence.
C         May need changing.
          IF ((iteration .LT. numiteration) .AND.
     &    (lamda .LT. 1.0E+8)) THEN
             GOTO 2030
          ENDIF

C         Then onto next channel for optimisation
        END DO


C       END LEVENBURG MARQUART FIT
C       ~~~~~~~~~~~~~~~~~~~~~~~~~~

 2050 chisqu=oldchisqu

        ELSE

C Set parameters to refined parameters from last fit
        dparam(1)=param(1)
        dparam(2)=param(2)
        dparam(3)=param(3)

C ***********************************************************************************
C Implement a different Levenburg-Marquart fitting routine of R. Denny: M.W. Shotton
C ***********************************************************************************

C On first call set lamda<0 (which then sets lamda=0.001: see subroutine mrqmin)
      lamda = -1.0
C Number of coefficients (param) for sigmoid is 3
      ncof=3
C Set target chi
      chitst=0.1
C Set up lista that flags values to be adjusted
      lista(1)=1
      lista(2)=2
      lista(3)=3

C Set up data arrays xtemp and ytemp containing points between lim1 and lim2

      i=1
      DO 2040 ij=lim1,lim2
          xtemp(i)=xdata(ij)
          ytemp(i)=ydata(ij)
          stdev(i)=1.0
          i=i+1
 2040 CONTINUE
      npoints=i-1

      good=.FALSE.

C Loop over maximum number of iterations and do fitting
      DO 2045 ij=1,numiteration

      CALL mrqmin(xtemp,ytemp,stdev,npoints,dparam,ncof,lista,ncof,
     &            dcovar,dalpha,ncof,chisqu,sigmoid,lamda)

C Convergence Test
      IF(ij.GT.1)THEN
        IF(lamda.LT.0.0)THEN
          good=.TRUE.
          temp = (oldchisqu-chisqu)*FLOAT(lim2-lim1)/chisqu
C          chisqu = oldchisqu
          lamda=-lamda
          IF(ABS(temp).lt.chitst)THEN
            WRITE(6,1245)ij
            GOTO 2046
          ENDIF
        ENDIF
      ENDIF

      oldchisqu = chisqu

 2045 CONTINUE

C Not reached convergence
      CALL showprompts(prompts(9))

 2046 CONTINUE

C Calculate standard error
      IF(npoints.GT.ncof)THEN
        estd = SQRT(chisqu/FLOAT(npoints-ncof))
      ELSE
        estd=0.0
      ENDIF

      bestchannel=lim1

      param(1)=dparam(1)
      param(2)=dparam(2)
      param(3)=dparam(3)
      if(param(3).lt.0)then
        param(3)=-param(3)
      endif

C ***********************************************************************************
C End of Levenburg-Marquart fitting
C ***********************************************************************************

      ENDIF

C       Output results to screen.
        IF (realflag .EQ. 0 .OR. realtime(4) .EQ. 1) THEN
C         single frame
          WRITE(6,*)
          CALL showprompts(prompts(10))
          CALL showprompts(prompts(11))

          WRITE(6,1090)
          WRITE(6,1100)'Background',param(1)
          WRITE(6,1100)'K',param(2)
          WRITE(6,1100)'Sigma',param(3)
          WRITE(6,1110)chisqu
          WRITE(6,1247)estd
          if(good)then
            WRITE(6,2047)temp
          endif
          WRITE(6,1120)bestchannel
 2047 format(" 100: Reduced Chi squared  ",E12.6)

        ELSE
C         More than one frame.
          best(nframe)=bestchannel
          DO i=1,3
            results(nframe,i)=param(i)
          END DO
          results(nframe,4)=chisqu
          WRITE(6,1160)nframe,'Final',param(1),param(2),
     &    param(3),chisqu,estd,temp

        ENDIF

C     End choice of models
      ENDIF

C     End loop through frames
      END DO


C     END NUMBER CRUNCHING
C     ~~~~~~~~~~~~~~~~~~~~

      IF (realflag .NE. 0 .OR. realtime(4) .NE. 1) THEN
C       Create graphs of backgr. vs. frame, K vs. frame,
C       and sigma vs. frame

C       Get filenames.
        DO i=1,4
          CALL strippath(filename,grname(i))
        END DO
        CALL swapexten(grname(1),'BAK')
        CALL swapexten(grname(2),'POR')
        CALL swapexten(grname(3),'SIG')
        CALL swapexten(grname(4),'FAX')


C       Create otoko header - xdata
C       Change header name to give data filename
        CALL changeotok(grname(4),othername,nerr)
        IF (nerr .EQ. 1) THEN
          CALL showprompts(prompts(7))
        STOP
        ENDIF
        grlong(4)=grname(4)
        OPEN(UNIT=9,FILE=grlong(4),STATUS='unknown',ERR=5020)
        WRITE(9,*,ERR=5020)
        WRITE(9,*,ERR=5020)
        notok(1)=realtime(4)
        notok(2)=1
        notok(3)=1
        CALL endian(notok(4))
        DO i=5,10
          notok(i)=0
        END DO
        WRITE(9,1070,ERR=5020)notok
        WRITE(9,1080,ERR=5020)othername(1:10)
        CLOSE(9)
C       Header done

C       Now do x data
C       Is RECL correct? For 512 .XAX it was 2048
        nrecl=MaxDim
        outname=othername
        OPEN(UNIT=10,FILE=outname,STATUS='unknown',
     &  ACCESS='direct',RECL=nrecl/lword,ERR=5020)
        WRITE(10,REC=1,ERR=5020)
     &  (FLOAT(i),i=realtime(1),realtime(2),realtime(3))
        CLOSE(10)

        write(6,1240)grlong(4)

C       create otoko header files - ydata
        IF (sigmodel .EQ. 'off' .OR. sigmodel .EQ. 'OFF') THEN
          nend=2
        ELSE
          nend=3
        ENDIF
        DO n=1,nend
C         n=1 ~ backgr
C         n=2 ~ Kporod
C         n=3 ~ sigma

C         Change header name to give data filename
          CALL changeotok(grname(n),othername,nerr)
          IF (nerr .EQ. 1) THEN
            CALL showprompts(prompts(7))
            STOP
          ENDIF
          grlong(n)=grname(n)
          OPEN(UNIT=9,FILE=grlong(n),STATUS='unknown',ERR=5020)
          WRITE(9,*,ERR=5020)
          WRITE(9,*,ERR=5020)
          WRITE(9,1070,ERR=5020)notok
          WRITE(9,1080,ERR=5020)othername(1:10)
          CLOSE(9)
C         Header done

C         Now do y data
          outname=othername
          OPEN(UNIT=10,FILE=outname,STATUS='unknown',
     &    ACCESS='direct',RECL=nrecl/lword,ERR=5020)
          WRITE(10,REC=1,ERR=5020)
     &    (results(i,n),i=realtime(1),realtime(2),realtime(3))
          CLOSE(10)

          if(n.eq.1)then
            write(6,1210)grlong(n)
          else if(n.eq.2) then
            write(6,1220)grlong(n)
          else if(n.eq.3)then
            write(6,1230)grlong(n)
          endif
        END DO

      ENDIF

      IF (realflag .NE. 0 .OR. realtime(4) .NE. 1) THEN

C       Output the thermal background data in ASCII format
       do i=realtime(1),realtime(2),realtime(3)       
           asciiresults(i)=results(i,1)
        end do
        CALL writeascii(grname(1)(1:3)//"BAK."//"TXT",notok,asciiresults,
     &                  irc,1)
        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(23))
        ELSE
          WRITE(6,1215)grname(1)(1:3)//"BAK."//"TXT"
        ENDIF

C       Output the Porod background data in ASCII format
       do i=realtime(1),realtime(2),realtime(3)       
           asciiresults(i)=results(i,2)
        end do
        CALL writeascii(grname(2)(1:3)//"POR."//"TXT",notok,asciiresults,
     &                  irc,1)
        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(23))
        ELSE
          WRITE(6,1225)grname(2)(1:3)//"POR."//"TXT"
        ENDIF

C       Output the second moment data in ASCII format
       do i=realtime(1),realtime(2),realtime(3)       
           asciiresults(i)=results(i,3)
        end do
        CALL writeascii(grname(3)(1:3)//"SIG."//"TXT",notok,asciiresults,
     &                  irc,1)
        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(23))
        ELSE
          WRITE(6,1235)grname(3)(1:3)//"SIG."//"TXT"
        ENDIF

C       Output the frame axis in ASCII format
       do i=realtime(1),realtime(2),realtime(3)       
           asciiresults(i)=i
        end do
        CALL writeascii(grname(4)(1:3)//"FAX."//"TXT",notok,asciiresults,
     &                  irc,1)
        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(23))
        ELSE
          WRITE(6,1241)grname(4)(1:3)//"FAX."//"TXT"
        ENDIF

      endif

      CALL showprompts(prompts(20))
      RETURN
      STOP


C     Error messages

C     Error reading in intensities or x data
5010  CALL showprompts(prompts(2))
      STOP

C     Error writing output otoko files
5020  CALL showprompts(prompts(5))
      STOP

C     Static image - can't get realtime info
5040  realflag=0
      realtime(1)=1
      realtime(2)=1
      realtime(3)=1
      realtime(4)=1

      RETURN
      END



      SUBROUTINE mrqcof2(xdata,ydata,lim1,lim2,param,
     &alpha,beta,chisqu)
C     Subroutine involved in Levenburg Marquart non linear fitting.
C     Adapted from Numerical Recipes, Flannery and Press.
C     Fills the matrix "alpha" with values. Alpha is symmetric.
C     Also fills vector beta.
C     Finally returns a chi squared value for the current fit.

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

C      DIMENSION xdata(MaxDim),ydata(MaxDim)
      REAL*4 :: xdata(MaxDim)
      REAL*4 :: ydata(MaxDim)
C      DIMENSION alpha(3,3),beta(3),dyda(3),param(3)
      REAL*4 :: alpha(3,3)
      REAL*4 :: beta(3)
      REAL*4 :: dyda(3)
      REAL*4 :: param(3)
      INTEGER lim1,lim2
      REAL chisqu

      DO j=1,3
        DO k=1,j
          alpha(j,k)=0.
        END DO
        beta(j)=0.
      END DO

C     calc chisqu
      chisqu=0.

C     loop over data points
      DO i=lim1,lim2
        CALL sigmoid2(xdata(i),y,dyda,param)
        dy=ydata(i)-y
        DO j=1,3
          wt=dyda(j)
          DO k=1,j
            alpha(j,k)=alpha(j,k)+wt*dyda(k)
          END DO
          beta(j)=beta(j)+dy*wt
        END DO
        chisqu=chisqu+dy*dy
      END DO
C     Normalise chisqu to no. of data points
      chisqu=chisqu/(1.*(lim2-lim1+1))

C     Fill remaining places in alpha by symmetry.
C     At the moment alpha is small so just do it 'by hand'.
      alpha(1,2)=alpha(2,1)
      alpha(1,3)=alpha(3,1)
      alpha(2,3)=alpha(3,2)

      RETURN
      END



      SUBROUTINE sigmoid(x,param,y,dyda,ma)
C     Levenburg Marquart subroutine.
C     Returns y and grad y (wrt fitting parameters)
C     at a point x, according to the sigmoid function
C     Iobs=B + (K/Q^4)*EXP(-sigma^2*Q^2)
C
C     Often leads to "IEEE" errors.
C     Debug in progress.
C

C      DOUBLE PRECISION x,y,param(3),dyda(3)
      REAL*8 x,y
      REAL*8 :: param(3)
      REAL*8 :: dyda(3)
      INTEGER ma

      dyda(1)=1.
      y=param(1)+(param(2)/(x**4))*EXP(-(param(3)**2)*(x**2))
      dyda(2)=EXP(-(param(3)**2)*(x**2))/(x**4)
      dyda(3)=-2.*param(3)*param(2)*EXP(-(param(3)**2)*(x**2))/(x**2)

      RETURN
      END



      SUBROUTINE sigmoid2(x,y,dyda,param)
C     Levenburg Marquart subroutine.
C     Returns y and grad y (wrt fitting parameters)
C     at a point x, according to the sigmoid function
C     Iobs=B + (K/Q^4)*EXP(-sigma^2*Q^2)
C
C     Often leads to "IEEE" errors.
C     Debug in progress.
C
C      DIMENSION dyda(3),param(3)
      REAL*4 :: dyda(3)
      REAL*4 :: param(3)
      REAL x,y

      dyda(1)=1.
      y=param(1)+(param(2)/(x**4))*EXP(-(param(3)**2)*(x**2))
      dyda(2)=EXP(-(param(3)**2)*(x**2))/(x**4)
      dyda(3)=-2.*param(3)*param(2)*EXP(-(param(3)**2)*(x**2))/(x**2)


      RETURN
      END


      SUBROUTINE gaussj(a,b,nerr)
C     Gauss-Jordan elimination algorithm
C     Solves linear equations.
C     From Numerical Recipes (Flannery and Press), pp.27
C     a is an n*n matrix, b an n-vector.
C     Gaussj inverts a and places soln. to equn. in b.
C     Max dimension of a is three.
C
C     There many be faster ways of solving linear equns. in 3
C     variables, but for a more complex background (eg. Vonk
C     instead of Bonart) the number of variables would increase.
C     Hence, for the moment, use gauss elimination.
C
C     5/8/94 Problems with singular matrices encountered:
C     Now when a sing. matrix is found, gaussj sets the output
C     vector b to zero ie. no chabge in parameters.
C     Hence chisqu won't change, and the program will behave as
C     if there was an increase in chisqu ie. new params discarded.
C
C     After experimentation it was decided that when sing.
C     matrices occur the user should be advised to re-run that frame
C     with a Porod approach.
C
C     Finally: built in auto default to sig=0. on failure.

C      DIMENSION a(3,3),b(3),ipiv(3),indxr(3),indxc(3)
      REAL*4 :: a(3,3)
      REAL*4 :: b(3)
      REAL*4 :: ipiv(3)
      REAL*4 :: indxr(3)
      REAL*4 :: indxc(3)
      INTEGER nerr

1000  FORMAT(/,1x,'100: Fitting error: singular alpha matrix!',/,1x,
     &'Re-run this frame using a Porod profile instead of',
     &' the sigmoid tail.',/,1x,'Error fatal (sorry)...',/)

C     Error handling:
      nerr=0

C     Three variables
      n=3

C     Initialise pivot array
      DO j=1,n
        ipiv(j)=0
      END DO

C     Find pivot elt.
      DO i=1,n
        big=0.
        DO j=1,n
          IF (ipiv(j) .NE. 1) THEN
            DO k=1,n
              IF (ipiv(k) .EQ. 0) THEN
                IF (ABS(a(j,k)) .GE. big) THEN
                  big=ABS(a(j,k))
                  irow=j
                  icol=k
                ENDIF

              ELSEIF (ipiv(k) .GT. 1) THEN
C               MTX is singular
C               PAUSE 'Singular matrix in subroutine gaussj'
C               WRITE(6,1000)
C               STOP
                nerr=1
                GOTO 2000
              ENDIF

            END DO
          ENDIF
        END DO
        ipiv(icol)=ipiv(icol)+1

C       We've got the pivot elt. so swap rows if nec. such that
C       pivot elt. on diagonal. Performed via column relabelling.
C       indxc(i) is the column of the ith pivotal elt. and is
C       ith to be reduced. indxr(i) is row of ith pivotal elt..
        IF (irow .NE. icol) THEN
          DO l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
          END DO
          dum=b(irow)
          b(irow)=b(icol)
          b(icol)=dum
        ENDIF

C       Now divide pivot row by pivot elt.
        indxr(i)=irow
        indxc(i)=icol
C       Check pivot elt. non-zero
        IF (a(icol,icol) .EQ. 0.) THEN
C         Singular MTX
C         PAUSE 'Singular matrix in subroutine gaussj'
C         WRITE(6,1000)
C         STOP
          nerr=1
          GOTO 2000
C         a(icol,icol)=1.E-8
        ENDIF
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        DO l=1,n
          a(icol,l)=a(icol,l)*pivinv
        END DO
        b(icol)=b(icol)*pivinv

C       Division complete.
C       Now reduce rows, subtracting row of pivot elt.
C       Remember not to change pivot
        DO ll=1,n
          IF (ll .NE. icol) THEN
C           Skip pivot
            dum=a(ll,icol)
            a(ll,icol)=0.
            DO l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
            END DO
            b(ll)=b(ll)-b(icol)*dum
          ENDIF
        END DO

      END DO
C     End of main loop.
C     Now unscramble solution from column interchanges.
C     ie. swap in reverse order.
      DO l=n,1,-1
        IF (indxr(l) .NE. indxc(l)) THEN
          DO k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
          END DO
        ENDIF
      END DO

2000  RETURN
      END



      SUBROUTINE showpos(newlim1,iteration,param,alpha,
     &beta,lamda,chisqu)
C     Outputs current fit location to screen.
C      DIMENSION param(3),alpha(3,3),beta(3)
      REAL*4 :: param(3)
      REAL*4 :: alpha(3,3)
      REAL*4 :: beta(3)
      REAL lamda,chisqu
      INTEGER newlim1,iteration

1000  FORMAT(///,1x,'CURRENT FIT LOCATION: iteration ',I3,
     &'. Tail starts at channel ',I3,'.')
1010  FORMAT(/,2x,'Parameter       Value        Gradient',/,
     &          '------------------------------------------',/)
1020  FORMAT(1x,A10,3x,E12.6,3x,E12.6)
1030  FORMAT(/,1x,'Alpha Matrix:')
1040  FORMAT(1x,E12.6,4x,E12.6,4x,E12.6)
1050  FORMAT(/,1x,'Lamda = ',E12.6,'.')
1060  FORMAT(/,1x,'Chi squared = ',E12.6,'.')

      WRITE(6,1000)iteration,newlim1
      WRITE(6,1060)chisqu
      WRITE(6,1010)
      WRITE(6,1020)'Background',param(1),beta(1)
      WRITE(6,1020)'K',param(2),beta(2)
      WRITE(6,1020)'Sigma',param(3),beta(3)
      WRITE(6,1030)
      DO i=1,3
        WRITE(6,1040)alpha(i,1),alpha(i,2),alpha(i,3)
      END DO
      WRITE(6,1050)lamda

      RETURN
      END



      SUBROUTINE lsqfit(x,y,n,a,b)
C     Linear least squares fit to y=a+bx
C     No chisquared or statistical errors given.

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)
      INTEGER n
      REAL a,b

C      DIMENSION x(MaxDim),y(MaxDim)
      REAL*4 :: x(MaxDim)
      REAL*4 :: y(MaxDim)

      sx=0.
      sy=0.
      st2=0.
      b=0.
      ss=float(n)

C     sigma x and sigma y
      DO i=1,n
        sx=sx+x(i)
        sy=sy+y(i)
      END DO
      sxoss=sx/ss

      DO i=1,n
        t=x(i)-sxoss
        st2=st2+t*t
        b=b+t*y(i)
      END DO

      b=b/st2
      a=(sy-sx*b)/ss

      RETURN
      END


C==============================================================================
C==============================================================================


      SUBROUTINE tailjoin 
     &(filename, qaxname, param, channel, ndata, qzero, backex, 
     &datastart, realtime, limit2, static_bn)

C     TMW NYE July 94
C     Correlation function analysis module.
C     Extrapolates experimental scattering data to "Q infinity",
C     actually an arbitrarily large Q value, according to the
C     parameters found in program tailfit
C
C     Then applies a Savitsky-Golay smoothing filter to the
C     transition between experimental and calculated data,
C     to ensure higher moments - and hence the region D ~ 2PI/(Qjoin) -
C     are not biased by the joining process.
C     REF: Numerical Recipes, Flannery and Press, p.646
C
C     Also performs back-extrapolation to Q=0 by least squares fitting
C     according to one of two models:
C     Vonk:    I(Q) = h1-h2*Q^2
C     Guinier: I(Q) = A * EXP(B*Q^2)
C     REF? No ref for vonk, but Guinier is standard:
C     X-ray scattering of synth. polymers; Balta-Calleja & Vonk.
C
C     6/8/94 Update:
C     Smoothing takes place in a Log(I) Q world, so problems arise
C     when tails are noisey and the background exceeds expt. data
C     points (attempted log of negative no.). Until now attempts
C     to remedy this have been poor so to simplify the situation
C     smoothing will take place on non backgr. subtracted data,
C     after which a backgr. will be subtracted.
C     Also testing effect at high R on gamma1(R).
C
C     9/8/94 Update:
C     Varible tailfit limits cf. frame no.
C
C     7/7/04 Update: S King
C     Added ASCII output of extrapolated data.
C
C     28/11/05 Update: S King
C     Redimensioned data arrays from 512 to MaxDim.
C     Added ASCII output of multi-frame data
C

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

      EXTERNAL corfunc_io, wrdlen, getfiletype, readascii, getpath
      EXTERNAL addstr, byte_order, mrqmin, strpath, tred2, tqli
      EXTERNAL syminv, covsrt, tbksb, writeascii, loqread

      CHARACTER*80 dirname,filename,qaxname,othername,xheader,yheader,
     &             yheader2
      CHARACTER*40 title,retrans,ascii,sigmodel
C      CHARACTER*80 prompts(28),fname,axisname,arbabs,graphics
      CHARACTER*80 :: prompts(28)
      CHARACTER*80 fname,axisname,arbabs,graphics
      CHARACTER*80 fname2,outname,stats,radname,radax
      CHARACTER*40 backex,user,idftoggle
      CHARACTER*4 statuss
      CHARACTER*1 letter
      CHARACTER*40 xlabel,ylabel
      CHARACTER*10 filetype
      LOGICAL static_bn
      INTEGER lim1,lim2,qzero,datastart,ndata
C      INTEGER realtime(4),realflag,start,chanend
      INTEGER :: realtime(4)
      INTEGER realflag,start,chanend,bestend,plotflag,lowerlim,upperlim,lword
C      INTEGER channel(MaxDim),bestend,plotflag,lowerlim,upperlim,lword
      INTEGER :: channel(MaxDim)
      INTEGER :: notok(10)
C      DIMENSION xdata(MaxDim),ydata(MaxDim)
      REAL*4 :: xdata(MaxDim)
      REAL*4 :: ydata(MaxDim)
      REAL*4 :: dataz(MaxDim,MaxDim)
C      DIMENSION param(MaxDim,5),coeff(61)
      REAL*4 :: param(MaxDim,5)
      REAL*4 :: coeff(61)
C      DIMENSION smooth(MaxDim),smoothed(MaxDim,MaxDim)
      REAL*4 :: smooth(MaxDim)
      REAL*4 :: smoothed(MaxDim,MaxDim)
C      DIMENSION yout(MaxDim),xout(MaxDim),radgyr(MaxDim)
      REAL*4 :: yout(MaxDim)
      REAL*4 :: xout(MaxDim)
      REAL*4 :: radgyr(MaxDim)
C      DIMENSION limit1(MaxDim),limit2(MaxDim)
      INTEGER :: limit1(MaxDim)
      INTEGER :: limit2(MaxDim)
C      REAL asciidata(MaxDim*MaxDim)
      REAL*4 :: asciidata(MaxDim*MaxDim)

C     Array "param"
C     Param(frame no., param no.)
C     Param 1 ~ Backgr.
C     Param 2 ~ porod
C     Param 3 ~ sigma
C     Param 4 ~ A/H1
C     Param 5 ~ B/H2

1000  FORMAT(A1)
1010  FORMAT(A80)
1020  FORMAT(2x,I3)
1030  FORMAT(10I8)
1040  FORMAT(A40)
1050  FORMAT(2x,I3,2x,I3,2x,I3,2x,I3)
1060  FORMAT(E12.6,2x,E12.6,2x,E12.6,2x,I3)
1070  FORMAT(/,1x,'100: Back extrap. results:')
1080  FORMAT(1x,'100: Frame: ',I3,'. Start: '
     &,I3,'. End: ',I3,'. A = ',E12.6,'. B = ',E12.6,'.')
1090  FORMAT(1x,'100: Start channel: ',I3,'. End channel: ',
     &I3,'. A = ',E12.6,'. B = ',E12.6,'.')
1100  FORMAT(1x,'100: Frame: ',I3,'. Start: ',I3,'. End: '
     &,I3,'. H1 = ',E12.6,'. H2 = ',E12.6,'.')
1110  FORMAT(1x,'100: Start channel: ',I3,'. End channel: ',
     &I3,'. H1 = ',E12.6,'. H2 = ',E12.6,'.')
1120  FORMAT(/,1x,'100: Smoothing to frame: ',I3,'.',//,1x,
     &'   Q value      Intensity   Smoothed Int  Data',/)
1130  FORMAT(1x,'100: ',E12.6,2x,E12.6,2x,E12.6,2x,A4)
1140  FORMAT(1x,'100: ',E12.6,2x,E12.6,2x,12x,2x,A4)
1150  FORMAT(10I8)
1160  FORMAT(A10)
1170  FORMAT(E12.6,2x,E12.6,2x,E12.6,2x,E12.6,2x,E12.6)
1180  FORMAT(1x,'ERROR: Guinier back extrapolation failed'
     &//' on frame ',I3,': FATAL')
1190  FORMAT(1x,'100: Frame ',I3,': Guinier radius of gyration = ',
     &E12.6,' Angstroms.')
1200  FORMAT(/,1x,'ERROR: Guinier back extrapolation failed!: FATAL')
1210  FORMAT(/,1x,'100: Guinier radius of gyration = '
     &,E12.6,' Angstroms.')
1220  FORMAT(A10)
1230  FORMAT(/,1x,'WARNING: Problem smoothing! Background exceeded ',
     &'experimental intensity for ',I3,' points.')
1240  FORMAT(2x,I3,2x,I3,2x,I3)
1250  FORMAT(2x,I3,2x,I3)
1260  FORMAT(1x,"100: Written extrapolated data: ",A10)
1270  FORMAT(1x,"100: Written extrapolated data q-axis: ",A10)
1280  FORMAT(1x,"100: Written Guinier radius of gyration v frame: ",A10)
1285  FORMAT(1x,"100: Written ASCII Guinier radius of gyration v frame: ",A10)
1290  FORMAT(1x,"100:")
C     S King, July 2004
1300  FORMAT(1x,"100: Written extrapolated ASCII q-axis data: ",A10)
1310  FORMAT(1x,"100: Written extrapolated ASCII data: ",A10)

C     Set up promptss, titles...
      prompts(1)='ERROR: Error reading parameter files: FATAL'
      prompts(2)='ERROR: Error reading otoko data files: FATAL'
      prompts(3)='ERROR: Expecting static Q axis,'
     &//' received dynamic: FATAL'
      prompts(4)='All necessary files correctly loaded...'
      prompts(5)='ERROR: Back extrapolate algorithm has failed '
     &//'to find start of genuine data: FATAL'
      prompts(6)='ERROR: Back extrapolate channel optimisation '
     &//'has failed: FATAL'
      prompts(7)='ERROR: Problem with lower tail limit: FATAL'
      prompts(8)='WARNING: Problem smoothing: background intensity '
     &//'greater than an expt. data point.'
      prompts(9)='Smoothing data: please wait...'
      prompts(10)='Smoothing complete...'
      prompts(11)='Creating extrapolated data files...'
      prompts(12)='ERROR: Software error with subroutine'
     &//' changeotok: FATAL'
      prompts(13)='ERROR: Error writing otoko format extrapolated '
     &//'data files: FATAL'
      prompts(14)='ERROR: Error with Q axis: FATAL'
      prompts(15)='Sorry: no graphics for realtime data...'
      prompts(16)='DISPLAYING GRAPHICS: PRESS MIDDLE MOUSE BUTTON '
     &//'TO CONTINUE...'
      prompts(17)='ERROR: Error writing file tailinfo.dat: FATAL'
      prompts(18)='100: Guinier radii of gyration [Angstroms]:'
      prompts(19)='ERROR: Software error with subroutine'
     &//' changeotok: FATAL'
      prompts(20)='ERROR: Error writing otoko file for graphics '
     &//'output: FATAL'
      prompts(21)='ERROR: Error with status file: FATAL'
      prompts(22)='ERROR: Data contains negative intensities: FATAL'
      prompts(23)='ERROR: Error with limitinfo.dat: FATAL'
      prompts(24)='Finished tailjoin'
      prompts(25)='100: EXTRAPOLATING TO ZERO AND INFINITY...'
C     S King, Jul 2004
      prompts(26)='WARNING: Problem writing the extrapolation in ASCII'
      prompts(27)=
     &'WARNING: Found negative intensities'
      prompts(28)='WARNING: Error writing ascii output'


      CALL WRDLEN(LWORD)
      realflag=1
      if(static_bn)then
        realflag=0
        realtime(1)=1
        realtime(2)=1
        realtime(3)=1
        realtime(4)=1
      endif

C     Initialise
      DO i=1,MaxDim
        radgyr(i)=0.
      END DO

C      WRITE(6,*)
C      title='Tailjoin: extrapolation to zero and infinity'
C      CALL showtitle(title)
      WRITE(6,1290)
      CALL showprompts(prompts(25))

      fname=filename
      axisname=qaxname
      CALL getpath(fname,dirname)

c     Read the data
      CALL getfiletype(filetype,fname)
      IF(filetype.EQ."ascii")THEN
          CALL readascii(fname,notok,asciidata,irc)
          IF(irc.EQ.0)THEN
              GOTO 5020
          ENDIF
          nndata=ndata*4
          ij=1
          DO nframe=realtime(1),realtime(2),realtime(3)
            DO i=1,ndata
              dataz(nframe,i)=asciidata(ij)
              ij=ij+1
            END DO
          END DO
      ELSE
C         First read intensity header.
          OPEN(UNIT=9,FILE=fname,STATUS='old',ERR=5020)
          READ(9,1000,ERR=5020)letter
          READ(9,1000,ERR=5020)letter
          READ(9,1030,ERR=5020)notok
          READ(9,1040,ERR=5020)othername
          CLOSE(9)

          nndata=4*ndata

C         READ intensities
          CALL addstrings(dirname,othername,fname2)
          OPEN(UNIT=9,FILE=fname2,STATUS='old',
     &    ACCESS='direct',RECL=nndata/lword,ERR=5020)
          DO nframe=realtime(1),realtime(2),realtime(3)
            READ(9,REC=nframe,ERR=5020)(dataz(nframe,i),i=1,ndata)
          END DO
          CLOSE(9)
      ENDIF

      IF(filetype.NE."ascii")THEN
C SMK   CALL swap(data,notok(4),512*512,4)
C        CALL swap(data,notok(4),MaxDim*MaxDim,4)
        CALL swap(dataz,notok(4),MaxDim*MaxDim,4)
      ENDIF

      nframe=1

C     BACKGROUND SUBTRACT!
      DO nframe=realtime(1),realtime(2),realtime(3)
        DO i=1,MaxDim
          dataz(nframe,i)=dataz(nframe,i)-param(nframe,1)
        END DO
      END DO

C     Open Q axis
      CALL getfiletype(filetype,axisname)
      IF(filetype.EQ."ascii")THEN
          CALL readascii(axisname,notok,asciidata,irc)
          IF(irc.EQ.0)THEN
              GOTO 5020
          ENDIF
C         Check static
          IF (notok(2) .NE. 1) THEN
            CALL showprompts(prompts(3))
            STOP
          ENDIF
          DO i=1,ndata
            xdata(i)=asciidata(i)
          END DO
      ELSE
          OPEN(UNIT=9,FILE=axisname,STATUS='old',ERR=5020)
          READ(9,1000,ERR=5020)letter
          READ(9,1000,ERR=5020)letter
          READ(9,1030,ERR=5020)notok
          READ(9,1040,ERR=5020)othername
          CLOSE(9)

C         Check static
          IF (notok(2) .NE. 1) THEN
            CALL showprompts(prompts(3))
            STOP
          ENDIF

C         Read x axis data
          CALL addstrings(dirname,othername,fname2)
          OPEN(UNIT=9,FILE=fname2,STATUS='old',
     &    ACCESS='direct',RECL=nndata/lword,ERR=5020)
          READ(9,REC=1,ERR=5020)(xdata(i),i=1,ndata)
2010      CLOSE(9)
      ENDIF

      IF(filetype.NE."ascii")THEN
C SMK   CALL swap(xdata,notok(4),512,4)
        CALL swap(xdata,notok(4),MaxDim,4)
      ENDIF

C     Ok - data fine
      CALL showprompts(prompts(4))


C     START NUMBER CRUNCH


C     Loop through frames
      WRITE(6,1070)
      DO nframe=realtime(1),realtime(2),realtime(3)

C       Load up data
        DO i=1,MaxDim
          ydata(i)=dataz(nframe,i)
        END DO

C       Back extrapolate.
C       Got our starting point from corfunc.dat.
C       Skip three points ~ beamstop scatter?
C       No points skipped! Beam up-turn is too small.
2020    start=datastart

C       Optimise end point.
C       Set up initial situation
        oldchisqu=1.0E+20
        h1store=0.
        h2store=0.
        astore=0.
        bstore=0.
        bestend=0

C       Loop through end points
        DO chanend=start+4,start+15

          IF (backex .EQ. 'vonk' .OR.
     &    backex .EQ. 'VONK') THEN
C           Vonk model
            sigy=0.
            sigxsqu=0.
            sigxsquy=0.
            sigx4=0.

            DO j=start,chanend
C             Get sums
              temp=xdata(j)**2
              sigy=sigy+ydata(j)
              sigxsqu=sigxsqu+temp
              sigxsquy=sigxsquy+temp*ydata(j)
              sigx4=sigx4+temp*temp
            END DO

C           Calculate parameters
            npts=(chanend-start+1)
            h2=(sigxsqu*sigy-npts*sigxsquy)/(npts*sigx4-(sigxsqu**2))
            h1=(sigy+h2*sigxsqu)/npts

C           Get chi squ.
            chisqu=0.
            DO j=start,chanend
              chisqu=chisqu+(ydata(j)-h1+h2*xdata(j)*xdata(j))**2
            END DO
            chisqu=chisqu/npts

          ELSE
C           Guinier model
            siglgy=0.
            sigxsqu=0.
            sigxsqulgy=0.
            sigx4=0.

            DO j=start,chanend
C             Get sums
              temp1=LOG(ydata(j))
              temp2=xdata(j)**2
              siglgy=siglgy+temp1
              sigxsqu=sigxsqu+temp2
              sigxsqulgy=sigxsqulgy+temp1*temp2
              sigx4=sigx4+temp2*temp2
            END DO

C           Calculate parameters
            npts=(chanend-start+1)
            b=(sigxsqu*siglgy-npts*sigxsqulgy)/((sigxsqu**2)-npts*sigx4)
            alg=(siglgy-b*sigxsqu)/npts

C           Get chi squ.
            a=EXP(alg)
            chisqu=0.
            DO j=start,chanend
              chisqu=chisqu+(ydata(j)-a*EXP(b*xdata(j)*xdata(j)))**2
            END DO
            chisqu=chisqu/npts
          ENDIF

C         Debug:
C           write(6,*)start,chanend,a,b,chisqu

C         Now compare chi squareds.
          IF (chisqu .LT. oldchisqu) THEN
            oldchisqu=chisqu
            h1store=h1
            h2store=h2
            astore=EXP(alg)
            bstore=b
            bestend=chanend
          ENDIF

C       End of "chanend" optimisation loop.
        END DO

C       Check it's worked
        IF (bestend .EQ. 0) THEN
          CALL showprompts(prompts(6))
          STOP
        ENDIF

C       Restore params.
        a=astore
        b=bstore
        h1=h1store
        h2=h2store

C       Output results.
        IF (backex .EQ. 'VONK' .OR.
     &  backex .EQ. 'vonk') THEN
          IF (realflag .EQ. 1) THEN
            WRITE(6,1100)nframe,start,bestend,h1,h2
          ELSE
            WRITE(6,1110)start,bestend,h1,h2
          ENDIF
          param(nframe,4)=h1
          param(nframe,5)=h2
        ELSE
          IF (realflag .EQ. 1) THEN
            WRITE(6,1080)nframe,start,bestend,a,b
          ELSE
            WRITE(6,1090)start,bestend,a,b
          ENDIF
          param(nframe,4)=a
          param(nframe,5)=b
        ENDIF

C       Back extrapolate over.

C     End loop through frames
      END DO

C     Output radii of gyration.
      IF (backex .EQ. 'guinier' .OR. backex .EQ. 'GUINIER') THEN
C       check realtime
        IF (realtime(4) .EQ. 1) THEN
C         Single frame
          nframe=realtime(1)
          IF (param(nframe,5) .GT. 0.) THEN
            WRITE(6,1200)
          ELSE
            rad=SQRT(-3.*param(nframe,5))
            WRITE(6,1210)rad
          ENDIF

        ELSE
C         More than one frame
          WRITE(6,*)
          CALL showprompts(prompts(18))
          WRITE(6,*)
          DO nframe=realtime(1),realtime(2),realtime(3)
            IF (param(nframe,5) .GT. 0.) THEN
              WRITE(6,1180)nframe
              radgyr(nframe)=0.
            ELSE
              rad=SQRT(-3.*param(nframe,5))
              WRITE(6,1190)nframe,rad
              radgyr(nframe)=rad
            ENDIF
          END DO
        ENDIF
      ENDIF


C     Start Savitsky Golay Smoothing.
C     The subroutines savgol ludcmp and lubksb are from numerical
C     recipes and simply calculate a set of numbers indept of the data.
C     Hence, once debugging is over, the smoothing routine will not
C     call these subroutines, but use a permament set of values.
C
C     Of course, once preliminary debugging is over, the only way
C     of telling whether the smoothing has "worked" is to perform
C     the Fourier transform and inspect that.
C
C     Subroutines debugged using results provided by Numerical Recipes.
C
C     Present Sav-Gol setup: 4th order polynomial,
C     10 points either side of window.
C
C     Smoothing performed in LOG(intensity-background) world.
C     ALTERATION!!!
C     Smoothing performed in LOG(intensity) world.

C     Tell user we're smoothing
      CALL showprompts(prompts(9))

      norder=4
      nsize=10
C     norder and nsize control all aspects of smoothing.
C     Max nsize=30.
      nl=nsize
      nr=nsize
      np=61
      CALL savgol(coeff,np,nl,nr,0,norder)
C     Debug: textoutput of coefficients and check sum
C      sum=0.
C      DO i=1,11
C        WRITE(6,*)coeff(i)
C        sum=sum+coeff(i)
C      END DO
C      sum=2.*(sum-coeff(1))+coeff(1)
C      write(6,*)sum


C     Loop through frames
      DO nframe=realtime(1),realtime(2),realtime(3)

C       Set up join channel limits
        lim1=channel(nframe)
        join=lim1+nsize
        IF ((lim1 .LT. qzero+nsize+1) .OR.
     &  (lim1 .GT. ndata+3*nsize)) THEN
          CALL showprompts(prompts(7))
          STOP
        ENDIF

C       Load expt. intensities
        DO i=lim1-nsize,lim1+nsize-1
C         Add on backgr.
          smooth(i)=dataz(nframe,i)+param(nframe,1)
C         Check for -ve log
          IF (smooth(i) .GT. 0.) THEN
            smooth(i)=LOG(smooth(i))
          ELSE
C           Fatal error: we can't be having negative intensities can we?
C           S King, July 2004
C            CALL showprompts(prompts(22))
C            STOP
            call showprompts(prompts(27))
C SMK           smooth(i)=LOG(1.0E-08)
            if (i.eq.1) then
              smooth(i)=-18.4
            else
              smooth(i)=smooth(i-1)
            end if
          ENDIF
        END DO

C       Load calc. intensities.
        DO i=lim1+nsize,lim1+(3*nsize)-1
          x=xdata(i)
          smooth(i)=((param(nframe,2)/(x**4))*
     &    EXP(-(param(nframe,3)**2)*(x**2)))
          smooth(i)=LOG(smooth(i)+param(nframe,1))
        END DO

C       Now smooth according to coefficients in "coeff".
        DO i=lim1,lim1+(2*nsize)-1
          sum=0.
          DO j=-nl,nr,1
            IF (j .LE. 0) THEN
              wt=coeff(1-j)
            ELSE
              wt=coeff(np+1-j)
            ENDIF
            sum=sum+wt*smooth(i+j)
          END DO
          smoothed(nframe,i)=EXP(sum)-param(nframe,1)
        END DO

C       Debug: text output
        noutflag=0
        IF (noutflag .EQ. 1) THEN
          WRITE(6,1120)nframe
          statuss='Expt'
          DO i=lim1-nsize,lim1-1
            y=dataz(nframe,i)
            WRITE(6,1140)xdata(i),y,statuss
          END DO
          DO i=lim1,lim1+(2*nsize)-1
            IF (i .LT. join) THEN
              statuss='Expt'
            ELSE
              statuss='Calc'
            ENDIF
            y=dataz(nframe,i)
            WRITE(6,1130)xdata(i),y,
     &      smoothed(nframe,i),statuss
          END DO
          statuss='Calc'
          DO i=lim1+(2*nsize),lim1+(3*nsize)-1
            x=xdata(i)
            y=((param(nframe,2)/(x**4))*
     &      EXP(-(param(nframe,3)**2)*(x**2)))
            WRITE(6,1140)x,y,statuss
          END DO
        ENDIF

C       Update data
        DO i=lim1,lim1+(2*nsize)-1
          dataz(nframe,i)=smoothed(nframe,i)
        END DO

C     End loop through frames
      END DO



C     END NUMBER CRUNCH


C     Tell user it's over
      CALL showprompts(prompts(10))

C     Rebuild: output extrapolated data to otoko file.
      CALL showprompts(prompts(11))

C     First build an x axis.
      DO i=1,ndata-qzero+1
        xout(i)=xdata(qzero+i-1)
      END DO
      IF(xdata(ndata).GE. 10.0)THEN
          CALL showprompts(prompts(14))
          STOP
      ENDIF
      deltaq=( (10.0 - xdata(ndata)) / FLOAT(MaxDim-ndata+qzero-1) )
      DO i=ndata-qzero+2,MaxDim
        xout(i)=xout(i-1)+deltaq
      END DO

C     Not sure what the record length should be: try -
      nrecl=8192

C     Open x axis otoko header.
      call strippath(qaxname,xheader)
      CALL swapexten(xheader,'FLX')
C     Change header name to give data filename
      CALL changeotok(xheader,othername,nerr)
      IF (nerr .EQ. 1) THEN
        CALL showprompts(prompts(12))
        STOP
      ENDIF
      fname=xheader
      OPEN(UNIT=9,FILE=fname,STATUS='unknown',ERR=5030)
      WRITE(9,*,ERR=5030)
      WRITE(9,*,ERR=5030)
      notok(1)=MaxDim
      notok(2)=1
      notok(3)=1
      CALL endian(notok(4))
      DO i=5,10
        notok(i)=0
      END DO
      WRITE(9,1150,ERR=5030)notok
      WRITE(9,1160,ERR=5030)othername(1:10)
      CLOSE(9)
C     X header done

C     Now do x data
C     Is RECL correct?
      fname=othername
      OPEN(UNIT=10,FILE=fname,STATUS='unknown',
     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)
      WRITE(10,REC=1,ERR=5030)
     &(xout(i),i=1,MaxDim)
      CLOSE(10)


C     Rebuild y data.
C     First open header.
      CALL strippath(filename,yheader)
      CALL swapexten(yheader,'FUL')
C     Change header name to give data filename
      CALL changeotok(yheader,othername,nerr)
      IF (nerr .EQ. 1) THEN
        CALL showprompts(prompts(12))
        STOP
      ENDIF
      fname=yheader
      yheader2=yheader(1:9)//'2'
      fname2=yheader2
      OPEN(UNIT=9,FILE=fname,STATUS='unknown',ERR=5030)
      OPEN(UNIT=11,FILE=fname2,STATUS='unknown',ERR=5030)
      WRITE(9,*,ERR=5030)
      WRITE(11,*,ERR=5030)
      WRITE(9,*,ERR=5030)
      WRITE(11,*,ERR=5030)
      notok(2)=realtime(4)
      notok(1)=MaxDim
      notok(3)=1
      CALL endian(notok(4))
      DO i=5,10
        notok(i)=0
      END DO
      WRITE(9,1150,ERR=5030)notok
      WRITE(11,1150,ERR=5030)notok
      WRITE(9,1160,ERR=5030)othername(1:10)
      WRITE(11,1160,ERR=5030)othername(1:9)//'2'
      CLOSE(9)
      CLOSE(11)
C     Y header done

C     Calculate and output Y data.
      fname=othername
      fname2=fname(1:9)//'2'
      OPEN(UNIT=10,FILE=fname,STATUS='unknown',
     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)
      OPEN(UNIT=11,FILE=fname2,STATUS='unknown',
     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)
      nframe=0
      DO nfr=realtime(1),realtime(2),realtime(3)
        lim1=channel(nfr)
        nframe=nframe+1
C       Calculate y data depending on region.
        DO i=1,MaxDim
          x=xout(i)
C         Regions -
          IF (i .LE. (datastart-qzero)) THEN

C           Build data using back extrapolate.
            IF (backex .EQ. 'VONK' .OR.
     &      backex .EQ. 'vonk') THEN
C             Vonk model
              y=param(nfr,4)-param(nfr,5)*x*x
            ELSE
C             Guinier model
              y=param(nfr,4)*EXP(param(nfr,5)*x*x)
            ENDIF

C           Debug test for F transform
C           y=0.

          ELSEIF (i .GT. lim1+20-qzero) THEN

C           Build data using sigmoid tail
            y=((param(nfr,2)/(x**4))*
     &      EXP(-(param(nfr,3)**2)*(x**2)))

          ELSE

C           Output expt. or smoothed data.
            y=dataz(nfr,i-1+qzero)

          ENDIF

C         Store:
          yout(i)=y
C         End rebuilding a frame
        END DO

C       Write to file
        WRITE(10,REC=nframe,ERR=5030)(yout(i),i=1,MaxDim)
        WRITE(11,REC=nframe,ERR=5030)(yout(i)+param(nframe,1),i=1,MaxDim)

C     End realtime loop.
      END DO
      CLOSE(10)
      CLOSE(11)

C     Otoko file output complete.

      WRITE(6,1260)yheader
      WRITE(6,1260)yheader2
      WRITE(6,1270)xheader

C     Change tailinfo.dat to include backextrap. params
C      OPEN(UNIT=9,FILE='tailinfo.dat',STATUS='unknown',ERR=5040)
C      DO nframe=realtime(1),realtime(2),realtime(3)
C        WRITE(9,1170,ERR=5040)param(nframe,1),param(nframe,2),
C     &  param(nframe,3),param(nframe,4),param(nframe,5)
C      END DO
C      CLOSE(9)

C     Also pass info to limitinfo.dat
C      OPEN(UNIT=9,FILE='limitinfo.dat',STATUS='unknown',ERR=5070)
C      DO nframe=realtime(1),realtime(2),realtime(3)
C        WRITE(9,1240,ERR=5070)limit1(nframe),limit2(nframe),
C     &  channel(nframe)
C      END DO
C      CLOSE(9)

C       Realtime output.
C       If guinier model create an otoko file containing radii of gyration.

      IF(realtime(4).NE.1) THEN
        IF (backex .EQ. 'guinier' .OR. backex .EQ. 'GUINIER') THEN
C         Get filename.
          CALL strippath(filename,radname)
          radax=radname
          CALL swapexten(radname,'RAD')
          CALL swapexten(radax,'FAX')
C         x-axis should already exist

C         Create otoko header
C         Change header name to give data filename
          CALL changeotok(radname,othername,nerr)
          IF (nerr .EQ. 1) THEN
            CALL showprompts(prompts(19))
          STOP
          ENDIF
          fname=radname
          OPEN(UNIT=9,FILE=fname,STATUS='unknown',ERR=5050)
          WRITE(9,*,ERR=5050)
          WRITE(9,*,ERR=5050)
          notok(1)=realtime(4)
          notok(2)=1
          notok(3)=1
          CALL endian(notok(4))
          DO i=5,10
            notok(i)=0
          END DO
          WRITE(9,1030,ERR=5050)notok
          WRITE(9,1220,ERR=5050)othername(1:10)
          CLOSE(9)
C         Header done

          nrecl=4*realtime(4)
          outname=othername
          OPEN(UNIT=10,FILE=outname,STATUS='unknown',
     &    ACCESS='direct',RECL=nrecl/lword,ERR=5050)
          WRITE(10,REC=1,ERR=5050)
     &    (radgyr(i),i=realtime(1),realtime(2),realtime(3))
          CLOSE(10)

          write(6,1280)fname

C         Output the radius of gyration data in ASCII format
          CALL writeascii(fname(1:3)//"RAD."//"TXT",notok,radgyr,
     &                  irc,1)
          IF(irc.EQ.0)THEN
            CALL showprompts(prompts(28))
          ELSE
            WRITE(6,1285)fname(1:3)//"RAD."//"TXT"
          ENDIF

        ENDIF
      ENDIF


C     S King, July 2004
C     Output in ASCII format
      open(12,file=fname(1:3)//'FLX.TXT',form='formatted',
     &status='unknown',err=5060)
        write(12,'(f12.5)',ERR=5060)(xout(i),i=1,MaxDim)
      close(12)
      write(6,1300)fname(1:3)//'FLX.TXT'

      open(12,file=fname(1:3)//'FUL.TXT',form='formatted',
     &status='unknown',err=5060)
        write(12,'(e12.5)',ERR=5060)(yout(i),i=1,MaxDim)
      close(12)
      write(6,1310)fname(1:3)//'FUL.TXT'

      open(12,file=fname(1:3)//'FU2.TXT',form='formatted',
     &status='unknown',err=5060)
        write(12,'(e12.5)',ERR=5060)(yout(i)+param(nframe,1),i=1,MaxDim)
      close(12)
      write(6,1310)fname(1:3)//'FU2.TXT'

C     End program.
      CALL showprompts(prompts(24))
      RETURN
      STOP


C     Error messages

C     Static data.
5010  realflag=0
      realtime(1)=1
      realtime(2)=1
      realtime(3)=1
      realtime(4)=1

C     Error reading otoko intensity or q data.
5020  CALL showprompts(prompts(2))
      STOP

C     Error writing otoko format rebuild files.
5030  CALL showprompts(prompts(13))
      STOP

C     Error writing radius of gyration output
5050  CALL showprompts(prompts(20))
      STOP

C     Error writing ASCII intensity or q data.
5060  CALL showprompts(prompts(26))

      RETURN
      END



      SUBROUTINE savgol(c,np,nl,nr,ld,m)
C     Savitsky Golay smoothing.
C     Taken from Numerical Recipes, Flannery and Press.
C     Calculates coefficients for linear smoothing,
C     such that higher moments are preserved.
C     Re. correlation function analysis, smoothing used
C     to ease the join between calc. and expt. data.
      INTEGER ld,m,nl,np,nr,mmax
C      REAL c(np)
      REAL*4 :: c(np)
      PARAMETER (mmax=6)
C      INTEGER imj,ipj,j,k,kk,mm,indx(mmax+1)
      INTEGER imj,ipj,j,k,kk,mm
      INTEGER :: indx(mmax+1)
C      REAL d,fac,sum,a(mmax+1,mmax+1),b(mmax+1)
      REAL d,fac,sums
      REAL*4 :: a(mmax+1,mmax+1)
      REAL*4 :: b(mmax+1)


      IF (np.LT.nl+nr+1 .OR. nl.LT.0 .OR. nr.LT.0 .OR. ld.GT.m .OR.
     &m.GT.mmax .OR. nl+nr.LT.m) then
          write(6,*) 'Bad arguments passed to subroutine savgol'
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
C      INTEGER n,np,indx(n),nmax
      INTEGER n,np,nmax
      INTEGER :: indx(n)
C      REAL d,a(np,np),tiny
      REAL d,teensy
      REAL*4 :: a(np,np)
      PARAMETER (nmax=500,teensy=1.E-8)
      INTEGER i,imax,j,k
C      REAL aamax,dum,sum,vv(nmax)
      REAL aamax,dum,sums
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
C      INTEGER n,np,indx(n)
      INTEGER n,np
      INTEGER :: indx(n)
C      REAL a(np,np),b(n)
      REAL*4 :: a(np,np)
      REAL*4 :: b(n)
      INTEGER i,ii,j,ll
      REAL sums
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

	  
C=============================================================================
C==============================================================================


      SUBROUTINE extrapolate
C     TMWN July 94
C     Gets information from user for tailfit and tailjoin subroutines.
C     SMK Nov 05
C     Redimensioned data arrays from 512 to MaxDim.

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

      EXTERNAL corfunc_io, wrdlen, getfiletype, readascii, getpath
      EXTERNAL addstr, byte_order, mrqmin, strpath, tred2, tqli
      EXTERNAL syminv, covsrt, tbksb, writeascii, loqread

      CHARACTER*80 filename,qaxname,dirname,exptname,stats
      CHARACTER*40 title,sigmodel,ascii,user,idftoggle,backex,retrans
      CHARACTER*80 :: prompts(44)
      CHARACTER*80 fname,axisname,fname2,arbabs
      CHARACTER*10 filetype
      CHARACTER*1 letter
C      INTEGER lim1,lim2,qzero,realtime(5),datastart,lword,best(MaxDim)
      INTEGER lim1,lim2,qzero,datastart,lword
      INTEGER :: best(MaxDim)
      INTEGER :: realtime(5)
C      DIMENSION notok(10)
      INTEGER :: notok(10)
C      DIMENSION xdata(MaxDim),ydata(MaxDim,MaxDim),param2(MaxDim,5)
      DIMENSION xdata(MaxDim)
      DIMENSION ydata(MaxDim,MaxDim)
      DIMENSION param2(MaxDim,5)
C      DIMENSION limit1(MaxDim),limit2(MaxDim),results(MaxDim,5),param(3)
      INTEGER :: limit1(MaxDim)
      INTEGER :: limit2(MaxDim)
      REAL*4 :: results(MaxDim,5)
      REAL*4 :: param(3)
      LOGICAL static_bn
C      INTEGER frames(4),bestchannel,channel2(MaxDim,2)
      INTEGER :: frames(4)
      INTEGER bestchannel
      INTEGER :: channel2(MaxDim,2)
C      REAL asciidata(MaxDim*MaxDim)
      REAL*4 :: asciidata(MaxDim*MaxDim)
C     S King, July 2004
      character*3 asciitype
      character*80 tit1,tit2
C      dimension c(MaxDim),e(MaxDim)
      real*4 :: c(MaxDim)
      real*4 :: e(MaxDim)

1000  FORMAT(A1)
1010  FORMAT(A80)
1020  FORMAT(2x,I3)
1030  FORMAT(10I8)
1040  FORMAT(A40)
1050  FORMAT(/,1x,'100: Q limits of tail: start ',E12.6,'.',/,
     &' 100: ',18x,'end   ',E12.6,'.')
1060  FORMAT(/,1x,'Data contains ',I3,' frames.')
1070  FORMAT(/,1x,'100: Analysis will be performed on ',I3,' frames.')
1080  FORMAT(2x,I3,2x,I3,2x,I3,2x,I3)
1090  FORMAT(E12.6)
1100  FORMAT(/,1x,'100: Limits for frame ',I3,':')
1110  FORMAT(1x,'100: Q start ',E12.6,'. Q end ',E12.6,'.')
1120  FORMAT(2x,I3,2x,I3)
1130  FORMAT(E12.6,2x,E12.6,2x,E12.6,2x,E12.6,2x,E12.6)
1140  FORMAT(2x,I3,2x,I3)

      title='Tail Fit User Input'
      CALL showtitle(title)
      prompts(1)='Enter data directory name'
      prompts(2)='Enter intensity data filename'
      prompts(3)='Enter Q axis filename'
      prompts(4)='ERROR: Error reading corfunc.txt: NON-FATAL'
      prompts(5)='ERROR: Problem with data: FATAL'
      prompts(6)='ERROR: X axis should be a static image, '
     &//'dynamic found: FATAL'
      prompts(7)='Data OK...'
      prompts(8)='Enter channel at start of tail'
      prompts(9)='Enter channel at end of tail'
      prompts(10)='Enter start channel optimisation range '
     &//'(type n for no opt.)'
      prompts(11)='Enter number of iterations for tailfit'
      prompts(12)='Do you want graphics [y/n]'
      prompts(13)='ERROR: Invalid input: FATAL'
      prompts(14)='Passing information to program tailfit...'
      prompts(15)='Data in arbitrary or absolute units [arb/abs]'
      prompts(16)='ERROR: Error writing corfunc.txt: FATAL'
      prompts(17)='ERROR: Error writing dummy file: FATAL'
      prompts(18)='CORFUNC.TXT: PARAMETER FILE USED BY '
     &//'CORRELATION PROGRAMS. AVOID DELETION.'
      prompts(19)='Enter start frame'
      prompts(20)='Enter end frame'
      prompts(21)='Enter increment'
      prompts(22)='static'
      prompts(23)='Vonk or Guinier model for back-extrapolation [v/g]'
      prompts(24)='Enter channel at start of genuine data'
      prompts(25)='Enter maximum D for Fourier transform [Angstroms]'
      prompts(26)='Enter step in D for Fourier transform [Angstroms]'
      prompts(27)='Re-transform correlation function [y/n]'
      prompts(28)='Apply sigmoid tail model (otherwise Porod) [y/n]'
      prompts(29)='Enter estimate of volume fraction crystallinity'
      prompts(30)='Do you want ascii output of '
     &//'correlation function [y/n]'
      prompts(31)='ERROR: Error with status file: FATAL'
      prompts(32)='Do you want to continue with this set-up [y/n]'
      prompts(33)='Do you want user control of extraction '
     &//'process [y/n]'
      prompts(34)='Same channel limits on each frame for tailfit [y/n]'
      prompts(35)='Enter tailfit start channel'
      prompts(36)='Enter tailfit end channel'
      prompts(37)='ERROR: Error writing limitinfo.txt: FATAL'
      prompts(38)='Calculate the interface distribution function [y/n]'
      prompts(39)='Start transform? [y/n]'
      prompts(40)='Start extraction? [y/n]'
C     Following added by S M King, July 2004
      prompts(41)='ERROR: READASCII_ returning IRC=0: FATAL'
      prompts(42)='ERROR: Error during file OPEN: FATAL'
      prompts(43)='ERROR: Error during file READ: FATAL'
      prompts(44)='ERROR: X axis has more/fewer points than data: FATAL'

      filename='A00000.XSH'
      qaxname='X00000.QAX'
      lim1=200
      lim2=400
      CALL WRDLEN(lword)

C     Got default info. Now quiz user.
      CALL defaultname(prompts(2),filename)
      CALL defaultname(prompts(3),qaxname)

C     Check files exist
      fname=filename
      axisname=qaxname

C     Open data
      CALL getfiletype(filetype,fname)
      IF(filetype.EQ."ascii")THEN
C       S King, July 2004
        asciitype='sca'
C       Is it a LOQ 1D file...
        irc=0
        open(12,file=fname,form='formatted',status='old',err=10)
C SMK         call loqread(tit1,tit2,512,ndata,xdata,c,e,12,irc)
         call loqread(tit1,tit2,MaxDim,ndata,xdata,c,e,12,irc)
10      close(12)
        if(irc.eq.1)then
          notok(1)=ndata
          notok(2)=1
          nndata=ndata*4
C         Assume single frame only (ie, not realtime)
         realtime(1)=notok(2)
          do i=1,ndata
            ydata(1,i)=c(i)
          end do
          asciitype='loq'
        else
C       ...or is it single column ascii (sca)?
          CALL readascii(fname,notok,asciidata,irc)
          IF(irc.EQ.0)THEN
              GOTO 5021
          ENDIF
          ndata=notok(1)
          nndata=ndata*4
          realtime(1)=notok(2)
          ij=1
          DO i=1,ndata
            DO nframe=1,realtime(1)
              ydata(nframe,i)=asciidata(ij)
              ij=ij+1
            END DO
          END DO
        endif
      ELSE
C SMK The program falls over at this next open under Cygwin... use MinGW instead!
C     The problem appears to be related to how the fname character string is handled in the different flavours
C     of Unix.
          OPEN(UNIT=9,FILE=fname,STATUS='old',ERR=5022)
          READ(9,1000,ERR=5023)letter
          READ(9,1000,ERR=5023)letter
          READ(9,1030,ERR=5023)notok
          READ(9,1040,ERR=5023)exptname
          CLOSE(9)
          ndata=notok(1)
          nndata=ndata*4

C         Check realtime
C         realtime(1)=number of frames
C         realtime(2, 3, and 4) are start,end,increm
C         realtime(5)=number of frames this setup uses
          realtime(1)=notok(2)

C         Check intensities
          CALL getpath(fname,dirname)
          CALL addstrings(dirname,exptname,fname2)
          OPEN(UNIT=9,FILE=fname2,STATUS='old',
     &    ACCESS='direct',RECL=nndata/lword,ERR=5022)
C         Loop through frames
          DO nframe=1,realtime(1)
            READ(9,REC=nframe,ERR=5023)(ydata(nframe,i),i=1,ndata)
          END DO
          CLOSE(9)
      ENDIF

      IF(filetype.NE."ascii")THEN
C SMK   CALL swap(ydata,notok(4),512*512,4)
        CALL swap(ydata,notok(4),MaxDim*MaxDim,4)
      ENDIF

C     reset frame counter
      nframe=1

C     Open Q axis
C     S King, July 2004
C     Don't need to do this if read a LOQ ASCII file
      if(asciitype.eq.'loq') goto 20

      CALL getfiletype(filetype,axisname)
      IF(filetype.EQ."ascii")THEN
          CALL readascii(axisname,notok,asciidata,irc)
          IF(irc.EQ.0)THEN
              GOTO 5021
          ENDIF
          IF (ndata .NE. notok(1)) THEN
C           x axis has more/fewer points than data
            GOTO 5024
          ENDIF
C         Check realtime
          IF (notok(2) .NE. 1) THEN
            CALL showprompts(prompts(6))
            STOP
          ENDIF
          DO i=1,ndata
            xdata(i)=asciidata(i)
          END DO
      ELSE
          OPEN(UNIT=9,FILE=axisname,STATUS='old',ERR=5022)
          READ(9,1000,ERR=5023)letter
          READ(9,1000,ERR=5023)letter
          READ(9,1030,ERR=5023)notok
          READ(9,1040,ERR=5023)exptname
          CLOSE(9)

          IF (ndata .NE. notok(1)) THEN
C           x axis has more/fewer points than data
            GOTO 5024
          ENDIF
C         Check realtime
          IF (notok(2) .NE. 1) THEN
            CALL showprompts(prompts(6))
            STOP
          ENDIF

C         check x axis data
          CALL getpath(fname,dirname)
          CALL addstrings(dirname,exptname,fname2)
          OPEN(UNIT=9,FILE=fname2,STATUS='old',
     &    ACCESS='direct',RECL=nndata/lword,ERR=5022)
          READ(9,REC=nframe,ERR=5023)(xdata(i),i=1,ndata)
2015      CLOSE(9)
      ENDIF

      IF(filetype.NE."ascii")THEN
C SMK   CALL swap(xdata,notok(4),512,4)
        CALL swap(xdata,notok(4),MaxDim,4)
      ENDIF

C     OK - the data's there
20    CALL showprompts(prompts(7))

C     Set up realtime defaults
      realtime(2)=1
      realtime(3)=realtime(1)
      realtime(4)=1
      realtime(5)=1
C     Get realtime info if necessary
2050  IF (realtime(1) .NE. 1) THEN
        WRITE(6,1060)realtime(1)
C       Get start of frames
        CALL defaultint(prompts(19),realtime(2))
C       check valid value
        IF (realtime(2) .LT. 1 .OR.
     &  realtime(2) .GT. realtime(1)) THEN
          CALL showprompts(prompts(13))
          GOTO 2050
        ENDIF
C       get end of frames
        CALL defaultint(prompts(20),realtime(3))
C       check valid
        IF (realtime(3) .LT. realtime(2) .OR.
     &  realtime(3) .GT. realtime(1)) THEN
          CALL showprompts(prompts(13))
          GOTO 2050
        ENDIF

C       get frame increment
        CALL defaultint(prompts(21),realtime(4))

C       work out number of frames this increment would use
        realtime(5)=0
        DO i=realtime(2),realtime(3),realtime(4)
          realtime(5)=realtime(5)+1
        END DO

C       check valid
        IF (realtime(5) .LT. 1 .OR.
     &  realtime(5) .GT. realtime(1)) THEN
          CALL showprompts(prompts(13))
          GOTO 2050
        ENDIF
        WRITE(6,1070)realtime(5)

      ENDIF
C     OK - realtime stuff over.


C     Find qzero
      qzero=1
      DO i=1,ndata
        IF (xdata(i-1) .LT. 0. .AND. xdata(i) .GE. 0.) THEN
          qzero=i
        ENDIF
      END DO

C     next get fit limits, no. iterations
C     altered 9/8/94
      numiteration=100
      arbabs='arb'
C     Is data arb or abs?
      CALL defaultname(prompts(15),arbabs)

C     Get tail limits
C     same for each frame
      letter='y'
      IF (realtime(5) .NE. 1) THEN
        CALL defaultletter(prompts(34),letter)
      ENDIF
      IF (letter .EQ. 'y' .OR. letter .EQ. 'Y') THEN
C       Same limits for every frame.
2020    CALL defaultint(prompts(8),lim1)
        CALL defaultint(prompts(9),lim2)
        IF (lim1 .LT. 1 .OR.
     &      lim1 .GT. lim2 .OR.
     &      lim2 .GT. MaxDim) THEN
          CALL showprompts(prompts(13))
          GOTO 2020
        ENDIF
C       write into array
        DO i=realtime(2),realtime(3),realtime(4)
          limit1(i)=lim1
          limit2(i)=lim2
        END DO
        WRITE(6,1050)xdata(lim1),xdata(lim2)

      ELSE
C       Different limits.
        DO i=realtime(2),realtime(3),realtime(4)
          WRITE(6,1100)i
2100      CALL defint(prompts(35),lim1)
          CALL defint(prompts(36),lim2)
          IF (lim1 .LT. 1 .OR.
     &      lim1 .GT. lim2 .OR.
     &      lim2 .GT. MaxDim) THEN
            CALL showprompts(prompts(13))
            GOTO 2100
          ENDIF
          limit1(i)=lim1
          limit2(i)=lim2
          WRITE(6,1110)xdata(lim1),xdata(lim2)
        END DO

      ENDIF

C     number of iterations?
2040  CALL defaultint(prompts(11),numiteration)
      IF (numiteration .LT. 1 .OR. numiteration .GT. 1000) THEN
        CALL showprompts(prompts(13))
        GOTO 2040
      ENDIF
C     Sigmoid model, or just straight Porod model?
      IF (sigmodel .EQ. 'off') THEN
        letter='n'
      ELSE
        letter='y'
      ENDIF
      CALL defaultletter(prompts(28),letter)
      IF (letter .EQ. 'n' .OR. letter .EQ. 'N') THEN
        sigmodel='off'
      ELSE
        sigmodel='on'
      ENDIF

C     Data for tailjoin.
C     back extrapolation model.
      IF (backex .EQ. 'vonk') THEN
        letter='v'
      ELSE
        letter='g'
      ENDIF
      CALL defaultletter(prompts(23),letter)
      IF (letter .EQ. 'V' .OR. letter .EQ. 'v') THEN
        backex='vonk'
      ELSE
        backex='guinier'
      ENDIF

C     Find start of genuine data.
C     Start of back extrap. is a point with intensity
C     greater than 10.*intensity at qzero and such that
C     intensity is decreasing with Q.
C     NB: will fail if semi transparent beamstop used - user must take over.
      compare=10.*ydata(realtime(2),qzero+1)
      DO i=qzero,qzero+50
        IF (ydata(realtime(2),i) .GT. compare .AND.
     &  ydata(realtime(2),i) .GT. ydata(realtime(2),i+1)) THEN
          GOTO 2060
        ENDIF
      END DO
      i=0
2060  datastart=i
      CALL defaultint(prompts(24),datastart)
      IF ((datastart .LT. qzero) .OR.
     &(datastart .GT. lim1)) THEN
        CALL showprompts(prompts(13))
        GOTO 2060
      ENDIF

C     User input over
      CALL showprompts(prompts(14))

C     Set up parameters for tailfit
      if(realtime(5).eq.1)then
        static_bn=.true.
      else
        static_bn=.false.
      endif

      frames(1)=realtime(2)
      frames(2)=realtime(3)
      frames(3)=realtime(4)
      frames(4)=realtime(5)

C     Call tailfit
      CALL tailfit(filename,qaxname,sigmodel,limit1,limit2,
     &             frames,ndata,numiteration,
     &             results,best,bestchannel,param,static_bn)

C     Set up parameters for tailjoin
      if(frames(4).eq.1)then
        best(1)=bestchannel
        results(1,1)=param(1)
        results(1,2)=param(2)
        results(1,3)=param(3)
      endif

C     Call tailjoin
      CALL tailjoin(filename,qaxname,results,best,ndata,qzero,
     &              backex,datastart,frames,limit2,static_bn)

C     Write corfunc.txt file
      open(unit=9,file='corfunc.txt',STATUS='unknown',ERR=5030)
      write(9,1010,ERR=5030)filename
      write(9,1010,ERR=5030)qaxname
      write(9,1080,ERR=5030)frames(1),frames(2),frames(3),
     &                     frames(4)
      write(9,1010,ERR=5030)arbabs
      write(9,1040,ERR=5030)sigmodel
      write(9,1040,ERR=5030)backex
      write(9,1020,ERR=5030)qzero
      DO i=frames(1),frames(2),frames(3)
          write(9,1130,ERR=5030)results(i,1),results(i,2),
     &                          results(i,3),results(i,4),
     &                          results(i,5)
      END DO
      DO i=frames(1),frames(2),frames(3)
          write(9,1140,ERR=5030)limit2(i),best(i)
      END DO
      close(9)

C     End program.
      STOP

C     Error statements

C     Error with data
5020  CALL showprompts(prompts(5))
      STOP
5021  CALL showprompts(prompts(41))
      STOP
5022  CALL showprompts(prompts(42))
      STOP
5023  CALL showprompts(prompts(43))
      STOP
5024  CALL showprompts(prompts(44))
      STOP

C     Error writing corfunc.dar
5030  CALL showprompts(prompts(16))
      STOP

      END

	  
C=============================================================================
C==============================================================================


      SUBROUTINE ftransform
C     T.M.W. Nye July 94.
C     Performs correlation function Fourier transform on given data.
C     Calculates 1D and 3D correlation functions and writes them
C     to files X??000.CF1 and X??000.CF3.
C
C     Reference: Strobl und Schneider; J.Polymer Sci., polymer Phys. Ed.;
C     Vol. 18, 1343-1359, 1980
C
C     Update 9/8/94
C     Tailfit channel limits vary with channel no.
C
C     Update 11/8/94
C     Interface distribution function added.
C
C     Update 2005
C     ASCII output added.
C     Gamma3All added.  Fdi added.  Notokidf added.
C
C     Update Nov 2005
C     Redimensioned data arrays from 512 to MaxDim.
C     Added ASCII output of multi-frame data & array asciiresults

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

      EXTERNAL corfunc_io, wrdlen, getpath, strpath, endian, writeascii

      CHARACTER*80 dirname,filename,qaxname,othername
      CHARACTER*80 momname,momax,ascname
      CHARACTER*40 title,ascii,sigmodel,user,idftoggle
      CHARACTER*80 :: prompts(25)
      CHARACTER*80 fname,axisname,fname2,arbabs,graphics
      CHARACTER*40 backex,retrans
      CHARACTER*80 fully,fullx,transhead,transname
      CHARACTER*80 cor1name,cor3name,cor1head,cor3head,dhead,dname
      CHARACTER*80 idfhead,idfname
      CHARACTER*80 storename,outname,stats
      CHARACTER*1 letter
      INTEGER qzero,plotflag,retransflag
C      INTEGER realtime(4),realflag
      INTEGER :: realtime(4)
      INTEGER realflag
      INTEGER datastart,counts,upperlim,lword
C      REAL moment(MaxDim,5),idf(MaxDim,MaxDim),fdi(MaxDim,MaxDim)
      REAL*4 :: moment(MaxDim,5)
      REAL*4 :: idf(MaxDim,MaxDim)
      REAL*4 :: fdi(MaxDim,MaxDim)
C      DIMENSION notok(10),notokcf1(10),xdata(MaxDim),ydata(MaxDim)
      INTEGER :: notok(10)
      INTEGER :: notokcf1(10)
      REAL*4 :: xdata(MaxDim)
      REAL*4 :: ydata(MaxDim)
C      DIMENSION param(MaxDim,5)
C      DIMENSION trans2(MaxDim)
C      DIMENSION trans3(MaxDim)
C      DIMENSION notokidf(10)
      REAL*4 :: param(MaxDim,5)
      REAL*4 :: trans2(MaxDim)
      REAL*4 :: trans3(MaxDim)
      INTEGER :: notokidf(10)
C      DIMENSION gamma1(MaxDim),gamma1all(MaxDim*MaxDim),gamma3(MaxDim),
C     &          gamma3all(MaxDim*MaxDim),xgamma(MaxDim)
      REAL*4 :: gamma1(MaxDim)
      REAL*4 :: gamma1all(MaxDim*MaxDim)
      REAL*4 :: gamma3(MaxDim)
      REAL*4 :: gamma3all(MaxDim*MaxDim)
      REAL*4 :: xgamma(MaxDim)
C      DIMENSION asciiresults(MaxDim)
      REAL*4 :: asciiresults(MaxDim)

C     Array "param"
C     Parameter 1: Bonart background
C     Parameter 2: K (interface per unit volume)
C     Parameter 3: sigma (diffuse boundary thickness
C     Parameter 4: A or H1
C     Parameter 5: B or H2

C     There are lots of arrays:
C     however, realtime data is not all stored in memory
C     at the same time. Instead, each frame is read in and
C     handled individually, eg. gamma1 instead of being
C     gamma1(nframe,channel) is simply gamma1(channel).


1000  FORMAT(A1)
1010  FORMAT(A80)
1020  FORMAT(2x,I3)
1030  FORMAT(10I8)
1040  FORMAT(A40)
1050  FORMAT(2x,I3,2x,I3,2x,I3,2x,I3)
1060  FORMAT(E12.6,2x,E12.6,2x,E12.6,2x,E12.6,2x,E12.6)
1070  FORMAT(A10)
1080  FORMAT(/,1x,'Working...')
1090  FORMAT(/,1x,'100: Transforming frame ',I3,'...')
1100  FORMAT(1x,'100: Re-transforming...')
1110  FORMAT(1x,'Correlation function analysis ',
     &'output file from program transform.f.')
1120  FORMAT(1x,'File: ',A40,' Frame: ',I3,'.')
1130  FORMAT(1x,'R [Angstroms]    Gamma1        Gamma3')
1140  FORMAT(1x,'----------------------------------------')
1150  FORMAT(1x,E12.6,2x,E12.6,2x,E12.6)
1160  FORMAT(2x,I3,2x,I3,2x,I3)
1170  FORMAT(1x,'100: Calculating interface distribution function...')
1180  FORMAT(1x,"100: Written Gamma1: ",A10)
1185  FORMAT(1x,"100: Written ASCII Gamma1: ",A10)
1190  FORMAT(1x,"100: Written Gamma3: ",A10)
1195  FORMAT(1x,"100: Written ASCII Gamma3: ",A10)
1200  FORMAT(1x,"100: Written x-axis: ",A10)
1205  FORMAT(1x,"100: Written ASCII x-axis: ",A10)
1210  FORMAT(1x,"100: Written re-transform: ",A10)
1215  FORMAT(1x,"100: Written ASCII re-transform: ",A10)
1220  FORMAT(1x,"100: Written second moment v frame: ",A10)
1225  FORMAT(1x,"100: Written ASCII second moment v frame: ",A10)
1230  FORMAT(1x,"100: Written interface distribution func: ",A10)
1235  FORMAT(1x,"100: Written ASCII interface distribution func: ",A10)
1240  FORMAT(E12.6)

      prompts(1)='ERROR: Error reading corfunc.txt file: FATAL'
      prompts(2)='ERROR: Error reading data file: FATAL'
      prompts(3)='ERROR: Expecting static Q axis,'
     &//' received dynamic: FATAL'
      prompts(4)='All necessary files correctly loaded...'
      prompts(5)='ERROR: Extrapolated Q data does not '
     &//'reach Q=0.6: FATAL'
      prompts(6)='ERROR: Software error with subroutine'
     &//' changeotok: FATAL'
      prompts(7)='ERROR: Error creating otoko output files: FATAL'
      prompts(8)='Transforms completed...'
      prompts(9)='Sorry: no graphics for realtime data...'
      prompts(10)='DISPLAYING GRAPHICS: PRESS MIDDLE MOUSE BUTTON '
     &//'TO CONTINUE...'
      prompts(11)='ERROR: Error reading cor. func. otoko files '
     &//'while preparing ascii output: FATAL'
      prompts(12)='WARNING: Error writing ascii output'
      prompts(13)='Preparing ascii output...'
      prompts(14)='Creating moments output file...'
      prompts(15)='ERROR: Error writing extract.txt file: FATAL'
      prompts(16)='ERROR: Error with status file: FATAL'
      prompts(17)='Finished transform'
      prompts(18)='ERROR: Exceeded maximum number of steps for F.T. '
     &//'= 511: FATAL'
      prompts(19)='Enter maximum D for Fourier transform [Angstroms]'
      prompts(20)='Enter step in D for Fourier transform [Angstroms]'
      prompts(21)='Re-transform correlation function [y/n]'
      prompts(22)='Enter estimate of volume fraction crystallinity'
      prompts(23)='Calculate the interface distribution function [y/n]'
      prompts(24)='100: TRANSFORMING...'
      prompts(25)='100:'

      dmax=200.
      dstep=1.
      retrans='off'
      idftoggle='off'
      CALL WRDLEN(lword)

C      WRITE(6,*)
C      title='Correlation Function Calculation'
C      CALL showtitle(title)
      CALL showprompts(prompts(25))
      CALL showprompts(prompts(24))

C     Read corfunc.txt
      open(unit=9,file='corfunc.txt',STATUS='old',ERR=5000)
      read(9,1010,ERR=5000)filename
      read(9,1010,ERR=5000)qaxname
      read(9,1050,ERR=5000)realtime(1),realtime(2),realtime(3),
     &                     realtime(4)
      read(9,1010,ERR=5000)arbabs
      read(9,1040,ERR=5000)sigmodel
      read(9,1040,ERR=5000)backex
      read(9,1020,ERR=5000)qzero
      DO i=realtime(1),realtime(2),realtime(3)
          read(9,1060,ERR=5000)param(i,1),param(i,2),param(i,3),
     &                         param(i,4),param(i,5)
      END DO
      close(9)

C     Input Fourier transform data
C     Maximum D
2070  CALL defaultreal(prompts(19),dmax)
C      IF (dmax .LT. 100. .OR. dmax .GT. 1000.) THEN
C        CALL showprompts(prompts(18))
C        GOTO 2070
C      ENDIF
C     Step in D
2080  CALL defaultreal(prompts(20),dstep)
C      IF (dstep .LT. 0.25 .OR. dstep .GT. 4.) THEN
C        CALL showprompts(prompts(18))
C        GOTO 2080
C      ENDIF
      IF ((dmax/dstep) .GT. 511.) THEN
C     too many steps
        CALL showprompts(prompts(18))
        GOTO 2070
      ENDIF
C     Re-transform data?
      letter='n'
      CALL defaultletter(prompts(21),letter)
      IF (letter .EQ. 'y' .OR. letter .EQ. 'Y') THEN
        retrans='on'
      ELSE
        retrans='off'
      ENDIF
C     interface distribution function?
      letter='n'
      CALL defaultletter(prompts(23),letter)
      IF (letter .EQ. 'y' .OR. letter .EQ. 'Y') THEN
        idftoggle='on'
      ELSE
        idftoggle='off'
      ENDIF

C     Calculate number of points in correlation functions
      numd=INT(dmax/dstep)+1

      CALL getpath(filename,dirname)
      CALL strippath(filename,fname2)
      filename=fname2

C     Change filenames to extrapolated data filenames
      fully=filename
      CALL swapexten(fully,'FUL')
      fullx=qaxname
      CALL swapexten(fullx,'FLX')

C     Read intensity header
      fname=fully
      axisname=fullx
      OPEN(UNIT=9,FILE=fname,STATUS='old',ERR=5020)
      
      READ(9,1000,ERR=5020)letter
      READ(9,1000,ERR=5020)letter
      READ(9,1030,ERR=5020)notok
      READ(9,1040,ERR=5020)othername
      CLOSE(9)

      nndata=4*MaxDim

C     READ intensities
      fname2=othername
      OPEN(UNIT=9,FILE=fname2,STATUS='old',
     &ACCESS='direct',RECL=nndata/lword,ERR=5020)
      DO nframe=1,realtime(4)
C       NB: this is just to check the data exists.
        READ(9,REC=nframe,ERR=5020)(ydata(i),i=1,MaxDim)
      END DO
      CLOSE(9)

C     Save name of intensity file for later
      storename=fname2

C     Open Q axis header
      OPEN(UNIT=9,FILE=axisname,STATUS='old',ERR=5020)
      READ(9,1000,ERR=5020)letter
      READ(9,1000,ERR=5020)letter
      READ(9,1030,ERR=5020)notok
      READ(9,1040,ERR=5020)othername
      CLOSE(9)

C     Check static
      IF (notok(2) .NE. 1) THEN
        CALL showprompts(prompts(3))
        STOP
      ENDIF

C     Read x axis data
      fname2=othername
      OPEN(UNIT=9,FILE=fname2,STATUS='old',
     &ACCESS='direct',RECL=nndata/lword,ERR=5020)
      READ(9,REC=1,ERR=5020)(xdata(i),i=1,MaxDim)
2010  CLOSE(9)

C     Ok - data fine
      CALL showprompts(prompts(4))

      WRITE(6,1080)

C     Limit on all calculations: integrate out this far.
      qmax=0.6
C     Qmax=0.6 correspond to fluctuations in cor. func. at D ~ 10 Angst.

C     Calculate x axis for transforms: ie. real space coordinate.
      DO i=1,MaxDim
        xgamma(i)=dstep*(i-1)
      END DO

C     get outputname for 1D and 3D cor. func.s
      cor1head=filename
      cor3head=filename
      CALL swapexten(cor1head,'CF1')
      CALL swapexten(cor3head,'CF3')
      CALL changeotok(cor1head,cor1name,nerr)
      IF (nerr .EQ. 1) THEN
        CALL showprompts(prompts(6))
        STOP
      ENDIF
      CALL changeotok(cor3head,cor3name,nerr)
      IF (nerr .EQ. 1) THEN
        CALL showprompts(prompts(6))
        STOP
      ENDIF

C     Also prepare real~space axis filenames
      dhead=cor1head
      CALL swapexten(dhead,'RXA')
      CALL changeotok(dhead,dname,nerr)
      IF (nerr .EQ. 1) THEN
        CALL showprompts(prompts(6))
        STOP
      ENDIF

C     Write otoko file headers.
C     First do the real space axis.
      fname=dhead
      OPEN(UNIT=9,FILE=fname,STATUS='unknown',ERR=5030)
      WRITE(9,*,ERR=5030)
      WRITE(9,*,ERR=5030)
      notok(1)=numd
      notok(2)=1
      notok(3)=1
      CALL endian(notok(4))
      DO i=5,10
        notok(i)=0
      END DO
      WRITE(9,1030,ERR=5030)notok
      WRITE(9,1070,ERR=5030)dname(1:10)
      CLOSE(9)
C     X header done

C     Now do x data
      nrecl=4*MaxDim
      fname=dname
      OPEN(UNIT=10,FILE=fname,STATUS='unknown',
     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)
      WRITE(10,REC=1,ERR=5030)
     &(xgamma(i),i=1,numd)
      CLOSE(10)

      write(6,1200)dhead

C     Output the real axis in ASCII format
      CALL writeascii(fname(1:3)//"RXA."//"TXT",notok,xgamma,irc,1)
      IF(irc.EQ.0)THEN
        CALL showprompts(prompts(12))
      ELSE
        WRITE(6,1205)fname(1:3)//"RXA."//"TXT"
      ENDIF

C     Now create headers for 1D and 3D cor. func.s.
C     1 dim.
      fname=cor1head
      OPEN(UNIT=9,FILE=fname,STATUS='unknown',ERR=5030)
      WRITE(9,*,ERR=5030)
      WRITE(9,*,ERR=5030)
      notok(1)=numd
      notok(2)=realtime(4)
      notok(3)=1
      CALL endian(notok(4))
      DO i=5,10
        notok(i)=0
      END DO
      WRITE(9,1030,ERR=5030)notok
      WRITE(9,1070,ERR=5030)cor1name(1:10)
      CLOSE(9)

C     For the correlation functions
      notokcf1(1)=notok(1)
      notokcf1(2)=notok(2)
      DO i=3,10
        notokcf1(i)=0
      END DO

C     For the interface distribution function
C     (Also see the comment below where the IDF is written out)
C     notokidf(1)=MaxDim
      notokidf(1)=notok(1)
      notokidf(2)=notok(2)
      DO i=3,10
        notokidf(i)=0
      END DO

C     1D header done
C     3 dim.
      fname=cor3head
      OPEN(UNIT=9,FILE=fname,STATUS='unknown',ERR=5030)
      WRITE(9,*,ERR=5030)
      WRITE(9,*,ERR=5030)
      notok(1)=numd
      notok(2)=realtime(4)
      notok(3)=1
      CALL endian(notok(4))
      DO i=5,10
        notok(i)=0
      END DO
      WRITE(9,1030,ERR=5030)notok
      WRITE(9,1070,ERR=5030)cor3name(1:10)
      CLOSE(9)
C     3D header done

C     Also open data files
      nrecl=numd*4
      fname=cor1name
      OPEN(UNIT=9,FILE=fname,STATUS='unknown',
     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)
      fname=cor3name
      OPEN(UNIT=10,FILE=fname,STATUS='unknown',
     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)

C     Finally create header and open file for re-transformed data.
      IF (retrans .EQ. 'on' .OR. retrans .EQ. 'ON') THEN

        npts=MaxDim-qzero+1

C       get filenames
        transhead=filename
        CALL swapexten(transhead,'SMO')
        CALL changeotok(transhead,transname,nerr)
        IF (nerr .EQ. 1) THEN
          CALL showprompts(prompts(6))
          STOP
        ENDIF
C       Create header
        fname=transhead
        fname2=transhead(1:9)//'2'
        OPEN(UNIT=8,FILE=fname,STATUS='unknown',ERR=5030)
        OPEN(UNIT=11,FILE=fname2,STATUS='unknown',ERR=5030)
        WRITE(8,*,ERR=5030)
        WRITE(11,*,ERR=5030)
        WRITE(8,*,ERR=5030)
        WRITE(11,*,ERR=5030)
        notok(1)=npts+qzero-1
        notok(2)=realtime(4)
        notok(3)=1
        CALL endian(notok(4))
        DO i=5,10
          notok(i)=0
        END DO
        WRITE(8,1030,ERR=5030)notok
        WRITE(11,1030,ERR=5030)notok
        WRITE(8,1070,ERR=5030)transname(1:10)
        WRITE(11,1070,ERR=5030)transname(1:9)//'2'
        CLOSE(8)
        CLOSE(11)

C       Open file
        nrecl=(npts+qzero-1)*4
        fname=transname
        fname2=transname(1:9)//'2'
        OPEN(UNIT=8,FILE=fname,STATUS='unknown',
     &  ACCESS='direct',RECL=nrecl/lword,ERR=5030)
        OPEN(UNIT=11,FILE=fname2,STATUS='unknown',
     &  ACCESS='direct',RECL=nrecl/lword,ERR=5030)
C       All done
      ENDIF

C     START NUMBER CRUNCH
C     ~~~~~~~~~~~~~~~~~~~

C     LOOP THROUGH FRAMES
      nfr=0
      DO nframe=realtime(1),realtime(2),realtime(3)
        nfr=nfr+1

        WRITE(6,1090)nframe

C       Get data.
        OPEN(UNIT=7,FILE=storename,STATUS='old',
     &  ACCESS='direct',RECL=8192/lword,ERR=5020)
          READ(7,REC=nfr,ERR=5020)(ydata(i),i=1,MaxDim)
        CLOSE(7)

C       Calculate moments.
C       ~~~~~~~~~~~~~~~~~~
        DO j=1,5
          moment(nframe,j)=0.
        END DO
C       xdata(1) is not necessarily zero: take care of first point.
        x1=xdata(1)
        y1=ydata(1)
        DO j=1,5
          moment(nframe,j)=moment(nframe,j)+
     &    x1*0.5*y1*(x1**(j-1))
        END DO
        counts=1

C       Then int. wrt. Q, summing at each Q value.
2020    counts=counts+1
        x2=xdata(counts)
        y2=ydata(counts)
        DO j=1,5
          moment(nframe,j)=moment(nframe,j)+
     &    (x2-x1)*0.5*(y1*(x1**(j-1))+y2*(x2**(j-1)))
        END DO
        x1=x2
        y1=y2

C       Check to see whether we've gone out far enough in Q
        IF (x1 .LT. qmax .AND. counts .LT. MaxDim) THEN
          GOTO 2020
        ELSEIF (x1 .LT. qmax .AND. counts .EQ. MaxDim) THEN
          CALL showprompts(prompts(5))
          STOP
        ENDIF

C       End calculating moments


C       Calculate transforms.
C       ~~~~~~~~~~~~~~~~~~~~~
C       Reset cor. func.s
        DO i=1,MaxDim
          gamma1(i)=0.
          gamma3(i)=0.
        END DO

C       Value at the ordinate is defined as 1.0
        gamma1(1)=1.
        gamma3(1)=1.

C       Check whether we are going to retransform data:
C       calculate correlation function as far out as this requires.
        IF (retrans .EQ. 'on' .OR. retrans .EQ. 'ON') THEN
          nend=MaxDim
        ELSE
          nend=numd
        ENDIF

C       Loop through the values of d
        DO i=2,nend
          d=xgamma(i)

C         Evaluate integral for this value of d
          x1=xdata(1)
          y1=ydata(1)
          gamma1(i)=0.5*x1*(x1*x1*y1*COS(d*x1))
          gamma3(i)=0.5*x1*(x1*y1*SIN(d*x1)/d)
          counts=1

2030      counts=counts+1
          x2=xdata(counts)
          y2=ydata(counts)

C         One dimensional correlation function:
          temp1=(x1*x1*y1*COS(d*x1))
          temp2=(x2*x2*y2*COS(d*x2))
          gamma1(i)=gamma1(i)+0.5*(x2-x1)*(temp1+temp2)

C         Three dimensional cor. func.:
          temp1=(x1*y1*SIN(d*x1))/d
          temp2=(x2*y2*SIN(d*x2))/d
          gamma3(i)=gamma3(i)+0.5*(x2-x1)*(temp1+temp2)

C         Update values:
          x1=x2
          y1=y2

C         Check to see whether we've looked sufficiently far out in Q.
          IF (x1 .LT. qmax .AND. counts .LT. MaxDim) THEN
            GOTO 2030
          ELSEIF (x1 .LT. qmax .AND. counts .EQ. MaxDim) THEN
            CALL showprompts(prompts(5))
            STOP
          ENDIF

C         Calculation for one value of D is complete.
C         Normalise to the second moment:
          gamma1(i)=gamma1(i)/moment(nframe,3)
          gamma3(i)=gamma3(i)/moment(nframe,3)

C       End loop through D
        END DO

C       Output cor. func.s into otoko files.
C       1 dim.
        WRITE(9,REC=nfr,ERR=5030)(gamma1(i),i=1,numd)
C       3 dim.
        WRITE(10,REC=nfr,ERR=5030)(gamma3(i),i=1,numd)

        DO i=1,numd
          gamma1all(((nfr-1)*notokcf1(1))+i)=gamma1(i)
          gamma3all(((nfr-1)*notokcf1(1))+i)=gamma3(i)
        END DO

C       Interface distribution function.
C       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        IF (idftoggle .EQ. 'on' .OR. idftoggle .EQ. 'ON') THEN
          WRITE(6,1170)

C         Reset array
          DO i=1,MaxDim
            idf(nframe,i)=0.
          END DO

C         Value at the ordinate is given by moments
          idf(nframe,1)=-moment(nframe,5)/moment(nframe,3)

C         Number of points
          nend=numd

C         Loop through the values of d
          DO i=2,nend
            d=xgamma(i)

C           Evaluate integral for this value of d
            x1=xdata(1)
            y1=ydata(1)
            idf(nframe,i)=0.5*x1*(-(x1**4)*y1*COS(d*x1))
            counts=1

2050        counts=counts+1
            x2=xdata(counts)
            y2=ydata(counts)

C           One dimensional correlation function:
            temp1=(-(x1**4)*y1*COS(d*x1))
            temp2=(-(x2**4)*y2*COS(d*x2))
            idf(nframe,i)=idf(nframe,i)+0.5*(x2-x1)*(temp1+temp2)

C           Update values:
            x1=x2
            y1=y2

C           Check to see whether we've looked sufficiently far out in Q.
            IF (x1 .LT. qmax .AND. counts .LT. MaxDim) THEN
              GOTO 2050
            ELSEIF (x1 .LT. qmax .AND. counts .EQ. MaxDim) THEN
              CALL showprompts(prompts(5))
              STOP
            ENDIF

C           Calculation for one value of D is complete.
C           Normalise to the second moment:
            idf(nframe,i)=idf(nframe,i)/moment(nframe,3)

C         End loop through D
          END DO

C       Calc. idf?
        ENDIF


C       Re-transform if necessary
C       ~~~~~~~~~~~~

        IF (retrans .EQ. 'on' .OR. retrans .EQ. 'ON') THEN

          WRITE(6,1100)

C         Reset re-transform data
          DO i=1,MaxDim
            trans2(i)=0.
          END DO

C         Calculate number of points required

C         Loop through the values of Q
C         Ignore the first point - artefacts?
          DO i=2,npts
            q=xdata(i)

C           Evaluate integral for this value of Q
            x1=xgamma(1)
            y1=gamma1(1)
            trans2(i)=0.5*x1*(y1*COS(q*x1))
            counts=1

2040        counts=counts+1
            x2=xgamma(counts)
            y2=gamma1(counts)

            temp1=(y1*COS(q*x1))
            temp2=(y2*COS(q*x2))
            trans2(i)=trans2(i)+0.5*(x2-x1)*(temp1+temp2)

C           Update values:
            x1=x2
            y1=y2

C           Check to see whether we've looked sufficiently far out in D.
C           Change these conditions later.
            IF (counts .LT. MaxDim) THEN
              GOTO 2040
            ENDIF

C           Calculation for one value of Q is complete.
C           Normalise to the second moment:

            trans2(i)=trans2(i)*moment(nframe,3)/(3.1416*0.5)

C           I think that there is also a factor of pi/2 in there.
C           We transformed IQQ originally, so we must now divide
C           by Q**2 to get back to the original data.
            trans2(i)=trans2(i)/(q*q)

C         End loop through Q
          END DO

C         Shift transformed data according to qzero, so that the
C         original x-ray data x axis can be use for displaying
C         the retransformed data.
          DO i=MaxDim,qzero,-1
            trans2(i)=trans2(i-qzero+1)
            trans3(i)=trans2(i)+param(nframe,1)
          END DO
          DO i=qzero-1,1,-1
            trans2(i)=0.
            trans3(i)=0.
          END DO

C         Output retrans. data
          WRITE(8,REC=nfr,ERR=5030)(trans2(i),i=1,npts+qzero-1)
          WRITE(11,REC=nfr,ERR=5030)(trans3(i),i=1,npts+qzero-1)

C       Output the retrans data in ASCII format
        CALL writeascii(cor1head(1:3)//"SMO."//"TXT",notokcf1,trans2,
     &                  irc,1)

        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(12))
        ELSE
          WRITE(6,1215)cor1head(1:3)//"SMO."//"TXT"
        ENDIF
        CALL writeascii(cor1head(1:3)//"SM2."//"TXT",notokcf1,trans3,
     &                  irc,1)

        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(12))
        ELSE
          WRITE(6,1215)cor1head(1:3)//"SM2."//"TXT"
        ENDIF

C       End check whether retrans needed
        ENDIF


C     End Loop through frames
      END DO
      write(6,1180)cor1head
      write(6,1190)cor3head
      IF(retrans.eq."on".or.retrans.eq."ON")then
        write(6,1210)transhead(1:10)
        write(6,1210)transhead(1:9)//'2'
        CLOSE(11)
      ENDIF
      CLOSE(8)
      CLOSE(9)
      CLOSE(10)


C       Output the 1-D Correlation function in ASCII format
        CALL writeascii(cor1head(1:3)//"CF1."//"TXT",notokcf1,gamma1all,
     &                  irc,1)

        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(12))
        ELSE
          WRITE(6,1185)cor1head(1:3)//"CF1."//"TXT"
        ENDIF

C       Output the 3-D Correlation function in ASCII format
        CALL writeascii(cor3head(1:3)//"CF3."//"TXT",notokcf1,gamma3all,
     &                  irc,1)

        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(12))
        ELSE
          WRITE(6,1195)cor3head(1:3)//"CF3."//"TXT"
        ENDIF

C     Output interface distribution function if necessary.
      IF (idftoggle .EQ. 'on' .OR. idftoggle .EQ. 'ON') THEN

C       get filenames
        idfhead=filename
        CALL swapexten(idfhead,'IDF')
        CALL changeotok(idfhead,idfname,nerr)
        IF (nerr .EQ. 1) THEN
          CALL showprompts(prompts(6))
          STOP
        ENDIF

C       Write otoko header
        fname=idfhead
        OPEN(UNIT=9,FILE=fname,STATUS='unknown',ERR=5030)
        WRITE(9,*,ERR=5030)
        WRITE(9,*,ERR=5030)
        notok(1)=numd
        notok(2)=realtime(4)
        notok(3)=1
        CALL endian(notok(4))
        DO i=5,10
          notok(i)=0
        END DO
        WRITE(9,1030,ERR=5030)notok
        WRITE(9,1070,ERR=5030)idfname(1:10)
        CLOSE(9)
C       Header done

C       Now do idf data
        nrecl=4*numd
        fname=idfname
        OPEN(UNIT=10,FILE=fname,STATUS='unknown',
     &  ACCESS='direct',RECL=nrecl/lword,ERR=5030)
        nfr=0
        DO nframe=realtime(1),realtime(2),realtime(3)
          nfr=nfr+1
          WRITE(10,REC=nfr,ERR=5030)
     &    (idf(nframe,i),i=1,numd)
        END DO
        CLOSE(10)

      WRITE(6,1230)idfhead
C     End idf output
      ENDIF

C     FORTRAN and C store their 2D arrays in a different order, so
C     need to resequqnce idf() in order to use WRITEASCII.C to output
C     the IDF in ASCII!
C     Note that as this stands, WRITEASCII.C will ouput trailing zeroes
C     from the maximum channel number plus 1, upto the array dimension.
C     If the array dimension is increased it would be smarter to introduce
C     a step equal to the maximum channel dimension each frame.
      do i=1,MaxDim
        do j=1,MaxDim
           fdi(j,i)=0.
        end do
      end do
C     Loop over channels
      do i=1,notokidf(1)
C     Loop over frames
        do j=1,notokidf(2)
           fdi(i,j)=idf(j,i)
        end do
      end do

C       Output the idf in ASCII format
        CALL writeascii(idfhead(1:3)//"IDF."//"TXT",notokidf,fdi,
     &                  irc,1)

        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(12))
        ELSE
          WRITE(6,1235)idfhead(1:3)//"IDF."//"TXT"
        ENDIF


C     END NUMBER CRUNCH
C     ~~~~~~~~~~~~~~~~~

      CALL showprompts(prompts(8))

C     Also consider re-transformed data
      IF (retrans .EQ. 'on' .OR. retrans .EQ. 'ON') THEN
        retransflag=1
      ELSE
        retransflag=0
      ENDIF

C     Check whether realtime or static output.
      IF (realtime(4) .NE. 1) THEN

C       Realtime.
C       Create an otoko file containing the second moment and plot.

C       Get filename.
        momname=filename
        momax=filename
        CALL swapexten(momname,'MO2')
        CALL swapexten(momax,'FAX')
C       x-axis should already exist

C       Create otoko header
C       Change header name to give data filename
        CALL changeotok(momname,othername,nerr)
        IF (nerr .EQ. 1) THEN
          CALL showprompts(prompts(6))
          STOP
        ENDIF
        fname=momname
        OPEN(UNIT=9,FILE=fname,STATUS='unknown',ERR=5030)
        WRITE(9,*,ERR=5030)
        WRITE(9,*,ERR=5030)
        notok(1)=realtime(4)
        notok(2)=1
        notok(3)=1
        CALL endian(notok(4))
        DO i=5,10
          notok(i)=0
        END DO
        WRITE(9,1030,ERR=5030)notok
        WRITE(9,1070,ERR=5030)othername(1:10)
        CLOSE(9)
C       Header done

        nrecl=4*realtime(4)
        outname=othername
        OPEN(UNIT=10,FILE=outname,STATUS='unknown',
     &  ACCESS='direct',RECL=nrecl/lword,ERR=5030)
        WRITE(10,REC=1,ERR=5030)
     &  (moment(i,3),i=realtime(1),realtime(2),realtime(3))
        CLOSE(10)

        write(6,1220)fname

C       Output the second moment data in ASCII format
       do i=realtime(1),realtime(2),realtime(3)       
           asciiresults(i)=moment(i,3)
        end do
        CALL writeascii(fname(1:3)//"MO2."//"TXT",notok,asciiresults,
     &                  irc,1)
        IF(irc.EQ.0)THEN
          CALL showprompts(prompts(12))
        ELSE
          WRITE(6,1225)fname(1:3)//"MO2."//"TXT"
        ENDIF

C     End realtime or static
      ENDIF

C     Write extract.txt
      open(unit=9,file='extract.txt',STATUS='unknown',ERR=5010)
      write(9,1240,ERR=5010)dmax
      write(9,1240,ERR=5010)dstep
      DO i=realtime(1),realtime(2),realtime(3)
          write(9,1060,ERR=5010)moment(i,1),moment(i,2),moment(i,3),
     &                          moment(i,4),moment(i,5)
      END DO
      close(9)

C     End Program.
      CALL showprompts(prompts(17))
      STOP


C     Error Catches:

C     Error reading corfunc.txt
5000  CALL showprompts(prompts(1))
      STOP

C     Error writing extract.txt
5010  CALL showprompts(prompts(15))
      STOP

C     Error reading in intensities or x data
5020  CALL showprompts(prompts(2))
      STOP

C     Error writing otoko output.
5030  CALL showprompts(prompts(7))
      STOP

C     Error reading in correlation function otoko files.
5040  CALL showprompts(prompts(11))
      STOP

      END



      SUBROUTINE convertname(fname,num)
C     Alters a file name X??000.* to X??nnn.*
C     where nnn is a frame number specified by "num".
      CHARACTER*80 fname
      CHARACTER*1 :: letter(3)
      CHARACTER*4 chunk
      INTEGER :: digit(3)
      INTEGER num

1000  FORMAT(/,1x,'ERROR: Software error in subroutine convertname: ',
     &'problem with filename: FATAL')

C     Get num as a string
      n=num
      DO i=0,9
        IF (n .LT. 100) THEN
          digit(1)=i
          GOTO 10
        ELSE
          n=n-100
        ENDIF
      END DO

10    DO i=0,9
       IF (n .LT. 10) THEN
          digit(2)=i
          GOTO 20
        ELSE
          n=n-10
        ENDIF
      END DO

20    DO i=0,9
       IF (n .LT. 1) THEN
          digit(3)=i
          GOTO 30
        ELSE
          n=n-1
        ENDIF
      END DO

C     Convert single figure integers to strings
30    DO i=1,3
        IF (digit(i) .EQ. 0) letter(i)='0'
        IF (digit(i) .EQ. 1) letter(i)='1'
        IF (digit(i) .EQ. 2) letter(i)='2'
        IF (digit(i) .EQ. 3) letter(i)='3'
        IF (digit(i) .EQ. 4) letter(i)='4'
        IF (digit(i) .EQ. 5) letter(i)='5'
        IF (digit(i) .EQ. 6) letter(i)='6'
        IF (digit(i) .EQ. 7) letter(i)='7'
        IF (digit(i) .EQ. 8) letter(i)='8'
        IF (digit(i) .EQ. 9) letter(i)='9'
      END DO

C     Finally, change the filename.
      DO i=2,77
        chunk=fname(i:i+3)
        IF (chunk .EQ. '000.') THEN
          fname=fname(1:i-1)//letter(1)//letter(2)//
     &    letter(3)//'.ASC'
          GOTO 40
        ENDIF
      END DO

C     Error message:
      WRITE(6,1000)
      STOP

40    RETURN
      END

	  
C=============================================================================
C==============================================================================


      SUBROUTINE extract_par
C     TMWN August 94
C     Program analyses correlation function and extracts information
C     based on a lamellar model.
C     Also outputs Porod and moments results.
C
C     Reference: Strobl und Schneider; J.Polymer Sci., polymer Phys. Ed.;
C     Vol. 18, 1343-1359, 1980
C
C     Updates: many attempts to get the position of the linear section
C     correct.
C     It's quite easy to write an algorithm that looks good on paper
C     and that works well most of the time. However, the big problem
C     is that the linear section occurs at low R values (high Q)
C     and so is greatly affected by the tail fitting. Of course,
C     tailfitting is a slightly dubious procedure, and any noise in
C     the tail can throw things completely.
C     Algorithms deciding where to place the linear section involve
C     polynomial interpolation of the data points, and differetiation
C     to determine some sort of gradient at each point.
C
C     8/8/94
C     Re-write linear section fit without interpolation or
C     differentiation. Algorithm has following form:
C     Loop through all possible start points for section.
C       Loop through all possible end points.
C         Lsqu fit to linear section.
C         Fit polynomials to curved sections at either end.
C         Calculate chisqu.
C       End loop.
C     End loop.
C
C     9/8/94 Update.
C     Handles channel limits for tailfit varying with frame no.
C
C     SMK Nov 05 Update.
C     Redimensioned data arrays from 512 to MaxDim.
C
      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

      EXTERNAL corfunc_io, wrdlen, getpath, strpath

      CHARACTER*80 dirname,filename,qaxname,filename2,resultname
      CHARACTER*40 title,ascii,sigmodel,user
      CHARACTER*80 :: prompts(25)
      CHARACTER*80 fname,arbabs,graphics,stats
      CHARACTER*60 text
      CHARACTER*40 backex,retrans,idftoggle
      CHARACTER*80 cor1name,cor3name,cor1head,cor3head,dhead,dname
      CHARACTER*1 letter
      INTEGER qzero
C      INTEGER realtime(4),realflag
C      INTEGER datastart,channel(MaxDim,2),lword
      INTEGER :: realtime(4)
      INTEGER :: channel(MaxDim,2)
      INTEGER datastart,lword,realflag
C      REAL moment(MaxDim,5),lp,intercept,volfrac
      REAL*4 :: moment(MaxDim,5)
      REAL lp,intercept,volfrac
C      DIMENSION notok(10),param(MaxDim,5),calc(MaxDim)
      INTEGER :: notok(10)
      REAL*4 :: param(MaxDim,5)
      REAL*4 :: calc(MaxDim)
C      DIMENSION gamma1(2048),xgamma(2048)
C      DIMENSION gamma1(MaxDim),xgamma(MaxDim)
      REAL*4 :: gamma1(MaxDim)
      REAL*4 :: xgamma(MaxDim)
C      DIMENSION reslam(MaxDim,11),resporod(MaxDim,7),nout(2)
      REAL*4 :: reslam(MaxDim,11)
      REAL*4 :: resporod(MaxDim,7)
      INTEGER :: nout(2)
      LOGICAL static_bn

C     Array "param"
C     Parameter 1: Bonart background
C     Parameter 2: K (interface per unit volume)
C     Parameter 3: sigma (diffuse boundary thickness
C     Parameter 4: A or H1
C     Parameter 5: B or H2

1000  FORMAT(A1)
1010  FORMAT(A80)
1020  FORMAT(2x,I3)
1030  FORMAT(10I8)
1040  FORMAT(A40)
1050  FORMAT(2x,I3,2x,I3,2x,I3,2x,I3)
1060  FORMAT(E12.6,2x,E12.6,2x,E12.6,2x,E12.6,2x,E12.6)
1070  FORMAT(A10)
1080  FORMAT(/,1x,'100: Frame ',I3,'. Gammamin ',E12.6,'.')
1090  FORMAT(2x,I3,2x,I3)
1100  FORMAT(/,1x,'Limits for frame ',I3,':')

      prompts(1)='ERROR: Error reading extract.txt file: FATAL'
      prompts(2)='ERROR: Expecting static Q axis,'
     &//' received dynamic: FATAL'
      prompts(3)='All necessary files correctly loaded...'
      prompts(4)='ERROR: Software error with subroutine'
     &//' changeotok: FATAL'
      prompts(5)='ERROR: Error reading corfunc.txt file: FATAL'
      prompts(6)='ERROR: Error reading correlation function '
     &//'otoko files: FATAL'
      prompts(7)='ERROR: Software error with subroutine'
     &//' changeotok: FATAL'
      prompts(8)='ERROR: Error with number of data points in '
     &//'correlation function otoko files: FATAL'
      prompts(9)='Please wait: working...'
      prompts(10)='ERROR: Slope on correlation function '
     &//'positive: FATAL'
      prompts(11)='ERROR: Error writing results file: FATAL'
      prompts(12)='Results coming up...'
      prompts(13)='Writing results to disk...'
      prompts(14)='Analysis session over: terminating as normal...'
      prompts(15)='ERROR: Error with status file: FATAL'
      prompts(16)='ERROR: Problem with limits: FATAL'
      prompts(17)='Enter R at start of linear section'
      prompts(18)='Enter R at end of linear section'
      prompts(19)='ERROR: Error with limitinfo.txt: FATAL'
      prompts(20)='Finished extract'
      prompts(21)='Do you want user control of extraction '
     &//'process [y/n]'
      prompts(22)='100: EXTRACTING STRUCTURAL PARAMETERS...'
      prompts(23)='100:'
      prompts(24)='Enter estimate of volume fraction crystallinity'
      prompts(25)='ERROR: Invalid input: FATAL'

      CALL WRDLEN(LWORD)
      realflag=1
      if(static_bn)then
        realflag=0
        realtime(1)=1
        realtime(2)=1
        realtime(3)=1
        realtime(4)=1
      endif

C      WRITE(6,*)
C      title='Extraction of information from corfunc'
C      CALL showtitle(title)
      CALL showprompts(prompts(23))
      CALL showprompts(prompts(22))

C     Read corfunc.txt
      open(unit=9,file='corfunc.txt',STATUS='old',ERR=5000)
      read(9,1010,ERR=5000)filename
      read(9,1010,ERR=5000)qaxname
      read(9,1050,ERR=5000)realtime(1),realtime(2),realtime(3),
     &                     realtime(4)
      read(9,1010,ERR=5000)arbabs
      read(9,1040,ERR=5000)sigmodel
      read(9,1040,ERR=5000)backex
      read(9,1020,ERR=5000)qzero
      DO i=realtime(1),realtime(2),realtime(3)
          read(9,1060,ERR=5000)param(i,1),param(i,2),param(i,3),
     &                         param(i,4),param(i,5)
      END DO
      DO i=realtime(1),realtime(2),realtime(3)
          read(9,1090,ERR=5000)channel(i,2),channel(i,1)
      END DO
      close(9)

C     Read extract.txt
      open(unit=9,file='extract.txt',STATUS='old',ERR=5010)
      read(9,*,ERR=5010)dmax
      read(9,*,ERR=5010)dstep
      DO i=realtime(1),realtime(2),realtime(3)
          read(9,1060,ERR=5010)moment(i,1),moment(i,2),moment(i,3),
     &                         moment(i,4),moment(i,5)
      END DO
      close(9)

      volfrac=0.5
C     Guess a volume fraction cryst.
2090  CALL defaultreal(prompts(24),volfrac)
      IF (volfrac .LT. 0. .OR. volfrac .GT. 1.) THEN
        CALL showprompts(prompts(15))
        GOTO 2090
      ENDIF

      letter='n'
      CALL defaultletter(prompts(21),letter)
      IF (letter .EQ. 'y' .OR. letter .EQ. 'Y') THEN
        user='user'
      ELSE
        user='auto'
      ENDIF

C     Calculate number of points in correlation functions
      numd=INT(dmax/dstep)+1

      CALL strippath(filename,fname)
      CALL getpath(filename,dirname)
      filename2=filename
      filename=fname

C     Read in D axis data.
C     Get filename.
      dhead=filename
      CALL swapexten(dhead,'RXA')
      fname=dhead
      OPEN(UNIT=10,FILE=fname,STATUS='old',ERR=5030)
      READ(10,1000,ERR=5030)letter
      READ(10,1000,ERR=5030)letter
      READ(10,1030,ERR=5030)notok
      READ(10,1040,ERR=5030)dname
      CLOSE(10)

C     Check static
      IF (notok(2) .NE. 1) THEN
        CALL showprompts(prompts(2))
        STOP
      ENDIF

C     Read D axis data
      nrecl=4*numd
      fname=dname
      OPEN(UNIT=9,FILE=fname,STATUS='old',
     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)
      READ(9,REC=1,ERR=5030)(xgamma(i),i=1,numd)
2010  CLOSE(9)

C     Open corfunc files: read headers.
C     Deal with filenames.
      cor1head=filename
      cor3head=filename
      CALL swapexten(cor1head,'CF1')
      CALL swapexten(cor3head,'CF3')
      CALL changeotok(cor1head,cor1name,nerr)
      IF (nerr .EQ. 1) THEN
        CALL showprompts(prompts(7))
        STOP
      ENDIF
      CALL changeotok(cor3head,cor3name,nerr)
      IF (nerr .EQ. 1) THEN
        CALL showprompts(prompts(7))
        STOP
      ENDIF

C     1 dim.
      fname=cor1head
      OPEN(UNIT=9,FILE=fname,STATUS='old',ERR=5030)
      READ(9,1000,ERR=5030)letter
      READ(9,1000,ERR=5030)letter
      READ(9,1030,ERR=5030)notok
      READ(9,1000,ERR=5030)letter
      CLOSE(9)
C     1D header done
C     Check c. realtime
      IF (notok(1) .NE. numd .OR. notok(2) .NE. realtime(4)) THEN
        CALL showprompts(prompts(8))
        STOP
      ENDIF

C     3 dim.
C      fname=cor3head
C      OPEN(UNIT=9,FILE=fname,STATUS='old',ERR=5030)
C      READ(9,1000,ERR=5030)letter
C      READ(9,1000,ERR=5030)letter
C      READ(9,1030,ERR=5030)notok
C      READ(9,1000,ERR=5030)letter
C      CLOSE(9)
C     Check c. realtime
C      IF (notok(1) .NE. numd .OR. notok(2) .NE. realtime(4)) THEN
C        CALL showprompts(prompts(8))
C        STOP
C      ENDIF
C     3D header done

C     Open data files:
      nrecl=4*numd
      fname=cor1name
      OPEN(UNIT=9,FILE=fname,STATUS='old',
     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)
C     At the moment gamma3 is not used in the analysis.
C      OPEN(UNIT=10,FILE=fname,STATUS='old',
C     &ACCESS='direct',RECL=nrecl/lword,ERR=5030)

C     Got all necessary data.
      CALL showprompts(prompts(3))


C     START PROCESSING
C     ~~~~~~~~~~~~~~~~

      CALL showprompts(prompts(9))

C     Loop through frames.
      nfr=0
      rstart=0.
      rend=0.
      DO nframe=realtime(1),realtime(2),realtime(3)
        nfr=nfr+1

C       Read in gamma1 data.
        READ(9,REC=nfr,ERR=5030)(gamma1(i),i=1,numd)

C       Search for the first minimum.
        DO i=2,numd-1
          IF ((gamma1(i) .LT. gamma1(i-1)) .AND.
     &    (gamma1(i) .LT. gamma1(i+1))) THEN
            IF (gamma1(i) .LT. 0.) THEN
              r3=xgamma(i)
              gammamin=-gamma1(i)
              minchannel=i
              lamflag=1
              GOTO 2020
            ENDIF
          ENDIF
        END DO

C       OK: we haven't found the first minimum
        lamflag=0
        minchannel=2

C       Search for local maxm.
2020    DO i=minchannel,(numd-1)
          IF ((gamma1(i) .GT. gamma1(i-1)) .AND.
     &     (gamma1(i) .GT. gamma1(i+1))) THEN
             IF (gamma1(i) .GT. 0.) THEN
               lp=xgamma(i)
               gammamax=gamma1(i)
               maxchannel=i
               lamflag=lamflag+1
               GOTO 2030
             ENDIF
           ENDIF
         END DO

C       OK: we haven't found the first maxm.
        lamflag=0
        maxchannel=2

C       Sort out flags
2030    IF (lamflag .NE. 2) THEN
          lamflag=0
        ELSE
          lamflag=1
        ENDIF

C       Check whether lammella information can be extracted.
        IF (lamflag .EQ. 1) THEN
C         Yes - we can go ahead.
C         Linear fit first.

C         8/8/94
C         New method for finding linear section.

C         Set up chisqu search
          oldchisqu=1.E+20
          linstart=1
          linend=minchannel
          grad=0.
          gammastar=0.

C         Loop through possible start points.
          DO i=2,minchannel-2

C           Loop through possible end points.
            DO j=i+1,minchannel-1

C             LSQ fit to linear section.
              sum1=0.
              sum2=0.
              sum3=0.
              sum4=0.
              ntemp=j-i+1
              DO k=i,j
                sum1=sum1+gamma1(k)
                sum2=sum2+xgamma(k)
                sum3=sum3+gamma1(k)*xgamma(k)
                sum4=sum4+xgamma(k)*xgamma(k)
              END DO
              gradient=(ntemp*sum3-sum1*sum2)/(ntemp*sum4-sum2*sum2)
              intercept=(sum1-gradient*sum2)/ntemp

C             Check whether polynom model applies:
              gr1=intercept+gradient*xgamma(i)
              gr2=intercept+gradient*xgamma(j)
              IF (intercept .GT. 1. .AND. gradient .LT. 0.
     &        .AND. gr1 .LT. 1. .AND. gr2 .GT. gammamin) THEN

C               Calculate polynomial parameters.
                r1=xgamma(i)
                r2=xgamma(j)
C               Low R end:
                p1=(-gradient*r1)/(1.-(gradient*r1)-intercept)
                a1=(1.-(gradient*r1)-intercept)/(r1**p1)
C               High R end:
                deltar=ABS(r2-xgamma(minchannel))
                p2=(-gradient*deltar)/(gammamin+intercept+(gradient*r2))
                a2=(gammamin+intercept+(gradient*r2))/(deltar**p2)

C               Fill array with calc. data
C               First low R curved section:
                DO k=1,i-1
                  x=xgamma(k)
                  calc(k)=1.-a1*(x**p1)
                END DO
C               Then linear section:
                DO k=i,j
                  x=xgamma(k)
                  calc(k)=intercept+gradient*x
                END DO
C               Then high R curved end:
                DO k=j+1,minchannel
                  x=ABS(xgamma(k)-xgamma(minchannel))
                  calc(k)=a2*(x**p2)-gammamin
                END DO

C               Calculate chisqu:
                chisqu=0.
                DO k=1,minchannel
                  chisqu=chisqu+(gamma1(k)-calc(k))**2
                END DO

C               Compare chisqu:
                IF (chisqu .LT. oldchisqu) THEN
                  oldchisqu=chisqu
                  linstart=i
                  linend=j
                  grad=gradient
                  gammastar=intercept
                ENDIF

C               Debug:
C                  write(6,*)gradient,intercept
C                 write(6,3000)i,j,a1,p1,a2,p2,chisqu
3000             FORMAT(I3,2x,I3,2x,E12.6,2x,E12.6,2x,
     &           E12.6,2x,E12.6,2x,E12.6)

              ENDIF

C           End loop through linear sections:
            END DO
          END DO

C         User control
          IF (user .EQ. 'user') THEN

C           Check whether first time
            IF (rstart .EQ. 0. .AND. rend .EQ. 0.) THEN
              rstart=xgamma(linstart)
              rend=xgamma(linend)
            ENDIF

C           prompts user
2040        WRITE(6,1080)nframe,xgamma(minchannel)
            WRITE(6,1100)nframe
            CALL defaultreal(prompts(17),rstart)
            CALL defaultreal(prompts(18),rend)

C           Turn R values into channel no.s
            linstart=0
            DO i=1,minchannel
              IF (xgamma(i) .LE. rstart .AND.
     &        xgamma(i+1) .GT. rstart) THEN
                linstart=i
              ENDIF
            END DO
            linend=0
            DO i=linstart,minchannel
              IF (xgamma(i) .LE. rend .AND.
     &        xgamma(i+1) .GT. rend) THEN
                linend=i
              ENDIF
            END DO

C           Check channel no.s
            IF (linstart .EQ. 0 .OR. linend .EQ. 0) THEN
              CALL showprompts(prompts(16))
              GOTO 2040
            ENDIF

C           Perform linear fit
C           LSQ fit to linear section.
            sum1=0.
            sum2=0.
            sum3=0.
            sum4=0.
            ntemp=linend-linstart+1
            DO i=linstart,linend
              sum1=sum1+gamma1(i)
              sum2=sum2+xgamma(i)
              sum3=sum3+gamma1(i)*xgamma(i)
              sum4=sum4+xgamma(i)*xgamma(i)
            END DO
            grad=(ntemp*sum3-sum1*sum2)/(ntemp*sum4-sum2*sum2)
            gammastar=(sum1-grad*sum2)/ntemp

C         End user control
          ENDIF

C          Debug:
C          IF (rstart .EQ. 0. .AND. rend .EQ. 0.) THEN
C            rstart=xgamma(linstart)
C            rend=xgamma(linend)
C          ENDIF
C          Write(6,7114)nframe,rstart,rend

C         Debug:
C         WRITE(6,*)grad,gammastar

C         That's everything! All measurements complete.
C         Now fill arrays with results.
C         Array "reslam" holds lamella results.

C         reslam(nframe,1)=Long period
C         (D position of first maxm).
C         reslam(nframe,2)=Average hard block thickness.
C         (D position of Xn. of extended linear section and base line).
C         reslam(nframe,3)=Average soft block thickness
C         (long period - hard block thickness).
C         reslam(nframe,4)=Bulk volume cryst. (ie. volume fraction)
C         (= 1. - (D where linear section crosses D axis)/(av. hard block)).
C         reslam(nframe,5)=Local cryst.
C         (= av. hard block / long period).
C         reslam(nframe,6)=Average core thickness
C         (= D at end of linear section).
C         reslam(nframe,7)=Interface thickness
C         (=D at start of linear section).
C         reslam(nframe,8)=Polydispersity
C         (= gamma at first min / gamma at first max)
C         reslam(nframe,9)=Electron density contrast
C         (= SQRT { Q / (bulk cryst)(1. - bulk cryst) }
C         reslam(nframe,10)=Specific inner surface
C         (= 2 * bulk cryst / av. hard block) ????
C         reslam(nframe,11)="Non-ideality" I made this one up myself!
C         (= (D at max - 2. * D at min)**2 / long period**2 ).
C         Non ideality has changed. And then changed back to the original.

          reslam(nframe,1)=xgamma(maxchannel)

          reslam(nframe,2)=(-gammamin-gammastar)/grad

          reslam(nframe,3)=reslam(nframe,1)-reslam(nframe,2)

          reslam(nframe,4)=1.+(gammastar/(grad*reslam(nframe,2)))

          reslam(nframe,5)=reslam(nframe,2)/reslam(nframe,1)

          reslam(nframe,6)=xgamma(linend)

          reslam(nframe,7)=xgamma(linstart)

          reslam(nframe,8)=gammamin/gammamax

C         Electron density contrast is the only result affected by
C         Absolute units of intensity.
          IF (arbabs .EQ. 'ABS' .OR. arbabs .EQ. 'abs') THEN
C           We're in abs units: convert invariant to cm^-4.
            qtemp=moment(nframe,3)*1.E+24
            reslam(nframe,9)=
     &      SQRT(qtemp*gammastar/((reslam(nframe,4))*(1.-reslam
     &           (nframe,4))))
          ELSE
C           We're in arb units
            reslam(nframe,9)=SQRT(moment(nframe,3)*gammastar/
     &      ((reslam(nframe,4))*(1.-reslam(nframe,4))))
          ENDIF

          reslam(nframe,10)=2*reslam(nframe,4)/reslam(nframe,2)

          reslam(nframe,11)=((reslam(nframe,1)-2.*xgamma(minchannel)
     &    )**2)/(reslam(nframe,1)**2)

C         End lamella calculations.

C         Now do Porod calculations.
C         NB!
C         Porod results use volume fraction hard block found from
C         the correlation function analysis.
          cryst=reslam(nframe,4)

C         But remember to check which model was used for tail fitting.
          IF (sigmodel .EQ. 'on' .OR. sigmodel .EQ. 'ON') THEN
C           Sigmoid model used for tail.

C           "resporod" holds results.
C           resporod(nframe,1)=Porod const.
C           resporod(nfarme,2)=sigma
C           resporod(nframe,3)=-1 (indicates sigmoid model)
C           resporod(nframe,4)=-1             "
C           resporod(nframe,5)=-1             "
C           resporod(nframe,6)=-1             "
C           resporod(nframe,7)=-1 (indicates lamella measurements successful)

            resporod(nframe,1)=param(nframe,2)

C           Take care of abs units (assume cm^-1)
            IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
C             Convert to cm^-5
              resporod(nframe,1)=resporod(nframe,1)*1.E+32
            ENDIF
            resporod(nframe,2)=param(nframe,3)
            resporod(nframe,3)=-1.
            resporod(nframe,4)=-1.
            resporod(nframe,5)=-1.
            resporod(nframe,6)=-1.
            resporod(nframe,7)=-1.

          ELSE
C           Ordinary Porod fit to tail.
C           "resporod" holds results.
C           resporod(nframe,1)=Porod const.
C           resporod(nfarme,2)=-1 (indicates Porod fit not sigmoid)
C           resporod(nframe,3)=Characteristic chord length
C           resporod(nframe,4)=Hard block chord length
C           resporod(nframe,5)=Soft block chord length
C           resporod(nframe,6)=surface to volume ratio
C           resporod(nframe,7)=-1 (indicates lamella measurements successful)

            resporod(nframe,1)=param(nframe,2)
C           Take care of abs units (assume cm^-1)
            IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
C             Convert to cm^-5
              resporod(nframe,1)=resporod(nframe,1)*1.E+32
            ENDIF
            resporod(nframe,2)=-1.
C           Chord length = 4Q / (pi * Kp)
            resporod(nframe,3)=4.*moment(nframe,3)/
     &      (param(nframe,2)*3.14)
            resporod(nframe,5)=resporod(nframe,3)/cryst
            resporod(nframe,4)=resporod(nframe,3)/(1.-cryst)
            resporod(nframe,6)=4.*cryst*(1.-cryst)/
     &                         resporod(nframe,3)
            resporod(nframe,7)=-1.

C         End check on which model for tail.
          ENDIF

        ELSE
C         Lamella information can't be extracted:
C         Just do porod results based on user's estimate of crystallinity.
          cryst=volfrac

C         First toggle lamella results
          DO i=1,11
            reslam(nframe,i)=-1.
          END DO

C         Remember to check which model was used for tail fitting.
          IF (sigmodel .EQ. 'on' .OR. sigmodel .EQ. 'ON') THEN
C           Sigmoid model used for tail.

C           "resporod" holds results.
C           resporod(nframe,1)=Porod const.
C           resporod(nfarme,2)=sigma
C           resporod(nframe,3)=-1 (indicates sigmoid model)
C           resporod(nframe,4)=-1             "
C           resporod(nframe,5)=-1             "
C           resporod(nframe,6)=-1             "
C           resporod(nframe,7)=Electron density contrast based on
C                              user supplied crystallinity.

            resporod(nframe,1)=param(nframe,2)
C           Take care of abs units (assume cm^-1)
            IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
C             Convert to cm^-5
              resporod(nframe,1)=resporod(nframe,1)*1.E+32
            ENDIF
            resporod(nframe,2)=param(nframe,3)
            resporod(nframe,3)=-1.
            resporod(nframe,4)=-1.
            resporod(nframe,5)=-1.
            resporod(nframe,6)=-1.

C           Electron density contrast is affected by abs. units
            IF (arbabs .EQ. 'ABS' .OR. arbabs .EQ. 'abs') THEN
C             We're in abs units: convert invariant to cm^-4.
              qtemp=moment(nframe,3)*1.E+24
              resporod(nframe,7)=SQRT(qtemp/
     &        ((cryst)*(1.-cryst)))
            ELSE
C             We're in arb units
              resporod(nframe,7)=SQRT(moment(nframe,3)/
     &        ((cryst)*(1.-cryst)))
            ENDIF

          ELSE
C           Ordinary Porod fit to tail.
C           "resporod" holds results.
C           resporod(nframe,1)=Porod const.
C           resporod(nfarme,2)=-1 (indicates Porod fit not sigmoid)
C           resporod(nframe,3)=Characteristic chord length
C           resporod(nframe,4)=Hard block chord length
C           resporod(nframe,5)=Soft block chord length
C           resporod(nframe,6)=surface to volume ratio
C           resporod(nframe,7)=Electron density contrast based on
C                              user supplied crystallinity.

            resporod(nframe,1)=param(nframe,2)
C           Take care of abs units (assume cm^-1)
            IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
C             Convert to cm^-5
              resporod(nframe,1)=resporod(nframe,1)*1.E+32
            ENDIF
            resporod(nframe,2)=-1.
C           Chord length = 4Q / (pi * Kp)
            resporod(nframe,3)=4.*moment(nframe,3)/
     &      (param(nframe,2)*3.14)
            resporod(nframe,5)=resporod(nframe,3)/cryst
            resporod(nframe,4)=resporod(nframe,3)/(1.-cryst)


            resporod(nframe,6)=4.*cryst*(1.-cryst)/
     &                         resporod(nframe,3)

C           Electron density contrast is affected by abs. units
            IF (arbabs .EQ. 'ABS' .OR. arbabs .EQ. 'abs') THEN
C             We're in abs units: convert invariant to cm^-4.
              qtemp=moment(nframe,3)*1.E+24
              resporod(nframe,7)=SQRT(qtemp/
     &        ((cryst)*(1.-cryst)))
            ELSE
C             We're in arb units
              resporod(nframe,7)=SQRT(moment(nframe,3)/
     &        ((cryst)*(1.-cryst)))
            ENDIF


C         End check on which model for tail.
          ENDIF


C       End lamella results or just porod
        ENDIF


C     End loop through frames.
      END DO

      CLOSE(9)
C      CLOSE(10)


C     END EXTRACTING RESULTS
C     ~~~~~~~~~~~~~~~~~~~~~~


C     START DISPLAYING RESULTS
C     ~~~~~~~~~~~~~~~~~~~~~~~~

C     NB differences between single frame and realtime.
C     Results:-
C     1. Title, filename, frames.
C     2. Moments.
C     3. Lamella results (if there are any).
C     4. Porod results (depending on tail model).
C
C     nout(1) = 6 is screen; nout(2) = 10 is file.
C
C     Remember to keep checking for abs units.

      nout(1)=6
      nout(2)=10

C     Output file: get filename.
      resultname=filename
      CALL swapexten(resultname,'LIS')
      fname=resultname
C     Open file.
      OPEN(UNIT=nout(2),FILE=fname,STATUS='unknown',ERR=5040)

C     Loop through screen and file.
      DO l=1,2
        nunit=nout(l)
        IF (l.EQ.1) THEN
C         screen.
          CALL showprompts(prompts(12))
          WRITE(6,*)
        ELSE
C         file.
          CALL showprompts(prompts(13))
        ENDIF


C       OK: start output with a title.
        WRITE(nunit,7000,ERR=5040)
7000    FORMAT(1x,'CORRELATION FUNCTION ANALYSIS OUTPUT',/,
     &  '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',/)
C       Next directory and filename.
        IF(l.EQ.1)THEN
            WRITE(nunit,7010,ERR=5040)filename2
        ELSE
            WRITE(nunit,7011,ERR=5040)filename2
        ENDIF
7010    FORMAT(1x,'100: Analysis performed on ',/,
     &  1x,'100: Data file: ',A40)
7011    FORMAT(1x,'Analysis performed on ',/,
     &  1x,'Data file: ',A40)
C       Realtime info.
        IF (realflag .EQ. 0) THEN
          IF(l.EQ.1) THEN
            WRITE(nunit,7020,ERR=5040)
          ELSE
            WRITE(nunit,7021,ERR=5040)
          ENDIF
7020      FORMAT(/,1x,'100: Realtime status: Static.')
7021      FORMAT(/,1x,'Realtime status: Static.')
        ELSE
          IF(l.EQ.1) THEN
            WRITE(nunit,7030,ERR=5040)realtime(1),realtime(2),
     &                                realtime(3),realtime(4)
          ELSE
            WRITE(nunit,7031,ERR=5040)realtime(1),realtime(2),
     &                                realtime(3),realtime(4)
          ENDIF
7030      FORMAT(1x,'100: Start frame: ',I3,'.       End frame: '
     &,I3,'.',/,1x,'100: Frame increment: ',I3,'.   Total number '
     &,'of frames: ',I3,'.')
7031      FORMAT(1x,'Start frame: ',I3,'.       End frame: '
     &,I3,'.',/,1x,'Frame increment: ',I3,'.   Total number '
     &,'of frames: ',I3,'.')
        ENDIF
C       Arb or abs units?
        IF(l.EQ.1)THEN
          IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
            WRITE(nunit,7040,ERR=5040)'Absolute (cm^-1).'
          ELSE
            WRITE(nunit,7040,ERR=5040)'Arbitrary.       '
          ENDIF
        ELSE
          IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
            WRITE(nunit,7041,ERR=5040)'Absolute (cm^-1).'
          ELSE
            WRITE(nunit,7041,ERR=5040)'Arbitrary.       '
          ENDIF
        ENDIF
7040    FORMAT(1x,'100: Intensity units: ',A17)
7041    FORMAT(1x,'Intensity units: ',A17)
C       Q axis
        IF(l.EQ.1)THEN
          WRITE(nunit,7050,ERR=5040)qaxname
        ELSE
          WRITE(nunit,7051,ERR=5040)qaxname
        ENDIF
7050    FORMAT(1x,'100: Q axis filename: ',A40)
7051    FORMAT(1x,'Q axis filename: ',A40)
C       Information on the tailfit.
        IF(l.EQ.1)THEN
          WRITE(nunit,7060,ERR=5040)
        ELSE
          WRITE(nunit,7061,ERR=5040)
        ENDIF
7060    FORMAT(//,1x,'100: EXTRAPOLATION PARAMETERS:')
7061    FORMAT(//,1x,'EXTRAPOLATION PARAMETERS:')
        IF(l.EQ.1)THEN
          IF (sigmodel .EQ. 'on' .OR. sigmodel .EQ. 'ON') THEN
            WRITE(nunit,7070,ERR=5040)'Sigmoid.'
          ELSE
            WRITE(nunit,7070,ERR=5040)'Porod.  '
          ENDIF
        ELSE
          IF (sigmodel .EQ. 'on' .OR. sigmodel .EQ. 'ON') THEN
            WRITE(nunit,7071,ERR=5040)'Sigmoid.'
          ELSE
            WRITE(nunit,7071,ERR=5040)'Porod.  '
          ENDIF
        ENDIF
7070    FORMAT(1x,'100: Tail model applied: ',A8)
7071    FORMAT(1x,'Tail model applied: ',A8)
C       Update 9/8/94
C       Check realtime
        IF(l.EQ.1)THEN
          IF (realtime(4) .EQ. 1) THEN
            WRITE(nunit,7080,ERR=5040)
     &      channel(realtime(1),1),channel(realtime(1),2)
7080        FORMAT(1x,'100: Tailfit start channel ',I3,
     &      '.',/,1x,'100: End channel ',I3,'.')
          ELSE
C           Realtime data
            WRITE(nunit,7081,ERR=5040)
7081        FORMAT(/,1x,'100: Frame    Start Channel    End Channel',/,
     &  '---------------------------------------',/)
7082        FORMAT(1x,'100: ',2x,I3,10x,I3,13x,I3)
            DO i=realtime(1),realtime(2),realtime(3)
              WRITE(nunit,7082,ERR=5040)i,channel(i,1),channel(i,2)
            END DO
            WRITE(nunit,*,ERR=5040)
          ENDIF
        ELSE
          IF (realtime(4) .EQ. 1) THEN
            WRITE(nunit,7083,ERR=5040)
     &      channel(realtime(1),1),channel(realtime(1),2)
7083        FORMAT(1x,'Tailfit start channel ',I3,
     &      '.',/,1x,'End channel ',I3,'.')
          ELSE
C           Realtime data
            WRITE(nunit,7084,ERR=5040)
7084        FORMAT(/,1x,'Frame    Start Channel    End Channel',/,
     &  '---------------------------------------',/)
7085        FORMAT(2x,I3,10x,I3,13x,I3)
            DO i=realtime(1),realtime(2),realtime(3)
              WRITE(nunit,7085,ERR=5040)i,channel(i,1),channel(i,2)
            END DO
            WRITE(nunit,*,ERR=5040)
          ENDIF
        ENDIF
        IF(l.EQ.1)THEN
          IF (backex .EQ. 'vonk') THEN
            WRITE(nunit,7090,ERR=5040)'Vonk.   '
          ELSE
            WRITE(nunit,7090,ERR=5040)'Guinier.'
          ENDIF
        ELSE
          IF (backex .EQ. 'vonk') THEN
            WRITE(nunit,7091,ERR=5040)'Vonk.   '
          ELSE
            WRITE(nunit,7091,ERR=5040)'Guinier.'
          ENDIF
        ENDIF
7090    FORMAT(1x,'100: Back extrapolation model: ',A8)
7091    FORMAT(1x,'Back extrapolation model: ',A8)
C       Fourier transform information.
        IF(l.EQ.1)THEN
          WRITE(nunit,7100,ERR=5040)
        ELSE
          WRITE(nunit,7101,ERR=5040)
        ENDIF
7100    FORMAT(//,1x,'100: FOURIER TRANSFORM PARAMETERS:')
7101    FORMAT(//,1x,'FOURIER TRANSFORM PARAMETERS:')
        IF(l.EQ.1)THEN
          WRITE(nunit,7110,ERR=5040)dmax,dstep
        ELSE
          WRITE(nunit,7111,ERR=5040)dmax,dstep
        ENDIF
7110    FORMAT(1x,'100: Transform performed up to D = ',E12.6,
     &' Angstroms.',/,1x,'100: Steps of ',E12.6,' Angstroms.')
7111    FORMAT(1x,'Transform performed up to D = ',E12.6,
     &' Angstroms.',/,1x,'Steps of ',E12.6,' Angstroms.')
7114    FORMAT(/,1x,'100: Frame',1x,I3,1x,'Start',1x,E12.6,
     &1x,'End',1x,E12.6)

C       End writing all the crap. Now the results.
C       Different output for single image.

        IF (realtime(4) .EQ. 1) THEN

C         Single frame.
          nframe=realtime(1)

C         First moments results
          IF(l.EQ.1)THEN
            CALL showprompts(prompts(23))
            WRITE(nunit,7120,ERR=5040)
          ELSE
            WRITE(nunit,7121,ERR=5040)
          ENDIF
7120      FORMAT(//,1x,'100: 1. MOMENTS RESULTS.',/)
7121      FORMAT(//,1x,'1. MOMENTS RESULTS.',/)
          IF (arbabs .EQ. 'ABS' .OR. arbabs .EQ. 'abs') THEN
            temp1=moment(nframe,1)*1.E+8
            temp2=moment(nframe,3)*1.E+24
            IF(l.EQ.1)THEN
              WRITE(nunit,7130,ERR=5040)temp1
            ELSE
              WRITE(nunit,7131,ERR=5040)temp1
            ENDIF
7130        FORMAT(1x,'100: Zeroth moment: ',E12.6,' cm^-2.')
7131        FORMAT(1x,'Zeroth moment: ',E12.6,' cm^-2.')
            IF(l.EQ.1)THEN
              WRITE(nunit,7140,ERR=5040)temp2
            ELSE
              WRITE(nunit,7141,ERR=5040)temp2
            ENDIF
7140        FORMAT(1x,'100: Invariant: ',E12.6,' cm^-4.')
7141        FORMAT(1x,'Invariant: ',E12.6,' cm^-4.')
          ELSE
            IF(l.EQ.1)THEN
              WRITE(nunit,7150,ERR=5040)moment(nframe,1)
            ELSE
              WRITE(nunit,7151,ERR=5040)moment(nframe,1)
            ENDIF
7150        FORMAT(1x,'100: Zeroth moment: ',E12.6,' [Arb.].')
7151        FORMAT(1x,'Zeroth moment: ',E12.6,' [Arb.].')
            IF(l.EQ.1)THEN
              WRITE(nunit,7160,ERR=5040)moment(nframe,3)
            ELSE
              WRITE(nunit,7161,ERR=5040)moment(nframe,3)
            ENDIF
7160        FORMAT(1x,'100: Invariant: ',E12.6,' [Arb.].')
7161        FORMAT(1x,'Invariant: ',E12.6,' [Arb.].')
          ENDIF

C         Next lamella results (if there are any).
          IF(l.EQ.1)THEN
            CALL showprompts(prompts(23))
            WRITE(nunit,7170,ERR=5040)
          ELSE
            WRITE(nunit,7171,ERR=5040)
          ENDIF
7170      FORMAT(//,1x,'100: 2. LAMELLA MORPHOLOGY RESULTS.',/)
7171      FORMAT(//,1x,'2. LAMELLA MORPHOLOGY RESULTS.',/)
C         Check whether cor. func. analysis successful.
          IF (reslam(nframe,1) .LT. 0.) THEN
C           No lamella results.
            WRITE(nunit,7180,ERR=5040)
7180        FORMAT(1x,'ERROR: Analysis of correlation function failed!'
     &      ,': NON-FATAL')
          ELSE
C           Output the whole batch of results.
7190        FORMAT(1x,'100: ',A60,E12.6,'.')
7191        FORMAT(1x,A60,E12.6,'.')
            IF(l.EQ.1)THEN
              text='Long period [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,1)
              text='Average hard block thickness [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,2)
              text='Average soft block thickness [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,3)
              text='Bulk volume crystallinity [No units]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,4)
              text='Local crystallinity [No units]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,5)
              text='Average hard block core thickness [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,6)
              text='Average interface thickness [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,7)
              text='Polydispersity [No units]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,8)
C             Check for abs units
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                text='Electron density contrast [cm^-2]'
                WRITE(nunit,7190,ERR=5040)
     &          text,reslam(nframe,9)
              ELSE
                text='Electron density contrast [Arb.]'
                WRITE(nunit,7190,ERR=5040)
     &          text,reslam(nframe,9)
              ENDIF
              text='Specific inner surface [Angstroms^-1]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,10)
              text='Non-ideality [No units]'
              WRITE(nunit,7190,ERR=5040)
     &        text,reslam(nframe,11)
            ELSE
              text='Long period [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,1)
              text='Average hard block thickness [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,2)
              text='Average soft block thickness [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,3)
              text='Bulk volume crystallinity [No units]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,4)
              text='Local crystallinity [No units]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,5)
              text='Average hard block core thickness [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,6)
              text='Average interface thickness [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,7)
              text='Polydispersity [No units]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,8)
C             Check for abs units
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                text='Electron density contrast [cm^-2]'
                WRITE(nunit,7191,ERR=5040)
     &          text,reslam(nframe,9)
              ELSE
                text='Electron density contrast [Arb.]'
                WRITE(nunit,7191,ERR=5040)
     &          text,reslam(nframe,9)
              ENDIF
              text='Specific inner surface [Angstroms^-1]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,10)
              text='Non-ideality [No units]'
              WRITE(nunit,7191,ERR=5040)
     &        text,reslam(nframe,11)
            ENDIF
C         End lamella results.
          ENDIF

C         Then Porod results.
C         Check where the crystallinity is from and which tail model was used.
          IF(l.EQ.1)THEN
            CALL showprompts(prompts(23))
            WRITE(nunit,7200,ERR=5040)
          ELSE
            WRITE(nunit,7201,ERR=5040)
          ENDIF
7200      FORMAT(//,1x,'100: 3. POROD RESULTS.',/)
7201      FORMAT(//,1x,'3. POROD RESULTS.',/)

          IF(l.EQ.1)THEN
            IF (sigmodel .EQ. 'on' .OR. sigmodel .EQ. 'ON') THEN
C             Tailfit used sigmoid model.
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                text='Porod constant [cm^-5]'
                WRITE(nunit,7190,ERR=5040)
     &          text,resporod(nframe,1)
              ELSE
                text='Porod constant [Arb.]'
                WRITE(nunit,7190,ERR=5040)
     &          text,resporod(nframe,1)
              ENDIF
              text='Sigma [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,resporod(nframe,2)
              IF (resporod(nframe,7) .GT. 0.) THEN
C               lamella morphology failed:-
                WRITE(6,7210)resporod(nframe,7)
7210            FORMAT(/,1x,'WARNING: Lamella-based interpretation of '
     &  ,'correlation function failed',/,1x,'100: User-supplied '
     &  ,'crystallinity suggests electron density contrast: ',E12.6,$)
C               Check abs units
                IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                  WRITE(nunit,7220,ERR=5040)' (cm^-2).'
                ELSE
                  WRITE(nunit,7220,ERR=5040)' (arb.).'
                ENDIF
7220            FORMAT(1x,'100: ',A10)
              ELSE
C               Lamella morphology OK - so no further results.
              ENDIF
            ELSE
C             Tailfit used a Porod model.
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                text='Porod constant [cm^-5]'
                WRITE(nunit,7190,ERR=5040)
     &          text,resporod(nframe,1)
              ELSE
                text='Porod constant [Arb.]'
                WRITE(nunit,7190,ERR=5040)
     &          text,resporod(nframe,1)
              ENDIF
              text='Porod chord length [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,resporod(nframe,3)
              IF (resporod(nframe,7) .GT. 0.) THEN
C               lamella model failed: user supplied cryst. used.
                WRITE(nunit,7230,ERR=5040)
7230  FORMAT(/,1x,'WARNING: Lamella-based interpretation of correlation'
     &,' function failed',/,1x,'100: Following results are based '
     &,'on the users crystallinity estimate.',/)
              ENDIF
              text='Crystalline chord length [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,resporod(nframe,4)
              text='Amorphous chord length [Angstroms]'
              WRITE(nunit,7190,ERR=5040)
     &        text,resporod(nframe,5)
              text='Surface to volume [Angstroms^-1]'
              WRITE(nunit,7190,ERR=5040)
     &        text,resporod(nframe,6)


              IF (resporod(nframe,7) .GT. 0.) THEN
C               lamella morphology failed:-
C               Check abs units
                IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                  text='Electron density contrast [cm^-2]'
                  WRITE(nunit,7190,ERR=5040)
     &            text,resporod(nframe,7)
                ELSE
                  text='Electron density contrast [Arb.]'
                  WRITE(nunit,7190,ERR=5040)
     &            text,resporod(nframe,7)
                ENDIF
              ENDIF

C           End choice of tailfit models
            ENDIF

          ELSE

            IF (sigmodel .EQ. 'on' .OR. sigmodel .EQ. 'ON') THEN
C             Tailfit used sigmoid model.
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                text='Porod constant [cm^-5]'
                WRITE(nunit,7191,ERR=5040)
     &          text,resporod(nframe,1)
              ELSE
                text='Porod constant [Arb.]'
                WRITE(nunit,7191,ERR=5040)
     &          text,resporod(nframe,1)
              ENDIF
              text='Sigma [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,resporod(nframe,2)
              IF (resporod(nframe,7) .GT. 0.) THEN
C               lamella morphology failed:-
                WRITE(6,7211)resporod(nframe,7)
7211            FORMAT(/,1x,'WARNING: Lamella-based interpretation of '
     &  ,'correlation function failed',/,1x,'User-supplied '
     &  ,'crystallinity suggests electron density contrast: ',E12.6,$)
C               Check abs units
                IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                  WRITE(nunit,7221,ERR=5040)' (cm^-2).'
                ELSE
                  WRITE(nunit,7221,ERR=5040)' (arb.).'
                ENDIF
7221            FORMAT(1x,A10)
              ELSE
C               Lamella morphology OK - so no further results.
              ENDIF
            ELSE
C             Tailfit used a Porod model.
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                text='Porod constant [cm^-5]'
                WRITE(nunit,7191,ERR=5040)
     &          text,resporod(nframe,1)
              ELSE
                text='Porod constant [Arb.]'
                WRITE(nunit,7191,ERR=5040)
     &          text,resporod(nframe,1)
              ENDIF

              text='Porod chord length [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,resporod(nframe,3)
              IF (resporod(nframe,7) .GT. 0.) THEN
C               lamella model failed: user supplied cryst. used.
                WRITE(nunit,7231,ERR=5040)
7231  FORMAT(/,1x,'WARNING: Lamella-based interpretation of correlation'
     &,' function failed',/,1x,'Following results are based '
     &,'on the users crystallinity estimate.',/)
              ENDIF
              text='Crystalline chord length [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,resporod(nframe,4)
              text='Amorphous chord length [Angstroms]'
              WRITE(nunit,7191,ERR=5040)
     &        text,resporod(nframe,5)
              text='Surface to volume [Angstroms^-1]'
              WRITE(nunit,7191,ERR=5040)
     &        text,resporod(nframe,6)


              IF (resporod(nframe,7) .GT. 0.) THEN
C               lamella morphology failed:-
C               Check abs units
                IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                  text='Electron density contrast [cm^-2]'
                  WRITE(nunit,7191,ERR=5040)
     &            text,resporod(nframe,7)
                ELSE
                  text='Electron density contrast [Arb.]'
                  WRITE(nunit,7191,ERR=5040)
     &            text,resporod(nframe,7)
                ENDIF
              ENDIF

C           End choice of tailfit models
            ENDIF

          ENDIF
C       More than one frame.
        ELSE

C         First moments results
          IF(l.EQ.1)THEN
            CALL showprompts(prompts(23))
            WRITE(nunit,7120,ERR=5040)
          ELSE
            WRITE(nunit,7121,ERR=5040)
          ENDIF
C         Write column headings depending on arb/abs units
          IF (arbabs .EQ. 'ABS' .OR. arbabs .EQ. 'abs') THEN
            IF(l.EQ.1)THEN
              WRITE(nunit,7240,ERR=5040)
            ELSE
              WRITE(nunit,7241,ERR=5040)
            ENDIF
7240        FORMAT(1x,'100: Frame    Moment(0) [cm^-2]    Moment(2)'
     &,' [cm^-4]'
     &,/,'-------------------------------------------------',/)
7241        FORMAT(1x,'Frame    Moment(0) [cm^-2]    Moment(2)'
     &,' [cm^-4]'
     &,/,'-------------------------------------------------',/)
            DO nframe=realtime(1),realtime(2),realtime(3)
C             Scale up
              moment(nframe,1)=moment(nframe,1)*1.E+8
              moment(nframe,3)=moment(nframe,3)*1.E+24
            END DO
          ELSE
            IF(l.EQ.1)THEN
              WRITE(nunit,7250,ERR=5040)
            ELSE
              WRITE(nunit,7251,ERR=5040)
            ENDIF
7250        FORMAT(1x,'100: Frame        Moment(0)            Moment(2)'
     &      ,/,'-------------------------------------------------',/)
7251        FORMAT(1x,'Frame        Moment(0)            Moment(2)'
     &      ,/,'-------------------------------------------------',/)
          ENDIF

          DO nframe=realtime(1),realtime(2),realtime(3)
            IF(l.EQ.1)THEN
              WRITE(nunit,7260,ERR=5040)nframe,
     &        moment(nframe,1),moment(nframe,3)
            ELSE
              WRITE(nunit,7261,ERR=5040)nframe,
     &        moment(nframe,1),moment(nframe,3)
            ENDIF
          END DO
7260      FORMAT(1x,'100: ',2x,I3,7x,E12.6,9x,E12.6)
7261      FORMAT(2x,I3,7x,E12.6,9x,E12.6)


C         Next lamella results (if there are any).
          IF(l.EQ.1)THEN
            CALL showprompts(prompts(23))
            WRITE(nunit,7170,ERR=5040)
          ELSE
            WRITE(nunit,7171,ERR=5040)
          ENDIF

C         Do block thicknesses and local cryst first
          IF(l.EQ.1)THEN
            WRITE(nunit,7270,ERR=5040)
          ELSE
            WRITE(nunit,7271,ERR=5040)
          ENDIF
7270      FORMAT(1x,'100: Frame    Long period    Av.Hard Block [A]'
     &,'    Av.Soft Block [A]    Local Cryst.',/,'-------'
     &,'--------------------------------------'
     &,'----------------------------------',/)
7271      FORMAT(1x,'Frame    Long period    Av.Hard Block [A]'
     &,'    Av.Soft Block [A]    Local Cryst.',/,'-------'
     &,'--------------------------------------'
     &,'----------------------------------',/)
          DO nframe=realtime(1),realtime(2),realtime(3)
            IF (reslam(nframe,1) .GT. 0.) THEN
              IF(l.EQ.1)THEN
                WRITE(nunit,7280,ERR=5040)nframe,reslam(nframe,1),
     &          reslam(nframe,2),reslam(nframe,3),reslam(nframe,5)
              ELSE
                WRITE(nunit,7281,ERR=5040)nframe,reslam(nframe,1),
     &          reslam(nframe,2),reslam(nframe,3),reslam(nframe,5)
              ENDIF
            ELSE
              WRITE(nunit,7290,ERR=5040)nframe
            ENDIF
          END DO
7280      FORMAT(1x,'100: ',2x,I3,5x,E12.6,5x,E12.6,9x,E12.6,7x,E12.6)
7281      FORMAT(2x,I3,5x,E12.6,5x,E12.6,9x,E12.6,7x,E12.6)
7290      FORMAT(1x,'WARNING: ',2x,I3,5x,
     &'Lamella interpretation failed')

C         Then do interfaces and bulk cryst.
          IF(l.EQ.1)THEN
            CALL showprompts(prompts(23))
            WRITE(nunit,7300,ERR=5040)
          ELSE
            WRITE(nunit,7301,ERR=5040)
          ENDIF
7300      FORMAT(//,1x,'100: Frame    Bulk Cryst.    Av.Core '
     &,'Thickness [A]    Av.Interface '
     &,'Thickness [A]',/,'-----------------------------'
     &,'-----------------------------------------------',/)
7301      FORMAT(//,1x,'Frame    Bulk Cryst.    Av.Core '
     &,'Thickness [A]    Av.Interface '
     &,'Thickness [A]',/,'-----------------------------'
     &,'-----------------------------------------------',/)
          DO nframe=realtime(1),realtime(2),realtime(3)
            IF (reslam(nframe,1) .GT. 0.) THEN
              IF(l.EQ.1)THEN
                WRITE(nunit,7310,ERR=5040)nframe,reslam(nframe,4),
     &          reslam(nframe,6),reslam(nframe,7)
              ELSE
                WRITE(nunit,7311,ERR=5040)nframe,reslam(nframe,4),
     &          reslam(nframe,6),reslam(nframe,7)
              ENDIF
            ELSE
              WRITE(nunit,7290,ERR=5040)nframe
            ENDIF
          END DO
7310      FORMAT(1x,'100: ',2x,I3,5x,E12.6,7x,E12.6,16x,E12.6)
7311      FORMAT(2x,I3,5x,E12.6,7x,E12.6,16x,E12.6)

C         Finally do the other things
          IF(l.EQ.1)THEN
            CALL showprompts(prompts(23))
            IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
              WRITE(nunit,7330,ERR=5040)
            ELSE
              WRITE(nunit,7320,ERR=5040)
            ENDIF
          ELSE
            IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
              WRITE(nunit,7331,ERR=5040)
            ELSE
              WRITE(nunit,7321,ERR=5040)
            ENDIF
          ENDIF
7320  FORMAT(//,1x,'100: Frame    Polydisp.   Elec.Dens.Contr.[Arb.]'
     &,'  Spec.Inner.Surf.[1/A]  Non-ideal',/,'---------------'
     &,'------------------'
     &,'-----------------------------------------------',/)
7330  FORMAT(//,1x,'100: Frame    Polydisp.   Elec.Dens.Contr.[cm^-2]'
     &,' Spec.Inner.Surf.[1/A]  Non-ideal',/,'---------------'
     &,'------------------'
     &,'-----------------------------------------------',/)
7321  FORMAT(//,1x,'Frame    Polydisp.   Elec.Dens.Contr.[Arb.]'
     &,'  Spec.Inner.Surf.[1/A]  Non-ideal',/,'---------------'
     &,'------------------'
     &,'-----------------------------------------------',/)
7331  FORMAT(//,1x,'Frame    Polydisp.   Elec.Dens.Contr.[cm^-2]'
     &,' Spec.Inner.Surf.[1/A]  Non-ideal',/,'---------------'
     &,'------------------'
     &,'-----------------------------------------------',/)
          DO nframe=realtime(1),realtime(2),realtime(3)
            IF (reslam(nframe,1) .GT. 0.) THEN
              IF(l.EQ.1)THEN
                WRITE(nunit,7340,ERR=5040)nframe,reslam(nframe,8),
     &          reslam(nframe,9),reslam(nframe,10),reslam(nframe,11)
              ELSE
                WRITE(nunit,7341,ERR=5040)nframe,reslam(nframe,8),
     &          reslam(nframe,9),reslam(nframe,10),reslam(nframe,11)
              ENDIF
            ELSE
              WRITE(nunit,7290,ERR=5040)nframe
            ENDIF
          END DO
7340      FORMAT(1x,'100: ',2x,I3,3x,E12.6,7x,E12.6,11x,E12.6,6x,E12.6)
7341      FORMAT(2x,I3,3x,E12.6,7x,E12.6,11x,E12.6,6x,E12.6)

C         End lamella results.

C         Then the Porod results - a right old dog's dinner (sorry).
          IF(l.EQ.1)THEN
            CALL showprompts(prompts(23))
            WRITE(nunit,7200,ERR=5040)
          ELSE
            WRITE(nunit,7201,ERR=5040)
          ENDIF
C         Check which model.
          IF (sigmodel .EQ. 'on' .OR. sigmodel .EQ. 'ON') THEN
C           Tailfit used sigmoid model.
            IF(l.EQ.1)THEN
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                WRITE(nunit,7350,ERR=5040)
              ELSE
                WRITE(nunit,7360,ERR=5040)
              ENDIF
            ELSE
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                WRITE(nunit,7351,ERR=5040)
              ELSE
                WRITE(nunit,7361,ERR=5040)
              ENDIF
            ENDIF
7350  FORMAT(1x,'100: Frame    Porod const.[cm^-5]    Sigma [A]    '
     &,'Cryst.    Elec.Dens.Contr.[cm^-2]',/,
     &'---------------------------------------------'
     &,'----------------------------------',/)
7360  FORMAT(1x,'100: Frame    Porod const. [Arb.]    Sigma [A]    '
     &,'Cryst.    Elec.Dens.Contr. [Arb.]',/,
     &'---------------------------------------------'
     &,'----------------------------------',/)
7351  FORMAT(1x,'Frame    Porod const.[cm^-5]    Sigma [A]    '
     &,'Cryst.    Elec.Dens.Contr.[cm^-2]',/,
     &'---------------------------------------------'
     &,'----------------------------------',/)
7361  FORMAT(1x,'Frame    Porod const. [Arb.]    Sigma [A]    '
     &,'Cryst.    Elec.Dens.Contr. [Arb.]',/,
     &'---------------------------------------------'
     &,'----------------------------------',/)

            DO nframe=realtime(1),realtime(2),realtime(3)
              IF (resporod(nframe,7) .GT. 0.) THEN
C             lamella morphology failed:-
                IF(l.EQ.1)THEN
                  WRITE(nunit,7370,ERR=5040)nframe,resporod(nframe,1),
     &            resporod(nframe,2),'  USER ',resporod(nframe,7)
                ELSE
                  WRITE(nunit,7371,ERR=5040)nframe,resporod(nframe,1),
     &            resporod(nframe,2),'  USER ',resporod(nframe,7)
                ENDIF
7370  FORMAT(1x,'100: ',2x,I3,8x,E12.6,6x,E12.6,3x,A7,8x,E12.6)
7371  FORMAT(2x,I3,8x,E12.6,6x,E12.6,3x,A7,8x,E12.6)
              ELSE
                IF(l.EQ.1)THEN
                  WRITE(nunit,7380,ERR=5040)nframe,resporod(nframe,1),
     &            resporod(nframe,2),'LAMELLA'
                ELSE
                  WRITE(nunit,7381,ERR=5040)nframe,resporod(nframe,1),
     &            resporod(nframe,2),'LAMELLA'
                ENDIF
7380  FORMAT(1x,'100: ',2x,I3,8x,E12.6,6x,E12.6,2x,A7)
7381  FORMAT(2x,I3,8x,E12.6,6x,E12.6,2x,A7)
              ENDIF
            END DO

          ELSE
C           Porod model - more complicated results.

C           First output the things that don't depend on cryst.
            IF(l.EQ.1)THEN
              WRITE(nunit,7390,ERR=5040)
            ELSE
              WRITE(nunit,7391,ERR=5040)
            ENDIF
7390        FORMAT(1x,'100: Results independent of crystallinity:-')
7391        FORMAT(1x,'Results independent of crystallinity:-')
            IF(l.EQ.1)THEN
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                WRITE(nunit,7400,ERR=5040)
              ELSE
                WRITE(nunit,7410,ERR=5040)
              ENDIF
            ELSE
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                WRITE(nunit,7401,ERR=5040)
              ELSE
                WRITE(nunit,7411,ERR=5040)
              ENDIF
            ENDIF
7400  FORMAT(/,1x,'100: Frame    Porod const.[cm^-5]    '
     &,'Chord Length [A]',/,
     &'---------------------------'
     &,'---------------------',/)
7410  FORMAT(/,1x,'100: Frame    Porod const. [Arb.]    '
     &,'Chord Length [A]',/,
     &'---------------------------'
     &,'---------------------',/)
7401  FORMAT(/,1x,'Frame    Porod const.[cm^-5]    '
     &,'Chord Length [A]',/,
     &'---------------------------'
     &,'---------------------',/)
7411  FORMAT(/,1x,'Frame    Porod const. [Arb.]    '
     &,'Chord Length [A]',/,
     &'---------------------------'
     &,'---------------------',/)

            DO nframe=realtime(1),realtime(2),realtime(3)
              IF(l.EQ.1)THEN
                WRITE(nunit,7420,ERR=5040)nframe,resporod(nframe,1),
     &          resporod(nframe,3)
              ELSE
                WRITE(nunit,7421,ERR=5040)nframe,resporod(nframe,1),
     &          resporod(nframe,3)
              ENDIF
7420          FORMAT(1x,'100: ',1x,I3,8x,E12.6,11x,E12.6)
7421          FORMAT(1x,I3,8x,E12.6,11x,E12.6)
            END DO

C           Then output the things that do depend on crystallinity.
            IF(l.EQ.1)THEN
              CALL showprompts(prompts(23))
              WRITE(nunit,7430,ERR=5040)
            ELSE
              WRITE(nunit,7431,ERR=5040)
            ENDIF
7430        FORMAT(//,1x,'100: Results dependent on crystallinity:')
7431        FORMAT(//,1x,'Results dependent on crystallinity:')

C           Check abs units
            IF(l.EQ.1)THEN
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                WRITE(nunit,7440,ERR=5040)
              ELSE
                WRITE(nunit,7450,ERR=5040)
              ENDIF
            ELSE
              IF (arbabs .EQ. 'abs' .OR. arbabs .EQ. 'ABS') THEN
                WRITE(nunit,7441,ERR=5040)
              ELSE
                WRITE(nunit,7451,ERR=5040)
              ENDIF
            ENDIF
7440  FORMAT(/,1x,'100: Frame  Cryst  HardChord[A]  SoftChord[A]  '
     &,'Surf:Vol[1/A]  Elec.Dens.Contr.[cm^-2]',/,
     &'------------------------------------------------'
     &,'----------------------------------',/)
7450  FORMAT(/,1x,'100: Frame  Cryst  HardChord[A]  SoftChord[A]  '
     &,'Surf:Vol[1/A]  Elec.Dens.Contr.[Arb.]',/,
     &'------------------------------------------------'
     &,'----------------------------------',/)
7441  FORMAT(/,1x,'Frame  Cryst  HardChord[A]  SoftChord[A]  '
     &,'Surf:Vol[1/A]  Elec.Dens.Contr.[cm^-2]',/,
     &'------------------------------------------------'
     &,'----------------------------------',/)
7451  FORMAT(/,1x,'Frame  Cryst  HardChord[A]  SoftChord[A]  '
     &,'Surf:Vol[1/A]  Elec.Dens.Contr.[Arb.]',/,
     &'------------------------------------------------'
     &,'----------------------------------',/)

            DO nframe=realtime(1),realtime(2),realtime(3)
              IF (resporod(nframe,7) .GT. 0.) THEN
C               Lamella interpretation failed
                IF(l.EQ.1)THEN
                  WRITE(nunit,7460,ERR=5040)nframe,'USER',
     &            resporod(nframe,4),resporod(nframe,5),
     &            resporod(nframe,6),resporod(nframe,7)
                ELSE
                  WRITE(nunit,7461,ERR=5040)nframe,'USER',
     &            resporod(nframe,4),resporod(nframe,5),
     &            resporod(nframe,6),resporod(nframe,7)
                ENDIF
              ELSE
C               Lamella interpretation worked
                IF(l.EQ.1)THEN
                  WRITE(nunit,7470,ERR=5040)nframe,'LAMLA',
     &            resporod(nframe,4),resporod(nframe,5),
     &            resporod(nframe,6)
                ELSE
                  WRITE(nunit,7471,ERR=5040)nframe,'LAMLA',
     &            resporod(nframe,4),resporod(nframe,5),
     &            resporod(nframe,6)
                ENDIF
              ENDIF
            END DO
7460  FORMAT(1x,'100: ',I3,3x,A5,2x,E12.6,2x,E12.6,2x,E12.6,8x,E12.6)
7470  FORMAT(1x,'100: ',I3,3x,A5,2x,E12.6,2x,E12.6,2x,E12.6)
7461  FORMAT(1x,I3,3x,A5,2x,E12.6,2x,E12.6,2x,E12.6,8x,E12.6)
7471  FORMAT(1x,I3,3x,A5,2x,E12.6,2x,E12.6,2x,E12.6)

C         End choice of model
          ENDIF

C       End check on single frame or realtime.
        ENDIF

C       End results file.
        IF(l.EQ.1)THEN
          CALL showprompts(prompts(23))
          WRITE(nunit,7480,ERR=5040)
        ELSE
          WRITE(nunit,7481,ERR=5040)
        ENDIF

7480    FORMAT(//,1x,'100: END OF CORRELATION FUNCTION RESULTS.')
7481    FORMAT(//,1x,'END OF CORRELATION FUNCTION RESULTS.')

C     End loop through output units
      END DO

      CLOSE(nout(2))

C     The end of the entire session...
      CALL showprompts(prompts(14))
      WRITE(6,*)

C     End program.
      CALL showprompts(prompts(20))
      STOP


C     Error Catches:

C     Error reading corfunc.txt
5000  CALL showprompts(prompts(5))
      STOP

C     Error reading extract.txt
5010  CALL showprompts(prompts(1))
      STOP

C     Error reading in correlation function otoko files.
5030  CALL showprompts(prompts(6))
      STOP

C     Error writing result file.
5040  CALL showprompts(prompts(11))
      STOP

      STOP
      END



      SUBROUTINE diffpoly(polynom)
C     Differentiates the polynomial whose coefficients are stored in
C     the array "polynom" algebraically.

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

C      DIMENSION polynom(MaxDim,5)
      REAL*4 :: polynom(MaxDim,5)

C     Loop through polynomials.
      DO i=1,MaxDim
C       Loop through coeffs
        DO j=1,4
C         Differentiate:
          polynom(i,j)=polynom(i,j+1)*FLOAT(j)
        END DO
        polynom(i,5)=0.
      END DO
      RETURN
      END



      SUBROUTINE evalpoly(polynom,x,y)
C     Evaluates the polynomial stored in the array "polynom"
C     at the point x.

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

C      DIMENSION polynom(MaxDim,5),x(MaxDim),y(MaxDim)
      REAL*4 :: polynom(MaxDim,5)
      REAL*4 :: x(MaxDim)
      REAL*4 :: y(MaxDim)

C     Loop through points
      DO npoint=1,MaxDim
        xx=x(npoint)
C       Evaluate polynomial.
        sum=0.
C       Loop through powers of x
        DO i=1,5
          IF (polynom(npoint,i) .NE. 0.) THEN
            sum=sum+polynom(npoint,i)*(xx**(i-1))
          ENDIF
        END DO
        y(npoint)=sum
      END DO
      RETURN
      END



      SUBROUTINE laginterp(xdata,ydata,nstart,nend,polynom)
C     Fits a 4th degree polynomial to the data points
C     specified by nstart and nend. The coefficients of the
C     polynomial are placed in the array "polynom".
C     Traditional Lagrangian interpolation.
C     Reference: any book on numerical analysis.

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

C      DIMENSION xdata(MaxDim),ydata(MaxDim),polynom(MaxDim,5),coeff(5)
      REAL*4 :: xdata(MaxDim)
      REAL*4 :: ydata(MaxDim)
      REAL*4 :: polynom(MaxDim,5)
      REAL*4 :: coeff(5)
      INTEGER nstart,nend

1000  FORMAT(/,1x,'ERROR: Error in subroutine laginterp: FATAL')

C     get number of points
      npts=nend-nstart+1

C     Check for errors
      IF (npts .LE. 0) THEN
        WRITE(6,1000)
        STOP
      ENDIF

C     Loop through points
      DO npoint=nstart,nend

C       Reset polynomial coefficients.
        DO i=1,5
          polynom(npoint,i)=0.
        END DO

C       Loop through the terms in the Lagrangian "formula".
        DO nterm=npoint-2,npoint+2

C         Get denominator.
          denom=1.
          DO i=npoint-2,npoint+2
            IF (i .NE. nterm) THEN
              factor=xdata(nterm)-xdata(i)
              denom=denom*factor
            ENDIF
          END DO

C         Reset  the array "coeff"
          DO i=1,5
            coeff(i)=0.
          END DO

C         Initialise "coeff"
          coeff(1)=1.

C         Loop through factors in the numerator, multiplying the polynomial
C         effectively stored in "coeff" by each factor.
          DO i=npoint-2,npoint+2
C           Check to see whether the current factor is actually in this term.
            IF (i .NE. nterm) THEN
C             Multiply by factor: update "coeff".
C             Multiply by x
              DO j=5,2,-1
                coeff(j)=coeff(j-1)
              END DO
              coeff(1)=0.
C             Multiply by constant term.
              DO j=1,4
                coeff(j)=coeff(j)-coeff(j+1)*xdata(i)
              END DO
            ENDIF
          END DO

C         Now add the polynomial stored in "coeff" onto "polynom".
          DO i=1,5
            polynom(npoint,i)=polynom(npoint,i)+
     &      coeff(i)*ydata(nterm)/denom
          END DO

C       End loop through terms in the formula.
        END DO

C     End loop through points.
      END DO

C     We're done!

      RETURN
      END


C=============================================================================
C==============================================================================


      SUBROUTINE TROPUS4

*********************************************************************************
*                                                                               *
*     PROGRAM TO INVERT I0(Q) TERM FROM SMALL-ANGLE NEUTRON SCATTERING DATA     *
*                                                                               *
*********************************************************************************
*                                                                               *
*     Based on program TLPHAS originally written for PDP/TSX                    *
*     (C)1984 Trevor Crowley, University of Bristol, UK                         *
*                                                                               *
*     Modified for use on VAX/VMS (Program TROPUS)                              *
*     (C)Terry Cosgrove, University of Bristol, UK                              *
*     Modified for use on IBM PCs (Program PCPHASE)                             *
*     (C)Terry Cosgrove, University of Bristol, UK                              *
*     PCPHASE revised for VAX/VMS (Program TROPUS2)                             *
*     (C)Steve King, ISIS Facility, Rutherford Appleton Laboratory, UK          *
*     TROPUS2 made a GENIE function (Program TROPUS3)                           *
*     (C)Steve King, ISIS Facility, Rutherford Appleton Laboratory, UK          *
*     TROPUS3 made a CORFUNC subroutine (Program TROPUS4)                       *
*     (C)Steve King, ISIS Facility, Rutherford Appleton Laboratory, UK          *
*                                                                               *
*********************************************************************************
*
*     S King, July 2004
*     Array bounds raised from 2048 --> 4096 because CORFUNC extrapolations have
*     2048 points and program falls over when filling yvec otherwise!
*
*********************************************************************************

      INTEGER      MaxDim
      PARAMETER    (MaxDim=4096)

      real*4       rescale2z,density,gamma,norm,atx,xstep,yis,dyis
      real*4       integral,rescale,delta,change
      real*4       spanint,span,boundfrac,z1,z2,rms,sigma
      real*8       dx,dk,rotx,rotk,ri
      integer*4    n,nqdata,spanchan,npower,j,limitq
      character*1  qisang
      character*3  ident
      character*80 :: prompts(7)
      character*80 fname

      real*4 :: work(MaxDim)
      real*4 :: x(MaxDim)
      real*4 :: y(MaxDim)
      real*4 :: e(MaxDim)
      real*4 :: xx(MaxDim)
      real*4 :: yy(MaxDim)
      real*8 :: yvec(2,MaxDim)

      common/fft/yvec,n0,m0

      prompts(1)='ERROR: Files have different numbers of channels: FATAL'
      prompts(2)='ERROR: More than 4096 channels in extrapolation: FATAL'
      prompts(3)='ERROR: 0.5 < density < 2.0: FATAL'
      prompts(4)='ERROR: 0 < adsorbed amount < 100: FATAL'
      prompts(5)='ERROR: 0.0 < difference < 1.0: FATAL'
      prompts(6)='ERROR: Axis units not specified.  Cannot normalise: FATAL'
      prompts(7)='WARNING: Unusual normalisation value.  Will use 1.0'

*     OPEN JOURNAL FILE
      open(1,file='tropus.txt',form='formatted',status='unknown')

*==============================================================================
*     GET DATA & PARAMETERS
*==============================================================================

*1000     format(' Extrapolated DATA filename [???FU2.TXT - enter first 3 characters]? ',$)
1000     format(' Extrapolated DATA filename [???FUL.TXT - enter first 3 characters]? ',$)
1001     format(a3,$)
1002     format(' Reading file: ',a10)
1003     format(' Data has',i5,' channels')
1004     format(' Is the Q-axis in units of PER ANGSTROM ........................[y]? ',$)
1005     format(a1,$)
1006     format(' Bulk density of adsorbed polymer .........................[g/cm^3]? ',$)
1007     format(' Adsorbed amount ..........................................[mg/m^2]? ',$)
1008     format(f12.5,$)
1009     format(' Normalisation factor is',f8.3)
1010     format(' Profile can differ from normalisation by what amount .............? ',$)
1015     format(' Integral is: ',f12.5,'      Difference from Normalisation factor: ',f12.5)
1016     format(' (NORM-INTEGRAL): ',f12.5,' RESCALE: ',f12.5,' CHANGE: ',f12.5)
1020     format(' Writing file: ',a10)

10    write(6,1000)
      read(5,1001,err=10)ident
C     SMK, 19/08/04
C     The CORFUNC FU2 file has a background added to it, the FUL file doesn't
*     fname=ident//'FU2.TXT'
      fname=ident//'FUL.TXT'
      write(6,1002)fname
      write(1,1002)fname
      call load_data(fname,n,y)

      fname=ident//'FLX.TXT'
      write(6,1002)fname
      write(1,1002)fname
      call load_data(fname,nqdata,x)

20     qisang='y'
      write(6,1004)
      read(5,1005,err=20)qisang
      if ((qisang.eq.' ').or.(qisang.eq.'y').or.(qisang.eq.'Y')) then
        rescale2z=10.0
      else if ((qisang.eq.'n').or.(qisang.eq.'N')) then
        rescale2z=1.0
      else
        rescale2z=0.0
      end if
      if (rescale2z.eq.0.0) call showprompts(prompts(6))
      if (rescale2z.eq.0.0) goto 9999
      write(1,'(a2)')qisang

30    write(6,1006)
      read(5,1008,err=30)density
      if ((density.gt.0.5).and.(density.lt.2.0)) goto 40
      call showprompts(prompts(3))
      goto 9999

40    write(6,1007)
      read(5,1008,err=40)gamma
      if ((gamma.gt.0.0).and.(gamma.lt.100.0)) goto 50
      call showprompts(prompts(4))
      goto 9999

*     DENSITY IS INPUT IN G/CM^3, GAMMA IS INPUT IN MG/M^2
*     WHICH MEANS THAT NORM=(GAMMA/DENSITY) IS A LENGTH; x10^-7 CM
*     SO HAVE TO SCALE BY x10 IF Z-AXIS IS ANGSTROMS
50    if ((qisang.eq.' ').or.(qisang.eq.'y').or.(qisang.eq.'Y')) then
        norm=10.0*gamma/density
      else if ((qisang.eq.'n').or.(qisang.eq.'N')) then
        norm=gamma/density
      else
        norm=1.0
        call showprompts(prompts(7))
      end if
      write(6,1009)norm
      write(1,'(f12.5)')density
      write(1,'(f12.5)')gamma
      write(1,1009)norm

60    write(6,1010)
      read(5,1008,err=60)delta
      if ((delta.gt.0.0).and.(delta.lt.1.0)) goto 70
      call showprompts(prompts(5))
      goto 9999

70    write(1,'(f12.5)')delta
      continue

*==============================================================================
*     PREPROCESSING OF DATA
*==============================================================================

1     if (nqdata.ne.n) goto 8001
      if (n.gt.4096) goto 8002
*     LET'S MAKE THINGS MORE MANAGEABLE!
      n=n/2
      if (mod(n,2).ne.0) n=n+1
      write(6,1003)n
      write(1,1003)n
      write(6,*)'Interpolating for equally-spaced q-intervals...'
*     INTERPOLATE THE DATA SET
*     S King, August 2004
*
*     Tropus requires that the input data has a q-axis with equally-spaced
*     points.  The data just read in though is only equally-spaced in each
*     of the low-q, original, and tail-fitted data ranges.  So need to run
*     an interpolation function through them all in one go.
*
*     But can't use the Savitsky-Golay routine as Corfunc does because
*     it is designed for equally-spaced X-data!  So here we use a generic
*     polynomial interpolation from Numerical Recipies.
*     Thanks to Richard Heenan for this suggestion.
*
*     XSTEP IS EQUALLY-SPACED X INCREMENT
      xstep=(x(n)-x(1))/float(n-1)
      atx=x(1)
*     NPOWER=4 FITS A CUBIC THROUGH 4 POINTS
*     NPOWER SHOULD BE AN EVEN NUMBER
      npower=4
      j=1
      do i=1,n
         do while ((j.lt.(n-npower+1)).and.(atx.gt.x(j+npower/2)))
            j=j+1
         end do
         yis=0.
         dyis=0.
         e(i)=0.
*        POLINT USES NPOWER CONSECUTIVE POINTS OF X()
*        I.E. X(J) TO X(J+NPOWER-1)
*        X() IS ASSUMED INCREASING BUT NOT IN EQUAL INTERVALS
*        ATX SHOULD BE BETWEEN THE MIDDLE TWO POINTS, **EXCEPT**
*        AT THE START AND END OF X()
         call polint(x(j),y(j),npower,atx,yis,dyis)
         xx(i)=atx
         yy(i)=yis
         e(i)=dyis
         atx=atx+xstep
      end do
*     NOW WRITE THE INTERPOLATED VALUES BACK INTO THE INITIAL DATA ARRAYS
      do i=1,n
        x(i)=xx(i)
        y(i)=yy(i)
        xx(i)=0.
        yy(i)=0.
      end do
*     WRITE INTERPOLATED DATA
      fname=ident//'INT.TXT'
      write(6,1020)fname
      write(1,1020)fname
      call save_data(fname,n,x,y)
*     SET ANY NEGATIVE OR ZERO INTENSITIES TO A VERY SMALL POSITIVE VALUE
      do i=1,n
        if (y(i).le.0.) y(i)=1.0e-08
        if (y(i).le.0.) then
           write(6,*)'Setting negative and zero intensities to 0.00000001...'
        endif
      end do
      write(6,*)'Multiplying intensities by Q-squared...'
*     MULTIPLY INTENSITES BY Q-SQUARED
      do i=1,n
        y(i)=y(i)*(x(i)*x(i))
      end do

*==============================================================================
*     MAIN PROGRAM
*==============================================================================

*     SHIFT FACTOR ROTX SET EQUAL TO ZERO FOR HALF POINT FORMAT
1030  n0=2*n
      n2=n-1
      rotk=0.d0
      rotx=-0.5d0
**     write(6,*)'Setting up function parameters...'
*     SET UP FUNCTION PARAMETERS
*     READ FTC ZERO PT.
      dk=dble((x(n)-x(1))/dfloat(n-1))
*     M Rodman, April 2004
*     INITIALISE YVEC ARRAY  
      do 1035 ll=1,MaxDim
        yvec(1,ll)=0d0
        yvec(2,ll)=0d0
1035  continue
*     CENTRAL POINT IN YVEC(N0/2+1) FROM Y(1), IE WHERE Q=0
      yvec(1,n+1)=dble(y(1))
      yvec(2,n+1)=0.d0
*     FIRST POINT YVEC(1) LEFT OVER.  TRY PUTTING EQUAL TO YVEC(2), 0.0, etc
      yvec(1,1)=dble(y(n))
      yvec(2,1)=0.d0
*     REST OF YVEC FILLED FROM Y BY REFLECTION ABOUT Q=0
      do 1040 ll=2,n
        yvec(1,n+ll)=dble(y(ll))
        yvec(2,n+ll)=0.d0
        yvec(1,n-ll+2)=dble(y(ll))
        yvec(2,n-ll+2)=0.d0
1040     continue
**     write(6,*)'Calling rephase...'
      call rephase(dk,dx,rotx)
      do 1050 ll=1,n
        xx(ll)=sngl(float(ll-1)*dx-rotx*dx)
*        yy(ll)=sngl(yvec(1,ll+n))
        yy(ll)=abs(sngl(yvec(1,ll+n)))
1050  continue

*==============================================================================
*     POST PROCESSING OF TRANSFORMED DATA
*==============================================================================

*     RESCALE SEGMENT DENSITY TO THE Z-AXIS UNITS
*     x10 IF PER ANGSTROM
*      x1 IF PER NANOMETRE
      write(6,*)'Rescaling segment density to the z-axis units...'
      do 1060 ll=1,n
        yy(ll)=yy(ll)*rescale2z
1060     continue
*     OUTPUT SEGMENT DENSITY PROFILE
      fname=ident//'SDP.TXT'
      write(6,1020)fname
      write(1,1020)fname
      call save_data(fname,n,xx,yy)
*     THEN RESCALE SO THAT THE AREA UNDER THE DISTRIBUTION
*     MATCHES THE KNOWN ADSORBED AMOUNT
      write(6,*)'Integrating under distribution...'
      rescale=1.0
      change=1.0
1070  integral=0.0
      call trapiz(xx,yy,n,1,integral)
*     write(1,1016)norm-integral,rescale,change
      if (abs(norm-integral).lt.delta) goto 1090
      if ((norm-integral).gt.delta) then
        rescale=rescale*(1+change)
        change=rescale/1.61803399
      else if ((norm-integral).lt.delta) then
        rescale=rescale/(1+change)
        change=rescale*3.14159265
      else
        rescale=-1.0
      end if
      if (rescale.eq.-1.0) goto 1090
      do 1080 ll=1,n
        yy(ll)=yy(ll)*rescale
1080     continue
      goto 1070
1090     write(6,*)'Normalised integral to adsorbed amount...'
      write(1,1015)integral,(norm-integral)

*==============================================================================
*     EXTRACT PARAMETERS
*==============================================================================

*     EXTENT OF PROFILE
      span=0.0
      spanchan=n
1100  spantint=0.0
      call trapiz(xx,yy,spanchan,1,spanint)
      if (spanint.le.(0.95*integral)) goto 1110
      spanchan=spanchan-1
      goto 1100
1110  span=xx(spanchan)

*     BOUND FRACTION
      if ((qisang.eq.'y').or.(qisang.eq.'Y').or.(qisang.eq.' ')) then
        ll=1
1120    if (xx(ll).gt.10.0) goto 1140
        ll=ll+1
        goto 1120
      else if ((qisang.eq.'n').or.(qisang.eq.'n')) then
        ll=1
1130    if (xx(ll).gt.1.0) goto 1140
        ll=ll+1
        goto 1130
      end if
1140  boundfrac=0.0
      call trapiz(xx,yy,ll,1,boundfrac)
      boundfrac=boundfrac/integral

*     MOMENT-WEIGHTED DISTRIBUTIONS
*     FOR <z^1>
      do ll=1,n
        work(ll)=yy(ll)*xx(ll)
      end do
      z1=0.0
      call trapiz(xx,work,spanchan,1,z1)
      z1=z1/integral
*     FOR <z^2>
      do ll=1,n
        work(ll)=yy(ll)*xx(ll)*xx(ll)
      end do
      z2=0.0
      call trapiz(xx,work,spanchan,1,z2)
      z2=z2/integral
      rms=sqrt(z2)
      sigma=sqrt(z2-(z1*z1))

7000  format(' Bound fraction ..........: ',f12.5)
7010  format(' Span .........[Angstroms]: ',f12.5)
7020  format(' Span ................[nm]: ',f12.5)
7030  format(' RMS thickness [Angstroms]: ',f12.5)
7040  format(' RMS thickness .......[nm]: ',f12.5)
7050  format(' Second moment [Angstroms]: ',f12.5)
7060  format(' Second moment .......[nm]: ',f12.5)

      write(6,7000)boundfrac
      write(1,*)' '
      write(1,7000)boundfrac
      if ((qisang.eq.'y').or.(qisang.eq.'Y').or.(qisang.eq.' ')) then
        write(6,7010)span
        write(6,7030)rms
        write(6,7050)sigma
        write(1,7010)span
        write(1,7030)rms
        write(1,7050)sigma
      else if ((qisang.eq.'n').or.(qisang.eq.'n')) then
        write(6,7020)span
        write(6,7040)rms
        write(6,7060)sigma
        write(1,7020)span
        write(1,7040)rms
        write(1,7060)sigma
      end if

*==============================================================================
*     WRITE DATA
*==============================================================================

*     DIFFERENCE BETWEEN INTERPOLATION ERROR ESTIMATE
      fname=ident//'ERR.TXT'
      write(6,1020)fname
      write(1,1020)fname
      call save_data(fname,n,x,e)

*     VOLUME FRACTION PROFILE
      fname=ident//'VFR.TXT'
      write(6,1020)fname
      write(1,1020)fname
      call save_data(fname,n,xx,yy)

      write(6,*)' '
      write(1,*)' '
      write(6,*)'END OF VOLUME FRACTION PROFILE ANALYSIS'
      write(1,*)'END OF VOLUME FRACTION PROFILE ANALYSIS'

      goto 9999

*==============================================================================
*     ERROR MESSAGES
*==============================================================================

8001  call showprompts(prompts(1))
      goto 9999

8002  call showprompts(prompts(2))
      goto 9999

9999  close(1)
      end


*==============================================================================
*     SUBROUTINES FOR PHASE RECOVERY
*==============================================================================


      subroutine rephase(dk,dx,rotx)
*     CARRIES OUT INVERSION FOR INTENSITY IN REAL PART OF YVEC
*     IN 0 PT. FORMAT (SHOULD BE SYMMETRICAL)
*     RETURNS DENSITY IN REAL PART OF YVEC
*     IN FORMAT SPECIFIED BY ROTX E.G. ROTX=-0.6,0.
*     REQUIRES INPUT OF DK,OUTPUTS DX
      common/fft/yvec,n0,m0
      real*8 :: yvec(2,4096)
      real*8 rst0,rotx,dk,dx

      rst0=yvec(1,n0/2+1)
c     write(6,*)'Calling tlz...'
      call tlz(1)
c     write(6,*)'Calling vlog...'
      call vlog
c     write(6,*)'Calling diff...'
      call diff
c     write(6,*)'Calling disp...'
      call disp
      rst0=0.5d0*dlog(rst0)
c     write(6,*)'Calling integ...'
      call integ(rst0)
c     write(6,*)'Calling vexp...'
      call vexp
c     write(6,*)'Calling tlfft0...'
      call tlfft0(-1.d0,dk,dx,0.0d0,rotx)
      return
      end


      subroutine diff
*     IN SITU DIFFERENTIATION(I.P.0 PT TO HALF PT)
*     ONLY NECESSARY TO DIFF. REAL PART AS IMAG PART ZEROED
      common/fft/yvec,n0,m0
      real*8 :: yvec(2,4096)
      real*8 yst

      yst=yvec(1,1)
      do 100 ll=1,n0-1
        yvec(1,ll)=yvec(1,ll+1)-yvec(1,ll)
100     continue
      yvec(1,n0)=yst-yvec(1,n0)
      return
      end


      subroutine integ(y0)
*     IN SITU INTEGRATION ROUTINE(I.P. HALF PT TO 0 PT)
*     INPUT Y0 AS VALUE AT ZERO
      common/fft/yvec,n0,m0
      real*8 yvec(2,4096)
      real*8 yst(2),yst1(2),y0

      n2=n0/2
      yst(1)=yvec(1,n2+1)
      yst(2)=yvec(2,n2+1)
*     SET INITIAL VALUE
      yvec(1,n2+1)=y0
      yvec(2,n2+1)=0.d0
*     INTEGRATE IN POSITIVE DIRECTION
*     TAKING CARE NOT TO OVERWRITE ARRAY
      do 100 ll=n2+2,n0
        do 101 ire=1,2
           yst1(ire)=yst(ire)
           yst(ire)=yvec(ire,ll)
           yvec(ire,ll)=yvec(ire,ll-1)+yst1(ire)
101     continue
100   continue
*     INTEGRATE IN NEGATIVE DIRECTION
      do 120 ll=n2,1,-1
        do 121 ire=1,2
           yvec(ire,ll)=yvec(ire,ll+1)-yvec(ire,ll)
121     continue
120   continue
      return
      end


      subroutine disp
*     DISPERSION INTEGRAL ROUTINE
*     HALF POINT I/O
      common/fft/yvec,n0,m0
      real*8 :: yvec(2,4096)
      real*8 risfl,dx1,dx2,rotx,rotk

      dx1=1.d0
*     DX1 ARBITRARY
      risfl=-1.d0
      rotx=-0.5d0
      rotk=rotx
*     INVERSE HALF-PT TRANSFORM
      call tlfft0(risfl,dx1,dx2,rotk,rotx)
*     MULTIPLY BY 2*HEAVISIDE STEP FUNCTION
      do 100 ll=1,n0/2
        yvec(1,ll)=0.0d0
        yvec(2,ll)=0.0d0
100   continue
      do 110 ll=n0/2+1,n0
        yvec(1,ll)=2.0d0*yvec(1,ll)
        yvec(2,ll)=2.0d0*yvec(2,ll)
110    continue
      risfl=-risfl
*     DIRECT HALF-PT. TRANSFORM
      call tlfft0(risfl,dx2,dx1,rotx,rotk)
      return 
      end


      subroutine tlz(iz)
*     TO ZERO REAL OR IMAG PART OF YVEC
*     IZ=0 FOR REAL,IZ=1 FOR IMAG
      common/fft/yvec,n0,m0
      real*8 :: yvec(2,4096)
      integer iz
      ire=iz+1
      do 100 ll=1,n0
        yvec(ire,ll)=0.d0
100   continue
      return
      end


      subroutine vlog
*     NEED ONLY OPERATE ON REAL PART
*     FACTOR OF 0.5 TO TAKE SQUARE ROOT OF INTENSITY
*     TO GIVE MODULUS OF FT
      common/fft/yvec,n0,m0
      real*8 :: yvec(2,4096)

      do 100 ll=1,n0
      yvec(1,ll)=0.5d0*dlog(yvec(1,ll))
100   continue
      return
      end


      subroutine vexp
*     COMPLEX EXPONENTIAL
      common/fft/yvec,n0,m0
      real*8 :: yvec(2,4096)
      real*8 yr,yi

      do 100 ll=1,n0
        yr=yvec(1,ll)
        yi=yvec(2,ll)
        yr=dexp(yr)
        yvec(1,ll)=yr*dcos(yi)
        yvec(2,ll)=yr*dsin(yi)
100   continue
      return
      end


      subroutine tlfft0(risfl,dx1,dx2,rot1,rot2)
      common/fft/y2,n0,m0
      real*8 :: y2(2,4096)
      real*8 rn0,risfl,dx1,dx2,rot1,rot2
      real*8 pi2,frac,gam

      pi2=6.28318530718002d0
      ialt=1
      m0=il2(n0)
      rn0=dble(float(n0))
      frac=-risfl*rot2/rn0
      call vphas0(frac,ialt)
      call tldft0(risfl)
      frac=-risfl*rot1/rn0
      call vphas0(frac,ialt)
      frac=risfl*((rot1*rot2)/rn0+0.5d0*(rot1+rot2))
      if(risfl.gt.0.d0)gam=dx1
      if(risfl.lt.0.d0)gam=dx1/pi2
      call vphc0(frac,gam)
      dx2=pi2/(rn0*dx1)
      return
      end


      subroutine tldft0(risfl)
*     HOME MADE FFT ROUTINE,USING COMPLEX ARITMETIC
*     INPUT COMPLEX YVEC(1-N0) ,N0,M0
*     WHERE M=2**N
      common/fft/y2,n0,m0
      real*8 :: y2(2,4096)
      real*8 :: yst1(2)
      real*8 :: ph(2)
      real*8 :: ph1(2)
      real*8 :: ph0(2)
      real*8 :: phs(2)
      real*8 :: dd(2)
      real*8 pi,one,zero,half,theta,risfl

      pi=3.1415926544d0
      one=1.0d0
      zero=0.0d0
      half=0.5d0
      call scr0
      ph0(1)=-one
      ph0(2)=zero
      ph1(1)=-one
      ph1(2)=zero
      theta=pi
      inc=1
      nbl2=n0
      do 100 nu=1,m0
        ph(1)=one
        ph(2)=zero
        nbl2=nbl2/2
        do 110 l1=1,inc
           lb1=l1
           lb2=lb1+inc
           do 120 l2=1,nbl2
              yst1(1)=y2(1,lb1)
              yst1(2)=y2(2,lb1)
              dd(1)=ph(1)*y2(1,lb2)-ph(2)*y2(2,lb2)
              dd(2)=ph(1)*y2(2,lb2)+ph(2)*y2(1,lb2)
              y2(1,lb1)=yst1(1)+dd(1)
              y2(2,lb1)=yst1(2)+dd(2)
              y2(1,lb2)=yst1(1)-dd(1)
              y2(2,lb2)=yst1(2)-dd(2)
              lb1=lb2+inc
              lb2=lb1+inc
120           continue
           phs(1)=ph(1)
           phs(2)=ph(2)
           ph(1)=ph0(1)*phs(1)-ph0(2)*phs(2)
           ph(2)=ph0(1)*phs(2)+ph0(2)*phs(1)
110        continue
        theta=half*theta
        ph1(1)=dcos(theta)
        ph1(2)=dsin(theta)
        ph0(1)=ph1(1)
        ph0(2)=risfl*ph1(2)
        inc=inc+inc
100     continue
      return
      end


      subroutine scr0
*     BITWISE SCRAMBLING OF YVEC -USED IN TLDFT0
      common/fft/y2,ns,ms
      real*8 :: y2(2,4096)
      real*8 :: yst(2)

      do 100 ll=1,ns
        j=ll-1
        ibitr=0
        do 200 i=1,ms
           j2=j/2
           ibitr=ibitr+ibitr+j-j2-j2
           j=j2
200        continue
        i=ibitr+1
        if(i.le.ll)goto 100
        yst(1)=y2(1,ll)
        yst(2)=y2(2,ll)
        y2(1,ll)=y2(1,i)
        y2(2,ll)=y2(2,i)
        y2(1,i)=yst(1)
        y2(2,i)=yst(2)
100     continue
      return
      end


      function il2(n0)
      integer n0
*     LOG2 OF N0
      il2=0
      n1=1
50    n1=n1+n1
      if (n1.gt.n0) return
      il2=il2+1
      goto 50
      end


      subroutine vphc0(frac,gam)
*     MULTIPLIES Y2 BY CONST GAM AND VPHASE(FRAC)
      common/fft/y2,ns,ms
      real*8 :: y2(2,4096)
      real*8 theta,pi2,frac,gam,deps
      real*8 :: ph(2)
      real*8 :: yst(2)

      deps=1.d-10
      pi2=6.28318530718002d0
      if (dabs(frac).lt.deps) goto 200
      theta=pi2*frac
      ph(1)=dcos(theta)
      ph(2)=dsin(theta)
      do 100 ll=1,ns
        yst(1)=y2(1,ll)
        yst(2)=y2(2,ll)
        y2(1,ll)=gam*(ph(1)*yst(1)-ph(2)*yst(2))
        y2(2,ll)=gam*(ph(1)*yst(2)+ph(2)*yst(1))
100     continue
      return
200     do 210 ll=1,ns
        y2(1,ll)=gam*y2(1,ll)
        y2(2,ll)=gam*y2(2,ll)
210     continue
      return
      end


      subroutine vphas0(frac,ialt)
*     MULTIPLIES Y2(LL) BY EXP(2*PI*FRAC*LL)
*     ALTERNATES SIGN IF IALT=1
      common/fft/y2,ns,ms
      real*8 :: y2(2,4096)
      real*8 theta,pi2,frac,deps,ph(2),yst(2),ph0(2),phst(2)
      real*8 :: ph(2)
      real*8 :: yst(2)
      real*8 :: ph0(2)
      real*8 :: phst(2)
      integer ialt

      deps=1.d-10
      if(dabs(frac).lt.deps)goto 200
      pi2=6.28318530718002
      theta=pi2*frac
      ph0(1)=dcos(theta)
      ph0(2)=dsin(theta)
      ph(1)=1.
      ph(2)=0.
      do 100 ll=2,ns
        yst(1)=y2(1,ll)
        yst(2)=y2(2,ll)
        phst(1)=ph(1)
        phst(2)=ph(2)
        ph(1)=ph0(1)*phst(1)-ph0(2)*phst(2)
        ph(2)=ph0(1)*phst(2)+ph0(2)*phst(1)
        y2(1,ll)=ph(1)*yst(1)-ph(2)*yst(2)
        y2(2,ll)=ph(1)*yst(2)+ph(2)*yst(1)
100     continue
200     if (.not.ialt.eq.1) return
      do 1010 ll=2,ns,2
        y2(1,ll)=-y2(1,ll)
        y2(2,ll)=-y2(2,ll)
1010     continue
      return
      end


      subroutine xft0(xx,ll,dx,rot,n0)
      common/ln/icln,k0
      real*8 dx,rot,xx,k0
      integer ll,n0

      xx=(dble(float(ll-n0/2-1))-rot)*dx
      if (icln.eq.1) xx=k0*dexp(xx)
      return
      end


*==============================================================================
*     SUBROUTINE FOR INTEGRATION
*==============================================================================


      subroutine trapiz(x_data_array,y_data_array,channels,start,sumz)
      real*4 :: x_data_array(4096)
      real*4 :: y_data_array(4096)
      real*4 sumz
      integer*4 channels,start

      h=0.0
      sumz=0.0

      h=(x_data_array(channels)-x_data_array(start))/(channels-start)

      do 1000 ll=start+1,channels-1
        sumz=sumz+2*y_data_array(ll)
1000  end do

      sumz=0.5*h*(y_data_array(start)+sumz+y_data_array(channels))

      return
      end


C=============================================================================
C==============================================================================
