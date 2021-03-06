*
* Fit tof spectrum to sum of N Voigt or Gaussian functions.
* Fitting parameters are heights and sigmas in y. Position fixed by mass.
* Individual widths or areas can also be fixed
      PROGRAM TFIT_FSE
      PARAMETER (NP=10,NDET=100)
      REAL CM(NP,NP),CMI(NP,NP),ERR(NP)
      REAL WORK(23000),F(1000),P(30),PS(30),PST(30)
      REAL TH(100),RL0(100),RL1(100),DT0A(100),
     $M,K0,K1,L0,SIGY(20)
      REAL COMP(20),DCOMP(20),WIDW(20),DWIDW(20)
    
      COMMON QA(1000),E0A(1000),YA(5,1000),TRES(5,1000),e1,NBACK
      COMMON M(20),DTG(20),XS(20),WID(5),NPAR
      COMMON XI(1000),YI(1000),EI(1000),V(30),DTE(20),TR0(20),NPEAKS

     
C  READ IN DATA USING GENIE ROUTINES.
      INCLUDE 'GENSOURCE:FUNCTCOM.CMN'
      CALL FUNCTION_IN

      WRITE(6,*) ' NUMBER OF INPUT DATA POINTS=',LPT
      YMAX=0.0
      DO I=1,LPT
       IF(YIN(I).GT.YMAX) YMAX=YIN(I)
      END DO
      DO I=1,LPT
       IF(YIN(I).EQ.0.0) THEN
        YIN(I)=YMAX*1E-6
        EIN(I)=1E6
       END IF
      END DO
      
      DO I=1,LPT
      XI(I)=XIN(I)
      YI(I)=YIN(I) 
      EI(I)=EIN(I)
      END DO

C  TAKE INPUT FROM COMMAND FILE
      N=ST_VAR(1) ! DETECTOR NUMBER
      ND=N
      E1=ST_VAR(2) ! ANALYSER ENERGY
      DE1=ST_VAR(3) ! Energy resolution
      IPNO=ST_VAR(4) ! IP File no.
      NPEAKS=ST_VAR(5) ! Number of peaks.
      MAXFUN=ST_VAR(6) ! Maximum no of iterations.
      NBACK=ST_VAR(7) ! Decides whether linear backgnd included.
      NPL=1
      IF(NPEAKS.GT.5) THEN
       WRITE(6,*) ' NUMBER OF PEAKS MUST BE .LE. 5'
       STOP
      END IF 
      PI=ACOS(-1.0)

C 2.  READ IN VALUES OF THETA AND L1 FROM FILE
      CALL PAR_READ(TH,DT0A,RL0,RL1,IPNO,NDET,NS)

C 3. DEFINE FIXED RESOLUTION COMPONENTS
      CALL RES_READ(DL0,DL1,DT0,DW)
      DL0=DL0/100.0
      DL1=DL1/100.0
      DW=DW/100.0
      DO I=1,NS
      TH(I)=TH(I)*PI/180.0
      END DO

      ISPEC=N
      L0=RL0(N)
      L1=RL1(N)
      THETA=TH(N)
      DTHETA=0.5*DW/RL1(N)
      THETAW=TH(N)*180.0/PI
      DTHW=DTHETA*180.0/PI
      T0=DT0A(N)

      WRITE(6,*) ' DETECTOR NUMBER=',N
      WRITE(6,*) ' L0=',L0,' +-',DL0,' METRES'
      WRITE(6,*) ' L1=',L1,' +-',DL1,' METRES'
      WRITE(6,*) ' THETA=',THETAW,' +-',DTHW,' DEGREES'      
      WRITE(6,*) ' E1=',E1,' +-',DE1,' meV'
      WRITE(6,*) ' T0=',T0,' +-',DT0, 'uSEC'
      IF(NBACK.NE.1) THEN
       WRITE(6,*) ' NO LINEAR BACKGROUND IN FIT'
      ELSE
       WRITE(6,*) ' LINEAR BACKGROUND INCLUDED'
      END IF
      WRITE(6,*) ' NUMBER OF PEAKS=',NPEAKS
*      WRITE(6,*) ' MAXFUN=',MAXFUN
      NXSP=0 ! Number of peak amplitude fitting parameters
      DO I=1,NPEAKS
      M(I)=ST_VAR(I+7)
      XS(I)=ST_VAR(I+11)
      IF(XS(I).EQ.0.0) NXSP=NXSP+1
      WID(I)=ST_VAR(I+15)
      IF(M(I).EQ.0.0) GO TO 130
      END DO
      IF(NXSP.EQ.NPEAKS) XS(1)=1.0 ! amp parameters not independent.

*  Calculate tof at at peak centre.
      DO I=1,NPEAKS
      TR0(I)=TREC(L0,L1,THETA,E1,T0,M(I))
      END DO
      
* 5. Calculate Q and E0 and y for each tof

       V1=SQRT(E1/5.2276E-6)
       T1=L1/V1 
      DO I=1,LPT
       T=XIN(I)
       V0=L0/(T*1E-6-T0*1E-6-T1)
       E0=5.2276E-6*V0**2
       K0=SQRT(E0/2.0717)
       K1=SQRT(E1/2.0717)
       QA(I)=SQRT(K0**2+K1**2-2.0*K0*K1*COS(THETA))
       E0A(I)=E0
      END DO

* Calculate y for each t and for each mass.
      DO IM=1,NPEAKS
      DO I=1,LPT
      W=E0A(I)-E1
      RM=M(IM)/1.00867 ! Mass in multiples of neutron mass.
      WR=2.0717*QA(I)**2/RM
      YA(IM,I)=0.2413*(RM/QA(I))*(W-WR)
      END DO
      END DO

      QAV=QA(LPT/2)

* 6. Calculate resolution function for each mass in time of flight.
      DO IM=1,NPEAKS
       DTG(IM)= TRESN(L0,DL0,L1,DL1,THETA,DTHETA,E1,0,T0,DT0,M(IM))
       DTE(IM)= TRESN(L0,0.0,L1,0.0,THETA,0.0,E1,DE1,T0,0.0,M(IM))
      END DO    
      DT=XIN(2)-XIN(1)
      DO IM=1,NPEAKS
       DO IT=1,LPT/2
        TV=(IT-1)*DT
        IF(E1.LT.6000) THEN ! Res fn is Voigt function.
         TRES(IM,IT)=DT*VOIGT(TV,DTG(IM),DTE(IM),0)
        ELSE ! Gaussian resolution.
         WID1=SQRT(DTG(IM)**2+DTE(IM)**2)
         TRES(IM,IT)=DT*VOIGT(TV,WID1,0,0)
        END IF
       END DO
       DO IT=LPT/2+1,LPT
        TV=(IT-1-LPT)*DT
        IF(E1.LT.6000) THEN ! Res fn is Voigt function.
         TRES(IM,IT)=DT*VOIGT(TV,DTG(IM),DTE(IM),0)
        ELSE ! Gaussian resolution.
         WID1=SQRT(DTG(IM)**2+DTE(IM)**2)
         TRES(IM,IT)=DT*VOIGT(TV,WID1,0,0)
        END IF
       END DO
      END DO


      WRITE(6,*) ' '
      WRITE(6,*) '   MASS     ','  POSITION ','    INTENSITY ',
     $'   WIDTH    ','    DTG     ','    DTE     '
      DO I=1,NPEAKS
      WRITE(6,50) M(I),TR0(I),XS(I),WID(I),DTG(I),DTE(I)
   50 FORMAT(' ',1P6E12.4)
      END DO

C Calculate estimates for amplitude.
      YMAX=0.0
      DO I=1,LPT
      IF(YMAX.LT.YI(I)) THEN
      YMAX=YI(I)
      IMAX=I
      END IF
      END DO

* Estimate of standard deviation in tof.
      DO I=1,NPEAKS
      WD=800/11.604 ! Estimated Debye energy in meV.
      SIGY(I)=SQRT(0.75*0.1196*M(I)*WD)
      SIGY(I)=SQRT(SIGY(I)**2)
      END DO

      NPAR=0 ! Number of fitting parameters other than scale and background.
      DO I=1,NPEAKS
       IF(XS(I).EQ.0.0) THEN 
        NPAR=NPAR+1
        P(NPAR)=0.1/SIGY(I) 
       END IF
       IF(WID(I).EQ.0.0) THEN 
        NPAR=NPAR+1
        P(NPAR)=SIGY(I) 
       END IF
      END DO
      N=NPAR+1
      P(N)=1
*      P(N)=YMAX ! Scale Factor.
      IF(NBACK.EQ.1) N=N+2
*      P(N-2)=1.0
* LINEAR BACKGROUND Y=P(N-1)*X+P(N)
      IF(NBACK.EQ.1) THEN
      X1=0.5*(XI(1)+XI(2))
      X2=0.5*(XI(LPT)+XI(LPT-1))
      Y1=0.5*(YI(1)+YI(2))
      Y2=0.5*(YI(LPT)+YI(LPT-1))
      P(N-1)=(Y2-Y1)/(X2-X1)
      P(N)=(X1*Y2-Y1*X2)/(X1-X2)
      END IF
      WRITE(6,*) ' Number of fitting parameters=',N
C SCALE PARAMETERS
      DO 70 I=1,N
      V(I)=1.0/P(I)
      PS(I)=P(I)*V(I)
      PST(I)=P(I)! Store start values of fit parameters.
  70  CONTINUE

C  CALCULATE DMAX USED IN VA05A.
      DMAX=3.0
C  H IS DISTANCE BETWEEN X PTS USED TO CALCULATE PARTIAL DERIVATIVES.
      H=DMAX*1E-2
C
      ACC=0.001
   80 IPRINT=0

      CALL VA05A(LPT,N,F,PS,H,DMAX,ACC,MAXFUN,IPRINT,WORK)

      IF(IPRINT.EQ.1.AND.NDMAX.NE.1) THEN 
*      TYPE *,' SUM OF SQUARES FAILED TO DECREASE'
*      TYPE *, ' DMAX REDUCED TO 0.1'
      NDMAX=1
      DMAX=0.1
      H=DMAX*1E-2
      GO TO 80
      END IF

      IF(IPRINT.EQ.2) TYPE *,' MAXIMUM NUMBER OF CALLS TO VAO5A MADE '

      CS=CHISQ(PS,N,LPT) ! Calculate chi-square
      CS=CS/(LPT-N) ! Reduced chi-square.
      WRITE(6,*) ' FIT COMPLETED CHI=',CS

      WRITE(6,*) 'ERRORS NOW BEING CALCULATED.'

* Calculate errors on fitted parameters.
      CALL ERRORS(N,LPT,PS,CM,CMI,ERR)

* Check calculated errors.
      CALL ERRCHECK(N,LPT,PS,ERR)

      DO 90 I=1,N
      P(I)=PS(I)/V(I)
      ERR(I)=ERR(I)/V(I)
   90 CONTINUE

      DO I=1,N
      WRITE(6,*) I, ' PST=',PST(I),' P=',P(I)
      END DO

      SUM=0.0 ! Total area of all peaks
      NP1=0
      DO I=1,NPEAKS
       IF(XS(I).EQ.0.0) THEN
        NP1=NP1+1
        SUM=SUM+ABS(P(NP1))
       ELSE
        SUM=SUM+ABS(XS(I)*P(NPAR+1))
       END IF
       IF(WID(I).EQ.0.0) NP1=NP1+1
      END DO

      NP1=0
      DO I=1,NPEAKS
       IF(XS(I).EQ.0.0) THEN
        NP1=NP1+1
        COMP(I)=ABS(P(NP1))/SUM
        DCOMP(I)=ABS(ERR(NP1))/SUM
       ELSE
        COMP(I)=ABS(XS(I)*P(NPAR+1)/SUM)
        DCOMP(I)=ABS(XS(I)*ERR(NPAR+1)/SUM)
       END IF
       IF(WID(I).EQ.0.0) THEN
        NP1=NP1+1
        WIDW(I)=P(NP1)
        DWIDW(I)=ERR(NP1)
       ELSE
        WIDW(I)=WID(I)
        DWIDW(I)=1E-6
       END IF
       WIDW(I)=ABS(WIDW(I)) ! Eliminate negative values
       DWIDW(I)=ABS(DWIDW(I)) ! Eliminate negative values
       
      END DO

      WRITE(6,*) ' '
      write(6,*) ' Relative wt is proportional to bound x-sect X',
     $' concentration'


      WRITE(6,*) ' PEAK NO','    MASS    ',' RELATIVE WT ',
     $'   ERROR    ','  SD OF J(Y) ','   ERROR'
      OPEN(UNIT=2,FILE='TFIT.OUT',STATUS='NEW')
      WRITE(2,*) ND
      DO I=1,NPEAKS
      WRITE(2,95) I,M(I),COMP(I),DCOMP(I),WIDW(I),DWIDW(I)
      WRITE(6,95) I,M(I),COMP(I),DCOMP(I),WIDW(I),DWIDW(I)
   95 FORMAT('      ',I3,1P5E12.4)
      END DO
      CLOSE(2)
      
      CHI=0.0
      CALL CTS(YOUT,LPT,NPEAKS,P)

  130 CONTINUE  

      CALL FUNCTION_OUT

      STOP
      END

      SUBROUTINE CALFUN(M,N,F,PS)
      REAL F(1000),PS(30),P(30)
      COMMON QA(1000),E0A(1000),YA(5,1000),TRES(5,1000),e1,NBACK
      COMMON AM(20),DYG(20),XS(20),WID(5),NPAR
      COMMON XI(1000),YI(1000),EI(1000),V(30),DYE(20),TR0(20),NPEAKS

      DO 10 I=1,N
      P(I)=PS(I)/V(I)
   10 CONTINUE      
      
      CHI=0.0
      CALL CTS(F,M,NPEAKS,P)
      DO 20 I=1,M
      X=XI(I)
      F(I)=(F(I)-YI(I))/EI(I)
      CHI=CHI+F(I)**2
   20 CONTINUE
      CHI=CHI/M
*      write(6,*) ' chi=',chi
   30 FORMAT(7E10.3)

      RETURN
      END

C Function to calculate tof at recoil peak in micro seconds
C L0 = incident flight path in metres.
C L1= final flight path in metres.
C TH = scattering angle in radians.
C E1 = analyser energy in meV.
C DT0 = time delay in microseconds.
C M=atomic mass in amu.
      FUNCTION TREC(L0,L1,TH,E1,DT0,M)
      REAL L0,L1,M
      M=M/1.00867 ! Convert to multiple of neutron mass.
      ARG=M**2-SIN(TH)**2
      IF(ARG.LT.0.0) ARG=0.0
      FACT=(COS(TH)+SQRT(ARG))/(M+1)
      M=M*1.00867
      V1=SQRT(E1/5.2276E-6)
      V0=V1/FACT
      TREC=L0/V0+L1/V1+DT0*1E-6
      TREC=TREC*1E6
      RETURN
      END
 
C Function to calculate y.
C L0 = incident flight path in metres.
C L1= final flight path in metres.
C TH = scattering angle in radians.
C T = time of flight in microsec.
C E1 = analyser energy in meV.
C DT0 = time delay in microseconds.
C M=atomic mass in amu.
       FUNCTION Y(L0,L1,TH,E1,T,DT0,M)
       REAL L0,L1,M,K0,K1
       M=M/1.00867 ! Convert to multiple of neutron mass.
       V1=SQRT(E1/5.2276E-6)
       T1=L1/V1 
       T0=T*1E-6-DT0*1E-6-T1
       V0=L0/T0
       E0=5.2276E-6*V0**2
       W=E0-E1
       K0=SQRT(E0/2.0717)
       K1=SQRT(E1/2.0717)
       Q=SQRT(K0**2+K1**2-2.0*K0*K1*COS(TH))
       WR=2.0717*Q**2/M
       Y=0.2413*M*(W-WR)/Q
       M=M*1.00867
       RETURN
       END

C
C    Voigt function centred at X=X0
C    sigma is gaussian standard deviation
c    DYE is Lorentzian DYE.
C    Peak area normalised to 1.
C    WRITTEN BY WIFD
C    Modified by JM
C
      FUNCTION VOIGT(X,SIGMA,DYE,X0)

      DOUBLE PRECISION WR,WI,XX,YY
C
      
      GAMMA=DYE*2.0
      XS=X
      X=XS-X0
      OVRTPI=0.564189584
      OVRT2=0.707106781
      BTEM=OVRT2/SIGMA
      ATEM=OVRTPI*BTEM
      XTEM=X*BTEM
      YTEM=0.5*GAMMA*BTEM
      XX= DBLE(XTEM)
      YY= DBLE(YTEM)
      CALL WERF(WR,WI,XX,YY)
      SWR=SNGL(WR)
      SWI=SNGL(WI)
      CTEM=ATEM*BTEM
      VOIGT=ATEM*SWR
      X=XS
      DWRDX=-2.*(XTEM*SWR-YTEM*SWI)
      DWRDY= 2.*(YTEM*SWR+XTEM*SWI-OVRTPI)
      DERX=CTEM*DWRDX
      DERS=-ATEM*(SWR+DWRDX*XTEM+DWRDY*YTEM)/SIGMA
      DERG=0.5*CTEM*DWRDY
C
      RETURN
      END
c
c
	subroutine WERF(rs1,rs2,xx,yy)
c	W.I.F.David 25-May-84
	implicit real*8		(a-h,o-z)
	real*8			lambda
	logical			b
	x=dabs(xx)
	y=dabs(yy)
	if (y .lt. 4.29 .and. x .lt. 5.33) go to 1
	h= 0.
	nc= 0
	nu= 8
	lambda= 0.
	b= .true.
	go to 2
 1	s=(1.0-y/4.29)*dsqrt(1.0-x**2/28.41)
	h=1.6*s
	h2=2.0*h
	nc=6+idint(23.0*s)
	nu=9+idint(21.0*s)
	lambda=h2**nc
	b= .false.
	if (lambda .eq. 0.) b= .true.
 2	r1=0.
	r2=0.
	s1=0.
	s2=0.
	n=nu+1
 3	n=n-1
	fn=n+1
	t1=y+h+fn*r1
	t2=x-fn*r2
	c=0.5/(t1**2+t2**2)
	r1=c*t1
	r2=c*t2
	if (h .le. 0.0 .or. n .gt. nc) go to 4
	t1= lambda+s1
	s1=r1*t1-r2*s2
	s2=r2*t1+r1*s2
	lambda=lambda/h2
 4	if (n .gt. 0) go to 3
	if (b) go to 6
	rs1=s1
	rs2=s2
	go to 7
 6	rs1=r1
	rs2=r2
 7	rs1= 1.12837916709551*rs1
	if (y .eq. 0.0) rs1= dexp(-x**2)
	rs2= 1.12837916709551*rs2
	if (xx .lt. 0) rs2= -rs2
	return
	end


      FUNCTION CHISQ(PS,N,M)
      
      COMMON QA(1000),E0A(1000),YA(5,1000),TRES(5,1000),e1,NBACK
      COMMON AM(20),DYG(20),XS(20),WID(5),NPAR
      COMMON XI(1000),YI(1000),EI(1000),V(30),DYE(20),TR0(20),NPEAKS
      REAL PS(30),P(30),F(1000)

      
      DO I=1,N
      P(I)=PS(I)/V(I)
      END DO

      CHISQ=0.0
      CALL CTS(F,M,NPEAKS,P)
      DO 20 I=1,M
      X=XI(I)
      F(I)=(F(I)-YI(I))/EI(I)
      CHISQ=CHISQ+F(I)**2
   20 CONTINUE

      END 


      SUBROUTINE ERRCHECK(N,M,PS,ERR)
      COMMON XI(1000),YI(1000),EI(1000),V(30)
      REAL PS(N),ERR(N)
     
      DO I=1,N
      CS=CHISQ(PS,N,M)
      PS(I)=PS(I)+ERR(I)
      CSP=CHISQ(PS,N,M)
      PS(I)=PS(I)-ERR(I)
      DCS=CSP-CS
      IF(DCS.LE.0.0) THEN
      ERR(I)=PS(I)
      WRITE(6,*) ' FOR PARAMETER',I, 
     $' CHISQ DECREASES AWAY FROM MINIMUM'
      ELSE IF(DCS.LT.1.0) THEN
      ERR(I)=ERR(I)/DCS
      WRITE(6,*) ' PROBLEM WITH ERROR ON PARAMETER',I
      WRITE(6,*) ' CHANGE IN CHISQ AT P(I)+DP(I)=',DCS
      END IF
      END DO

      END

      SUBROUTINE ERRORS(N,M,PS,CM,CMI,ERR)

      REAL PS(N),ERR(N),CM(N,N),PD(30),CMI(N,N),DP(30)
      COMMON XI(1000),YI(1000),EI(1000),V(30)

      DO I=1,N
      PD(I)=PS(I)
      END DO

* Find increments used in calculation of curvature matrix.
      CALL INCREMENT(PS,N,M,DP)

*  Calculate curvature matrix.
      CALL CURVATURE(PS,DP,N,M,CM)

*  Invert curvature matrix
      CALL MINV(N,CM,CMI,SM)
     
*  CALCULATE ERRORS
      DO I=1,N
      IF(SM.EQ.0.0) ERR(I)=SQRT(ABS(CMI(I,I)))
      IF(SM.EQ.1) ERR(I)=PS(I) ! Singular error matrix.
      END DO

      END

*  Subroutine to invert matrix. Uses NAG subroutine F01AAF.
*  A is N x N matrix. AI is inverse of A.
      SUBROUTINE MINV(N,A,AI,SM)
      REAL A(N,N),AI(N,N),W(50)

* Check for singularity.
      SM=0.0
      DO I=1,N
      ZERO=0.0
      DO J=1,N
      IF(A(I,J).NE.0.0) ZERO=1
      END DO
      IF(ZERO.EQ.0.0) THEN
      WRITE(6,*) ' ALL ELEMENTS IN ROW',I,
     $'OF ERROR MATRIX ARE ZERO'
      SM=1.0
      RETURN
      END IF
      END DO

      IFAIL=0
      CALL F01AAE(A,N,N,AI,N,W,IFAIL)
      IF(IFAIL.EQ.1) THEN 
      WRITE(6,*) ' ERROR MATRIX IS SINGULAR'
      DO I=1,N
      DO J=1,N
      AI(I,J)=0.0
      END DO
      END DO
      RETURN
      END IF

      END

* Calculate increments for calculation of curvature
* take increment which increases chisq by BETWEEN 1 AND 2.
      SUBROUTINE INCREMENT(PS,N,M,DP)
      REAL PS(N),DP(N),PD(30)
      COMMON XI(1000),YI(1000),EI(1000),V(30)

      DO I=1,N
      PD(I)=PS(I)
      END DO

      CS=CHISQ(PS,N,M)
      DO I=1,N
*  Find value of DCS >2.
      DP(I)=PS(I) ! Start value for increment
      DPMIN =0
      NINC=0 
   10 PD(I)=PS(I)+DP(I)

      NINC=NINC+1 ! Check for no minimum.
      IF(NINC.GT.10) THEN
      DP(I)=PS(I)
      GO TO 30
      END IF 
     
      CSP=CHISQ(PD,N,M)
      DCS=CSP-CS ! Calculate change in chisq in step dp from min.
      IF(DCS.LT.1.0) THEN
       DPMIN=DP(I)
       DP(I)=DP(I)*2.0
       GO TO 10
      ELSE IF(DCS.GT.2.0) THEN 
       DPMAX=DP(I)
      ELSE
       GO TO 30
      END IF

* Now find value of dp(i) which changes chisq by between 1 and 2.
   20 DP(I)=(DPMAX+DPMIN)/2.0
* Calculate change in chi-sq corresponding to dp(i).
      PD(I)=PS(I)+DP(I)
      CSP=CHISQ(PD,N,M)
      DCS=CSP-CS
      IF(DCS.GT.2.0) THEN
      DPMAX=DP(I)
      GO TO 20
      ELSE IF(DCS.LT.1.0) THEN
      DPMIN=DP(I)
      GO TO 20
      END IF
   30 CONTINUE
      PD(I)=PS(I) ! Reset pd to value at minimum.

      END DO
      END


      SUBROUTINE CURVATURE(PS,DP,N,M,CM)
      REAL PS(N),DP(N),CM(N,N),PD(30)
      COMMON XI(1000),YI(1000),EI(1000),V(30)

      DO I=1,N
      PD(I)=PS(I) ! Initialise pd.
      END DO

      C=CHISQ(PS,N,M) ! chisq at minimum.

      DO I=1,N
      DO J=1,N
      IF(J.LT.I) THEN
      CM(I,J)=CM(J,I)
      TIJ=0
      GO TO 5
      END IF

      DPI=DP(I)
      DPJ=DP(J)

      IF(I.NE.J) THEN
      PD(I)=PS(I)+DPI
      PD(J)=PS(J)+DPJ
      CIJ=CHISQ(PD,N,M)
      PD(I)=PS(I)
      PD(J)=PS(J)+DPJ
      CJ=CHISQ(PD,N,M)
      PD(I)=PS(I)+DPI
      PD(J)=PS(J)
      CI=CHISQ(PD,N,M)
      PD(I)=PS(I)
      PD(J)=PS(J)
      CM(I,J)=0.5*(CIJ-CI-CJ+C)/(DPI*DPJ)
      TIJ=(CIJ-C)/(DPI*DPJ)

      ELSE IF(I.EQ.J) THEN
      PD(I)=PS(I)+2.0*DPI
      CII=CHISQ(PD,N,M)
      PD(I)=PS(I)+DPI
      CI=CHISQ(PD,N,M)
      CM(I,I)=0.5*(CII-2.0*CI+C)/(DPI**2)
      TIJ=(CII-C)/DPI**2
      END IF
    5 CONTINUE    
  
      END DO
      END DO

      END

* Reads in resolution parameters.
      SUBROUTINE RES_READ(DL0,DL1,DT0,DW)
      OPEN(UNIT=3,FILE='EVS$DISK:[EVSMGR.CALIB.DATA]RESOLUTION.DAT',
     $STATUS='OLD',READONLY,SHARED)
      READ(3,*) DL0 ! Uncertainty in incident flight path (cm)
      READ(3,*) DL1 ! Uncertainty in scattered flight path(cm)
      READ(3,*) DT0 ! Tof uncertainty (usec)
      READ(3,*) DW  ! Detector width (cm)
      CLOSE(3)
      END
* IPNO is IP run number,NDET is maximum no of detectors,NS is
* number of detectors in IPNRUN.dat.
      SUBROUTINE PAR_READ(TH,DT0,L0,L1,IPNO,NDET,NS)
      REAL TH(NDET),L0(NDET),L1(NDET),DT0(NDET)
      CHARACTER RUN*4,FIN*40

* Define file name
      WRITE(RUN(1:4),'(I4.4)')IPNO
      FIN='EVS$disk0:[EVSMGR.CALIB.PAR]IP'//RUN
      WRITE(6,'(40A)') ' Instrument parameters read from file',FIN

      OPEN(UNIT=3,FILE=FIN,STATUS='OLD',READONLY,SHARED)
      NS=1
   10 READ(3,*,END=20) I,TH(NS),DT0(NS),L0(NS),L1(NS)
      NS=NS+1
      GO TO 10
   20 CLOSE(3)
      NS=NS-1

      END


* Calculate time of flight spectrum.
* NM=number of masses, M=mass, XS=amplitude, SIGY=width in y.
* TRES=resolution fn in tof, L0,L1, lengths in metres, TH angle in rad.
* T0 time delay in usec, E1 analyser energy in meV.
      SUBROUTINE CTS(CT,NPTS,NPEAKS,P)
      COMMON QA(1000),E0A(1000),YA(5,1000),TRES(5,1000),e1,NBACK
      COMMON M(20),DTG(20),XS(20),WID(5),NPAR
      REAL CT(NPTS),M
      REAL RES(512),CTM(512),P(30),AMP(10),WIDTH(10)

      NP=0
      DO I=1,NPEAKS
       IF(XS(I).EQ.0.0) THEN
        NP=NP+1
        AMP(I)=ABS(P(NP))
       ELSE
        AMP(I)=ABS(P(NPAR+1)*XS(I))
       END IF
       IF(WID(I).EQ.0.0)THEN
        NP=NP+1
        WIDTH(I)=P(NP)
       ELSE 
        WIDTH(I)=WID(I)
       END IF
      END DO


      DO IT=1,NPTS
       CT(IT)=0.0
      END DO

      DO IM=1,NPEAKS
       IF(M(IM).LT.3) THEN ! Harmonic Oscillator
        DSQV=12.0*4.18036*WIDTH(IM)**4/M(IM)
       ELSE ! Debye Approximation
        DSQV=12.8*4.18036*WIDTH(IM)**4/M(IM)
       END IF 
       DO IT=1,NPTS
        YV=YA(IM,IT)
        E0=E0A(IT)
        Q=QA(IT) 
        CTM(IT)=M(IM)*AMP(IM)*CP(YV,Q,M(IM),WIDTH(IM),DSQV)*E0**0.1/Q
        RES(IT)=TRES(IM,IT)
       END DO

       CALL C06EKE(1,CTM,RES,NPTS,IFAIL) ! Convolve with resn fn.
       DO IT=1,NPTS
        CT(IT)=CT(IT)+CTM(IT)
       END DO
      END DO

      END




C Function to calculate resolution in time.
C L0 +- DL0 = incident flight path in metres.
C L1 +- DL1 = final flight path in metres.
C TH +- DTH = scattering angle in radians.
C E1 +- DE1 = analyser energy in meV.
C T0 +- DT0 = time delay in microseconds.
C M=atomic mass in amu.
* Calls TREC.
      FUNCTION TRESN(L0,DL0,L1,DL1,TH,DTH,E1,DE1,T0,DT0,M)
      REAL L0,L1,M

      TR=TREC(L0,L1,TH,E1,T0,M)
*      write(6,*) ' tr=',tr,' t0=',t0
      DTL0=TREC(L0+DL0,L1,TH,E1,T0,M)-TR
      DTL1=TREC(L0,L1+DL1,TH,E1,T0,M)-TR
      DTTH=TREC(L0,L1,TH+DTH,E1,T0,M)-TR
      DTE1=TREC(L0,L1,TH,E1+DE1,T0,M)-TR
      DTT0=TREC(L0,L1,TH,E1,T0+DT0,M)-TR
      TRESN=SQRT(DTL0**2+DTL1**2+DTTH**2+DTE1**2+DTT0**2)
*      WRITE(6,*) ' L0=',L0,' +-',DL0,' L1=',L1,' +-',DL1
*      WRITE(6,*) ' TH=',TH,' +-',DTH,' E1=',E1,' =-',DE1
*      WRITE(6,*) ' T0=',T0,' +-',DT0, ' M=',M
*      WRITE(6,*) ' TRESN=',TRESN

      END   
      
* Calculates J(y) including first Sears Correction terms
      FUNCTION CP(Y,Q,M,SIG1,D2V)
      REAL Y,M,J0,J3

      SIG=ABS(SIG1)
      X=Y/SIG
      J0=EXP(-X**2/2)/(2.506628*SIG) ! 2.5066=SQRT(2PI)
      D3J=X*(3-X**2)*J0/SIG**3
      J3=-M*D2V*D3J/(150.49*Q)
      CP=J0+J3

      END


