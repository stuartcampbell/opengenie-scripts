CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ------------------------------------------------------------
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    N        I*4    I       -      Number of data.
C    X        R*4    I       N      X-coordinate of data.
C    Y        R*4    I       N      Y-coordinate of data.
C    E        R*4    I       N      Error-bar for data.
C    XW,YW,EW R*4    -       N      Work arrays.
C    PARS     R*4    I       5      Parameter-values that are to be 
C                                   fixed or used for the initial guess.
C    PARS     R*4    O       5      Best estimates of the parameters.
C    SIGPAR   R*4    O       5      1-sigma error-bars for PARS; can be
C                                   negative if IGOOD < 0!
C    IPROPT   I*4    I       5      Option for treatment of parameters:
C                                     0 = no initial estimate;
C                                     1 = use as "      "    ;
C                                     2 = fix at "      "    .
C    IGOOD    I*4    O       -      Negative if a problem occurs:
C                                     < -1 = serious; -1 = may be OK.
C
C Globals
C    PKCOM1
C    PKCOM2
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PKINIT     Gives crude non-parametric estimate of peak properties.
C     CHPOLY     Calculates "chi-squared" misfit for peak.
C
C History
C   Modifed by F. A. Akeroyd from original by D. S. Sivia
C-----------------------------------------------------------------------
C
C INPARS = max number of parameters for a peak
C
      SUBROUTINE PKMULTI(X,Y,E,N,XW,YW,EW,PARS,SIGPAR,
     1                   IPROPT,TYPES,NPEAKS,INPARS,IGOOD)
      EXTERNAL        CHMULTI
      REAL            CHMULTI
      INTEGER 	      NPEAKS, INPARS, K, NTYPE
      REAL            X(N),Y(N),E(N),XW(N),YW(N),EW(N)
      REAL            PARS(INPARS,NPEAKS), SIGPAR(INPARS,NPEAKS)
      INTEGER         IPROPT(INPARS,NPEAKS), TYPES(NPEAKS)
      INCLUDE 'pkdefs.inc'
      INTEGER ITYPES(MAXPEAKS)
      COMMON /PKMULTBLK/ ITYPES,NTYPE
      DATA            SMALL /1.0E-20/
      SAVE /PKMULTBLK/
C
      NPKPARS=INPARS
      IGOOD=-10
      NPARS=0
      NTYPE=NPEAKS
      IF (NPKPARS*NPEAKS .GT. MAXPKPARS) THEN
          WRITE(6,*) 'PEAK: Too many parameters - max allowed is ', MAXPKPARS
          RETURN
      ENDIF
      IF (NPEAKS .GT. MAXPEAKS) THEN
          WRITE(6,*) 'PEAK: Too many peaks - max allowed is ', MAXPEAKS
          RETURN
      ENDIF
      write(6,*) 'inpars,peaks', NPKPARS, NPEAKS
      DO I=1,NTYPE
          ITYPES(I)=TYPES(I)
      ENDDO
C *** we pass a fixed size array, but not all peaks have the same number of variable parameters
C *** this we set to "2" any we do not need
      DO J=1,NPEAKS
C GAUSS or LORNZ
          IF ((TYPES(J) .EQ. 1) .OR. (TYPES(J) .EQ. 3)) THEN 
              K=5
C VEXP
          ELSEIF (TYPES(J) .EQ. 6) THEN
              K=7
C POLY
          ELSEIF (TYPES(J) .EQ. 7) THEN
              K=NPKPARS
          ELSE
              K=6
          ENDIF
          DO I=K+1,NPKPARS
              IPROPT(I,J) = 2
              PARS(I,J) = 1.0
          ENDDO
      ENDDO
      DO J=1,NPEAKS
        DO I=1,NPKPARS
          K = I + (J-1) * NPKPARS
          IF (IPROPT(I,J) .EQ. 2) THEN
            PKPARS(K)=PARS(I,J)
            SIGPAR(I,J)=0.0
          ELSE
            NPARS=NPARS+1
            IPARS(NPARS)=K
            IF (IPROPT(I,J) .EQ. 1) PKPARS(K)=PARS(I,J)
          ENDIF
        ENDDO
      ENDDO
      write(6,*) 'fitting ',npars
      IF (NPARS.EQ.0 .OR. N.LE.NPARS) THEN
          WRITE(6,*) 'PEAK: Invalid number of parameters: ', npars
          RETURN
      ENDIF
      YMAX = Y(1)
      YMIN = Y(1)
      XAV=0.0
      DO I=1,N
          YMAX=MAX(YMAX,Y(I))
          YMIN=MIN(YMIN,Y(I))
          XAV=XAV+X(I)
      ENDDO
      XAV=XAV/FLOAT(N)
      DO I=1,N
          XW(I) = X(I)
          YW(I) = Y(I)
          IF (E(I) .GT. SMALL) THEN
              EW(I) = 1.0 / E(I)
          ELSE
              EW(I) = 20.0 / YMAX
          ENDIF
      ENDDO
C Estimate PARS
      CALL MULTIGUESS(X, Y, E, N, PARS, IPROPT, TYPES, NPEAKS, NPKPARS, 
     +               YMIN, YMAX, XAV, IGOOD)
      IF (IGOOD.LT.-2) THEN
        CALL VCOPY(PKPARS,PARS,NPKPARS*NPEAKS)
        RETURN
      ENDIF
      DO I=1,NPARS
        J=IPARS(I)
        IF (IPROPT(J,1).EQ.0) PKPARS(J)=PARS(J,1)
      ENDDO
      DO I=1,NPKPARS*NPEAKS
          DPAR(I) = MAX(ABS(PKPARS(I)),1.0)
      ENDDO
      DO I=1,NPARS
        J=IPARS(I)
        PARS(I,1)=PKPARS(J)
        DPAR(I)=DPAR(J)
      ENDDO
      write(6,*) 'pkpars',(PKPARS(K),k=1,NPEAKS*NPKPARS)
      write(6,*) 'pars',(PARS(K,1),k=1,NPARS)
      IGOOD=ISMPLX(PARS,DPAR,COVAR,NPARS,.FALSE.,CHMULTI,XW,YW,EW,N)
      DO I=1,NPARS
        J=IPARS(I)
        PKPARS(J)=PARS(I,1)
        SIG2=COVAR((I-1)*NPARS+I)
        SIGPAR(J,1)=SQRT(ABS(SIG2)+SMALL)
        IF (SIG2.LT.SMALL) THEN
          SIGPAR(J,1)=-SIGPAR(J,1)
          IF (IGOOD .GE. 0) IGOOD=-1
        ENDIF
      ENDDO
      DO I=1,NPKPARS*NPEAKS
          PARS(I,1)=PKPARS(I)
      ENDDO
      END
C
      REAL FUNCTION CHMULTI(PARS,NPARS,X,D,EW,N,IOFF)
C     ------------------------------------
C
      REAL      PARS(NPARS),X(N),D(N),EW(N)
      INCLUDE 'pkdefs.inc'
      INTEGER ITYPES(MAXPEAKS)
      REAL WORK1(10000),WORK2(10000)
      COMMON /PKMULTBLK/ ITYPES,NTYPE
      SAVE /PKMULTBLK/
      CHMULTI=1.0E10
      DO I=1,NPARS
          PKPARS(IPARS(I)) = PARS(I)
      ENDDO
      J=IOFF
      DO K=1,N
          WORK2(K) = 0.0
      ENDDO
      DO I=1,NTYPE
          IF (ITYPES(I) .EQ. 1) THEN
              CALL FNGAUS(X,WORK1,N,PKPARS(1+J))
          ELSE IF (ITYPES(I) .EQ. 2) THEN
              CALL FNGEXP(X,WORK1,N,PKPARS(1+J))
          ELSE IF (ITYPES(I) .EQ. 3) THEN
              CALL FNLORZ(X,WORK1,N,PKPARS(1+J))
          ELSE IF (ITYPES(I) .EQ. 4) THEN
              CALL FNLEXP(X,WORK1,N,PKPARS(1+J))
          ELSE IF (ITYPES(I) .EQ. 5) THEN
              CALL FNVOGT(X,WORK1,N,PKPARS(1+J))
          ELSE IF (ITYPES(I) .EQ. 6) THEN
              CALL FNVEXP(X,WORK1,N,PKPARS(1+J))
          ELSE IF (ITYPES(I) .EQ. 7) THEN
              CALL FNPOLY(X,WORK1,N,PKPARS(1+J),NPKPARS)
          ENDIF
          DO K=1,N
              WORK2(K) = WORK2(K) + WORK1(K)
          ENDDO
          J = J + NPKPARS
      ENDDO
      CHI=0.0
      DO K=1,N
         CHI=CHI+((WORK2(K)-D(K))*EW(K))**2 
      ENDDO
      CHMULTI=CHI
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MULTIGUESS(X, Y, E, N, PARS, IPROPT, TYPES, NPEAKS, 
     1     NPKPARS, YMIN, YMAX, XAV, IGOOD)
      INTEGER N, NPEAKS, NPKPARS, IGOOD
      REAL            PARS(NPKPARS,NPEAKS)
      INTEGER         IPROPT(NPKPARS,NPEAKS), TYPES(NPEAKS)
      REAL            X(N),Y(N),E(N)
      REAL    M, C
C *** must be given guess at peak centre
      IGOOD = 0
      YDIFF = YMAX - YMIN
C *** Guess background by drawing straight line
      M = (Y(N)-Y(1)) / (X(N)-X(1))
      C = Y(1) - M * X(1)
      DO J=1,NPEAKS
          XI=1.0
          DO I=1,NPKPARS
	      write(6,*) IPROPT(I,J),XI
              IF (IPROPT(I,J) .EQ. 0) THEN
                  IF (TYPES(J) .EQ. 7) THEN
                      PARS(I,J) = YDIFF / XI
                  ELSE
                      PARS(I,J) = 1.0
                  ENDIF
              ENDIF
              XI=XI*XAV
          ENDDO
          IF (TYPES(J) .NE. 7) THEN
              IF (IPROPT(1,J) .EQ. 0) PARS(1,J) = M
              IF (IPROPT(2,J) .EQ. 0) PARS(2,J) = C
              IF (IPROPT(4,J) .EQ. 0) THEN
                  WRITE(6,*) 'MULTIPEAK: No peak centre given'
              ENDIF
          ENDIF
      ENDDO
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
