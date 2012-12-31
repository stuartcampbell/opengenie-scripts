      SUBROUTINE SMOOGEN (pars_get, pars_put)

      INCLUDE 'genie_modules.inc'
	include 'genie_sources:transcom.cmn'
	EXTERNAL pars_get, pars_put

      real c(10000)
      integer k,wb,nww

c      CALL TRANSFORM_IN

c     dal file di comando

 	call module_get_int(pars_get,'order',m)
      call module_get_real(pars_get,'perc',ww)
      ns = 1

c      NS=ST_VAR_IN(1)
c      ww=ST_VAR_IN(2) 
c       M=ST_VAR_IN(3)
      icon=LPTIN+1
      ww=abs(ww)
      ww=icon*ww
      nww=int(ww/100)
      if(nww.lt.1)nww=1
      if(mod(nww,2).eq.0)nww=nww+1


c  INFORMAZIONI

      WRITE(6,*) ' No OF DATA POINTS=',LPTIN+1
      WRITE(6,*) ' No OF WINDOW POINTS=',nww
      WRITE(6,*) ' POL. ORDER=',M

c      calcolo coefficienti

       ld=0
       nr=(nww-1)/2
       nl=nr
       call savgol(c,icon,nl,nr,ld,m)
       

c      smoothing

       do j=1,icon
          yout(j)=0.
          do i=-nl,nr
          k=i+j
          if(k.lt.1)k=icon+k
          if(k.gt.icon)k=k-icon
          if(i.eq.0)wb=1
          if(i.lt.0)wb=-i+1
          if(i.gt.0)wb=icon+1-i
          yout(j)=yout(j)+yin(k)*c(wb)
          enddo
       enddo

c     estremi in media mobile

       do j=1,nl
       yout(j)=0.
       net=0
       do i=-nl,nr
       k=i+j
       if(k.ge.1)then
          yout(j)=yout(j)+yin(k)
          net=net+1
       endif
       enddo
       yout(j)=yout(j)/net
       enddo

       do j=icon-nr,icon
       yout(j)=0.
       net=0
       do i=-nl,nr
       k=i+j
       if(k.le.icon)then
         yout(j)=yout(j)+yin(k)
         net=net+1
       endif
       enddo
       yout(j)=yout(j)/net
       enddo
        


c      chiusura

       do i=1,icon
        xout(i)=xin(i)
        eout(i)=ein(i)
       enddo 
      

       XCAPTOUT=XCAPTIN
       YCAPTOUT=YCAPTIN

	call module_put_real_array(pars_put,'y',yout,icon)

c       CALL TRANSFORM_OUT

       STOP
       END






      SUBROUTINE savgol(c,np,nl,nr,ld,m)
      INTEGER ld,m,nl,np,nr,MMAX
      REAL c(np)
      PARAMETER (MMAX=6)
      INTEGER imj,ipj,j,k,kk,mm,indx(MMAX+1)
      REAL d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
      if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m
     *.or.m.gt.MMAX.or.nl+nr.lt.m) pause 'bad args in savgol'
      do ipj=0,2*m
         sum=0.
         if(ipj.eq.0)sum=1.
         do k=1,nr
            sum=sum+float(k)**ipj
         enddo
         do k=1,nl
            sum=sum+float(-k)**ipj
         enddo
         mm=min(ipj,2*m-ipj)
         do imj=-mm,mm,2
            a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
         enddo
      enddo
      call ludcmp(a,m+1,MMAX+1,indx,d)
      do j=1,m+1
         b(j)=0.
      enddo
      b(ld+1)=1.
      call lubksb(a,m+1,MMAX+1,indx,b)
      do kk=1,np
         c(kk)=0.
      enddo
      do k=-nl,nr
          sum=b(1)
          fac=1.
          do mm=1,m
             fac=fac*k
             sum=sum+b(mm+1)*fac
          enddo
          kk=mod(np-k,np)+1
          c(kk)=sum
      enddo
      return
      END



      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      INTEGER N,NP,INDX(N),NMAX
      REAL D,A(NP,NP),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER I,IMAX,J,K
      REAL AAMAX,DUM,SUM,VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix in ludcmp'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END


      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      INTEGER N,NP,INDX(N)
      REAL A(NP,NP),B(N)
      INTEGER I,II,J,LL
      REAL SUM
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END

       

