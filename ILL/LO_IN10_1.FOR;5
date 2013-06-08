	PROGRAM LO_IN10				!loads IN10 data into w-space
C
	include	'genie_sources:LOADCOM.CMN'

	character line*80, lin60*60,lin20*20
	integer nchan
	parameter (nchan=256)
      INTEGER D(nchan),M1(nchan)
      REAL MAX,MAXH,NUE
      DIMENSION DY(nchan),Y1(nchan),VE(nchan),YI(nchan),
     1 DF(10),KT(8),ang(8)
C ************
      DATA MDIV,CONV,YLEVEL/50,0.06595,0.3/

	call LOAD_IN
	num_sp=8
	user_par(1)=FLOAT(num_sp)			!number of spectra
	midch=nchan/2 +1
	ktr=1
	inst_name(1:3)='INX'
	run_user='IN10'
	xcode=5					!code for energy
	xcaption=' Energy micro-eV '
	ycaption='         Intensity          '

	INRM=1					!to normalise to mon
	read(data_unit,1102)line		!line of RRR
	read(data_unit,1103)inst_name(4:8)	!run number
	read(data_unit,1102)line		!line of AAA
	read(data_unit,1102)line		!line starting 80
	read(data_unit,1102)line		!line starting IN10
	read(data_unit,1102)line		!line of III
	read(data_unit,1104)jl			!#lines if I
	jline=jl/10
	do j=i,jline
	 read(data_unit,1102)line		!lines of 10 integers
	end do
	read(data_unit,1102)line		!line of AAA
	read(data_unit,1102)lin20		!line starting 512
	read(data_unit,1102)lin60
	long_title(1:60)=lin60(1:60)		!title line
	WRITE(6,1007)long_title(1:40)
c	do j=i,6
c	 read(data_unit,1102)line		!6 lines of text
c	end do
c	start_time(1:9)=line(21:29)
c	start_time(10:10)=' '
c	start_time(11:18)=line(11:18)
	read(data_unit,1102)line		!line of FFF
	read(data_unit,1104)jl			!#points
	jline=jl/5
	read(data_unit,1105)a,NUE
c	read(data_unit,1070)KTR,(KT(K),K=1,KTR)	!detector numbers
c	read(data_unit,1071)NUE,KON,IDR		!freq,xtal,drive
	write(6,1008)NUE
	do j=1,jline
	 read(data_unit,1105)a		!lines of 5 reals
	end do
	read(data_unit,1102)line		!line of FFF
	read(data_unit,1104)jl			!#points
	read(data_unit,1105)(ang(k),k=1,8)
	jline=jl/5 -1
	do j=1,jline
	 read(data_unit,1105)a		!lines of 5 reals
	end do
	do K=1,KTR				!loop over detectors
	 read(data_unit,1102)line		!line of SSS
	 read(data_unit,1104)j1,j2,jt		!#points
	 write(6,1009)j1,jt,ang(j1)
	 read(data_unit,1102)line		!line of III
	 read(data_unit,1104)lpt		!#points
	 read(data_unit,1104)(D(l),l=1,lpt)	!detector data
	end do
	do K=1,j2				!loop over detectors
	 read(data_unit,1102)line		!line of SSS
	 read(data_unit,1104)j1,j2,jt		!#points
	 write(6,1009)j1,jt
	 read(data_unit,1102)line		!line of III
	 read(data_unit,1104)nch		!#points
	 read(data_unit,1104)(M1(l),l=1,lpt)	!detector data
	end do
	lpt=nch
	idr=1
	kon=0
	if(IDR.EQ.0)then			!slow drive
	 DRV=19.656769
	else
	 if(IDR.EQ.1)then			!fast drive
	  DRV=18.706
	 else					!bad number
	  stop ' lo_inx> drive not defined'
	 endif
	endif
	VMAX=NUE*DRV				!max speed

	if(KON.EQ.0)then			!Si111 xtal
	 CONV=0.06595
c	 MDIV=50
	 MDIV=1
	  wav=2.*3.1354				!2 times d-spacing
	else
	 if(KON.EQ.1)then			!Si311 xtal
	  CONV=0.126284
	  MDIV=2
	 else					!xtal not defined
	  write(6,1006)				!ask for input
	  accept 1032,CONV
	 endif
	endif
	E0=81.787/(wav*wav)
	user_par(3)=E0
	user_par(4)=ang(ktr)			!2theta in deg
	tthet=ang(ktr)*3.14159265/180.
	user_par(5)=4.*3.14159265*SIN(0.5*tthet)/wav	!Q
	VMAX=VMAX*CONV				!max energy
	SUMI=0.
	DO I=1,nch				!calc total count
	 YI(I)=M1(I)
	 SUMI=SUMI+YI(I)
	end do
	SUMM=0.					!calc sum over flat part
	DO I=midch-25,midch+25
	 SUMM=SUMM+YI(I)
	end do
	FAKT1=50./SUMM
	FAKTM=1.0E-3
	if(INRM.EQ.0)then			!normalise to mon sum flat
	 DO I=1,nch
	  YI(I)=YI(I)*FAKT1
	 end do
	 FAK=FAKT1
	else					!norm to 1e-3
	 DO I=1,nch
	  YI(I)=YI(I)*FAKTM
	 end do
	 FAK=FAKTM
	end if
	WRITE(6,1016)SUMI			!print total count

	NH=lpt/2				!mid-point
	NHH=1					!search for max monitor -ve
	max=0.
	do n=1,NH
	 if(YI(n).gt.max) then
	  max=YI(n)
	  NU=n
	 endif
	end do
	NHH=NU
	MAXH=max/2.
  330 NHH=NHH-1
      IF( NHH .LT. 2 ) GO TO 320
      IF( YI(NHH) .GT. MAXH ) GO TO 330
      IF( YI(NHH) .LT. YLEVEL*MAX ) NHH=NHH+1
  320 NU=NHH

      NHH=lpt					!search for max mon +ve
	max=0.
	do n=NH,NHH
	 if(YI(n).gt.max) then
	  max=YI(n)
	  NO=n
	 endif
	end do
	NHH=NO
	MAXH=max/2.
  430 NHH=NHH+1
      IF( NHH .GE. lpt-1) GO TO 420
      IF( YI(NHH) .GT. MAXH ) GO TO 430
      IF( YI(NHH) .LT. YLEVEL*MAX ) NHH=NHH-1
  420 NO=NHH
      WRITE(6,1018)NU,NO
	NU=NU-1
	lpt=NO-NU

	DO I=1,lpt
	 YI(I)=YI(I+NU)
	end do
	TIFN=FLOAT(lpt)
	TIFN1=(TIFN-1.)*0.5
	DELTA=2.*VMAX/(TIFN-1.)
	DO I=1,lpt
	 x(I)=(I-midch+NU)*DELTA
	 VE(I)=x(I)
	end do

	SUME=0.					!calc total counts
	DO I=1,nch
	 SUME=SUME+D(I)
	end do
	FL=0.0
	DO I=1,lpt
	 D(I)=D(I+NU)				!move data down
	 x(I)=VE(I)*1e-3			!x=energy in meV
	 y(I)=D(I)/YI(I)			!y=d/mon
	 FL=FL+y(I)				!normalised sum
	end do
	FL=FL*FAK
	FL1=FL/FAKT1
	FL2=FL/FAKTM
	WRITE(6,1017)no_spectrum,SUME

	DO I=1,lpt
	 DYE=ABS(FLOAT(D(I)))
	 if(DYE.LT.0.00001)then
	  DYE=1.
	 endif
	 DYI=ABS(YI(I))*FAK/MDIV
         DY(I)=(1./YI(I)**2)*DYE + ((D(I)/YI(I)**2)**2)*DYI
	 e(I)=SQRT(DY(I))
	end do
C
	call LOAD_OUT
c
 	stop ' '
1005  FORMAT(10I7)
1006	FORMAT(' lo_inx> xtal conversion factor ?',$)
1007  FORMAT(' lo_inx> run title : ',a40)
1008  FORMAT(' lo_inx> Doppler frequency : ',f7.2)
1009  FORMAT(' lo_inx> reading spectrum ',i3,' of ',i3,
	1 ' angle : ',f7.2)
1016  FORMAT(' lo_inx> Total Monitor Count :',E14.6)
1017  FORMAT(' lo_inx> Detector',I2,' ; Total Sample Count :',E14.6)
1018  FORMAT(' lo_inx> Maximum Monitor Intensity in channels ',i3,
	1 ' & ',i3)
1032	FORMAT(8f)
 1070 FORMAT(I5,24X,8I3)
 1071 FORMAT(8X,F6.3,14X,I1,10X,I1)
 1073 FORMAT(10X,I2)
1100	FORMAT(' which spectrum ?',$)
1101	FORMAT(q,i)
1102	FORMAT(a)
1103	FORMAT(3x,a5)
1104	FORMAT(10i8)
1105	FORMAT(5e16.8)
	END
