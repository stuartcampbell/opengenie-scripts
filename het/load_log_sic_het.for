 	SUBROUTINE load_log(pars_get,pars_put)

	include 'genie_modules.inc'

	integer			data_unit
	parameter			(mn=33000,data_unit=4)
	real				L1,L2
	integer			run_number,xcode,ycode
	character*200		ws_history
	character*80		long_title
	character*40		infile,xcaption,ycaption
	character*20		run_user,start_time
	character*8			inst_name
	common	/loader/ 	infile,lpt,x(mn),y(mn),e(mn),
     +				xcaption,ycaption,long_title,
     +				delta,inst_code,inst_name,
     +				xcode,ycode,L1,L2,tthet,
     +				no_spectrum,run_number,run_user,
     +				run_duration,start_time,fi,
     +				user_par(30),st_var(20),ws_history
	REAL		YY(32000)
      INTEGER 	status,TAR$SETDTIM,VERIFY,i
      INTEGER 	qdword(2)
	CHARACTER 	from_time*23, to_time*23,tstring*50
	CHARACTER*8 block_name
      CHARACTER*4 ans    
	CHARACTER*1 units
	real*8	myy(50000)

	DO i=1, mn
		x(i) = 0.0
		y(i) = 0.0
		e(i) = 1.0
	ENDDO

	CALL MODULE_GET_STRING(pars_get,'file',infile)
	CALL MODULE_GET_STRING(pars_get,'block_name',block_name)
	CALL MODULE_GET_STRING(pars_get,'units',units)
	CALL MODULE_GET_INT(pars_get,'logcol',icol)

	WRITE(block_name,'(A8)') block_name
	WRITE(units,'(A1)') units
	WRITE(icol,5000) icol


	write(*,*)'SE Block Name :',block_name
	write(*,*)'units  = ',units
	write(*,*)'column = ',icol

	OPEN(unit=data_unit,file=infile,status='old',readonly,shared)

	to_time = ' '
	from_time = ' '	        	

*	TYPE '(A,$)', ' Which log column do you want, 2 or 3 (Def 3)?'
*	READ(5,5000)icol

	if (icol.ne.2) icol=3

*       Get times from the log file in seconds from the start time

	Call str$upcase( from_time, from_time )
	Call str$upcase( to_time, to_time )
      Call str$upcase( block_name , block_name )

	if (icol.eq.2) then
		CALL C_get_log_values( x, y, yy, from_time, to_time,
     +                         block_name, lpt, data_unit, icol )
	else
		CALL C_get_log_values( x, yy, y, from_time, to_time,
     +                         block_name, lpt, data_unit, icol )
	endif

*	Decide whether to scale in days, hours, minutes or seconds
	IF ( x(lpt) .LE. 3 * 60.0 ) THEN
	    xcaption = 'Seconds from ' // from_time
	ELSEIF ( x(lpt) .LE. 3 * 3600.0 ) THEN
	    xcaption = 'Minutes from ' // from_time
	    DO I = 1, lpt
		x(i) = x(i) / 60.0
	    END DO
       	ELSEIF ( x(lpt) .LE. 3 * 3600.0 * 24.0 ) THEN
	    xcaption = 'Hours from ' // from_time
	    DO I = 1, lpt
		x(i) = x(i) / 3600.0
	    END DO
	ELSE
	    xcaption = 'Days from ' // from_time
	    DO I = 1, lpt
		x(i) = x(i) / ( 3600.0 * 24.0 )
	    END DO
	END IF
	    
	ycaption = 'Temperature (' // units // ')'
	Long_title = 'Plot of temperature from sample environment log file'

        lpt = lpt - 1

	  CALL MODULE_PUT_STRING(pars_put,'title',long_title)
	  CALL MODULE_PUT_STRING(pars_put,'ylabel',ycaption)
	  CALL MODULE_PUT_STRING(pars_put,'xlabel',xcaption)
	  CALL MODULE_PUT_STRING(pars_put,'spec_no',block_name)
	  CALL MODULE_PUT_INT(pars_put,'ntc',lpt)
	  CALL MODULE_PUT_REAL_ARRAY(pars_put,'x',x,lpt+1)
	  CALL MODULE_PUT_REAL_ARRAY(pars_put,'y',y,lpt)

	CLOSE(data_unit)

 5000	format(i1)

	END



        INTEGER FUNCTION VERIFY
        INTEGER x
        INTEGER y
        CHARACTER ans*4

        PRINT *, ' Is this what you want y/n ? '
        ACCEPT '(a)', ans
        CALL str$upcase(ans,ans)
        x = INDEX(ans,'Y')
        y = INDEX(ans,'N')
        DO WHILE ((x .EQ. 0).AND. (y .EQ. 0))
             x = INDEX(ans,'Y')
             y = INDEX(ans,'N')
             TYPE *, ' what (y/n) ? '
             ACCEPT *,ans
             CALL str$upcase(ans,ans)
        END DO
           
        IF ( x .EQ. 0) THEN
           VERIFY = 2
         ELSE 
           VERIFY = 1
        END IF
        RETURN
        END

c
c ------------------------------------------------------------------------------
	SUBROUTINE C_GET_LOG_VALUES(TIME,VALUE,VALUE1,START,FIN,KEY,TOTAL,IUNIT)
C
C This subroutine will extract numeric values associated with a specific key
C form the instrument log file.
C 
C If the start time string is blank then the subroutine uses the first entry
C as the intial value, and similarly if the finish string is blank then all
C values before the end of file are returned
C
C  modified Dec 85 to use new log file format that stores month as a string
C  instead of a numeric value
C bug fix Feb88 KJK: TOTAL=TOTAL+1 moved to prevent funny values when program
C		     detects an error in the input
C
	CHARACTER*(*) START,FIN,KEY
	CHARACTER IP*40
	CHARACTER ADATE*20,IDENT*5,SEBLOCK*8
	REAL*4 VALUE(*),TIME(*),VALUE1(*)
	INTEGER*4 STIME(2),FTIME(2),ETIME(2),YTIME(2)
	INTEGER*4 TOTAL,TIMCMP,NSTRT,NFIN,IDATE,IUNIT,ITIME,IVALUE,ST,FT,
     &	    ET,REM,temp
	INTEGER*4 ICOUNT,IHOURS,IMINS,ISECS,ID,STATUS,start_check,
     &	    end_check,IVALU1
	LOGICAL*4 DEFAULT_START,DEFAULT_FIN
c
	ID=10000000
	TOTAL=0
	start_check = 0 
	if (start .eq. '                  ') then
		start_check = 1
	end if
c
	if (fin .eq. '                  ') then
		 fin = '31-DEC-1999 23:59'
	end if
c
	CALL SYS$BINTIM(START,STIME)
	CALL SYS$BINTIM(FIN,FTIME)
c
9	STATUS = LIB$SUBX(FTIME,STIME,YTIME)
        STATUS = LIB$EDIV(ID, YTIME, FT, REM)
	ST=0
c
10	READ(IUNIT,3010,END=200)ADATE,IDENT,SEBLOCK,IP
3010	FORMAT(4X,A20,1x,a5,1x,A8,A40)
c
	IF (start_check .eq.  1) then
		start = adate 
		CALL SYS$BINTIM(START,STIME)
		start_check = 0
c 25-2-97: Following line removed so as first point in the log file gets
c into the plot...
c		goto 9
	end if	
c
	IF(IDENT.NE.'CAMAC') GOTO 10
	IF(SEBLOCK.NE.KEY) GOTO 10
	CALL SYS$BINTIM(ADATE,ETIME)
	status = LIB$SUBX(ETIME,STIME,YTIME)
	IF(BITEST(YTIME(2),31)) GOTO 10
        STATUS = LIB$EDIV(ID, YTIME, ET, REM)
c
c 25-2-97: The following line removed as it seemed only to cause problems.
c	IF(ET.GT.FT) GOTO 200
c
	READ(IP,3015,ERR=900) IVALUE,IVALUE1
3015	FORMAT(14X,I8,3X,I8)
	TOTAL=TOTAL+1
	VALUE(TOTAL)=FLOAT(IVALUE)/100.0	! convert out of engineering units
	VALUE1(TOTAL)=FLOAT(IVALUE1)/100.0	! convert out of engineering units
	TIME(TOTAL)=FLOAT(ET)
	GOTO 10
c
 900	CONTINUE
	TYPE *, 'Invalid format'
	GOTO 10
c
c 200	CLOSE(IUNIT)   OLD
 200	RETURN
	END
                                                                                                     