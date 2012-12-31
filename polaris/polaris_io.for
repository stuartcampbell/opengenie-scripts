*******************************************************************
* Stuart Campbell (ISIS/RAL) 28/7/98
*******************************************************************

	SUBROUTINE FILEGEN(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put

      INTEGER size
	PARAMETER(size=1000000)

	INTEGER np
	REAL x(size), y(size), e(size)

	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR) .EQ. 0) RETURN

      CALL MODULE_GET_INT(pars_get,'ntc',np)
	
      CALL PROSFILE(np,x,y,e,pars_get,pars_put)

	RETURN
	END


	SUBROUTINE PROSFILE(np,x,y,e,pars_get,pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put

	INTEGER i, unit, np, item, jtem, icount
	REAL x(np), y(np), e(np), tmin, tmax, xtmp
	CHARACTER filename*20, tmpstr*10

	CALL MODULE_GET_REAL_ARRAY(pars_get,'x',x,np)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'y',y,np)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'e',e,np)

	CALL MODULE_GET_REAL(pars_get,'tmin',tmin)		
	CALL MODULE_GET_REAL(pars_get,'tmax',tmax)		
	CALL MODULE_GET_STRING(pars_get,'filename',filename)		

	unit = 34
	item = 1
	jtem = 0
	icount = 0

	tmpstr = ' '
		
	OPEN(unit, FILE=filename, STATUS='NEW')
		
	DO i = 1, np
	
* MODIFIED BY JM 7-SEPT-1987 TO REMOVE ZEROS...AND MODIFIED AGAIN BY RIS
* 10-MARCH-1992 TO BEHAVE BETTER WITH ESD'S

          IF(y(i).LE.0.0) THEN
            y(i)=y(i-1)
            e(i)=e(i-1)
          ENDIF

          IF(y(i).LE.0.0) y(i)=0.001
          IF(e(i).LE.0.0) e(i)=0.001

		if (x(i).lt.tmin.or.x(i).gt.tmax) goto 999
		  icount=icount+1
		  xtmp= 0.5*(x(i)+x(i+1))

		  WRITE(unit,1200) xtmp,y(i),e(i),item,jtem
999	ENDDO

	CLOSE(unit)

	write(tmpstr,'(I5)') icount
*	write(tmpstr,'(A10)') icount

	CALL MODULE_INFORMATION('No. of data points  : '//tmpstr)

1000	FORMAT(A10)	
1200	FORMAT(X,F10.3,2(F10.4),2I5)

	RETURN
	END

*******************************************************************

	SUBROUTINE GMODULE_QUERY(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put

	INTEGER n
	PARAMETER (n=1)  
	CHARACTER*20 symbol(n),descrip(n),language(n),type(n),extent(n)
	
	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR) .EQ. 0) RETURN

	symbol(1)  = 'filegen'
	descrip(1) = 'Exports CCSL file'
	language(1) = 'FORTRAN'
	type(1) = 'TRANSFORM' 
	EXTENT(1) = 'DAT'

	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'symbol', symbol, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'descrip', descrip, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'language', language, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'type', type, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'extent', extent, n)
	RETURN	
	END

*******************************************************************
