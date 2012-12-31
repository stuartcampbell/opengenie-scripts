	SUBROUTINE GETPAR_3D(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules_f90.inc'
	EXTERNAL pars_get, pars_put

	INCLUDE 'canvas_f90.inc'

	INTEGER 		:: i, j, unit, stat, num_p, number
	REAL, ALLOCATABLE	:: dist(:),phi(:),theta(:),width(:),height(:)
	CHARACTER(LEN=80)	:: input_file

! Check for GENIE version...
	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR).EQ.0) THEN
	  CALL MODULE_ERROR('GSASGEN','ERROR - Open GENIE version mismatch'& 
     & ,'Please re-compile this module')
!		RETURN
	ENDIF

	input_file=' '

	CALL MODULE_GET_STRING(pars_get, 'par_file', input_file)
	CALL MODULE_GET_INT(pars_get, 'num_p', num_p)
	
	unit = 34
	stat = 0
	
	OPEN(unit,FILE=input_file,STATUS='OLD',readonly,shared,IOSTAT=stat)
	IF(stat.ne.0) THEN
		CALL module_error('GETPAR_3D','Unable to open ' &
     &           //input_file,'Check file exists and is readable')
		CALL MODULE_PUT_INT(pars_put,'stat',stat)
		RETURN
	ENDIF
	
	READ(unit,*) number
	
      IF(number.ne.num_p) THEN
	CALL module_error('GETPAR','PAR file does not match data', &
     &	'Check file name is correct')
		RETURN
	ENDIF
	
	CALL MODULE_INFORMATION('Opening parameter file... '//input_file)

 	ALLOCATE(dist(num_p),phi(num_p),theta(num_p),width(num_p),height(num_p))

	CALL READ_PAR(num_p,dist,phi,theta,width,height,unit,pars_put)
	
	close(unit,iostat=stat)

	CALL MODULE_PUT_INT(pars_put,'stat',stat)

	RETURN
	END
	
                            
	SUBROUTINE READ_PAR(n,dist,phi,theta,width,height,unit,pars_put)
	IMPLICIT NONE
	INCLUDE 'genie_modules_f90.inc'
	EXTERNAL pars_put

	INTEGER n, i, unit, iaddr(2)
!	INTEGER*8 ioff
	REAL dist(n), phi(n), theta(n), width(n), height(n)

	DO i = 1, n
	 READ(unit,*) dist(i),phi(i),theta(i),width(i),height(i)
	ENDDO

	CALL MODULE_PUT_REAL_ARRAY(pars_put,'dist',dist,n)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'phi',phi,n)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'theta',theta,n)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'width',width,n)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'height',height,n)

	RETURN
	END



	SUBROUTINE GMODULE_QUERY(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules_f90.inc'
	EXTERNAL pars_get, pars_put

	INTEGER n
	PARAMETER (n=1)  
	CHARACTER*20 symbol(n),descrip(n),language(n),type(n),extent(n)
	
	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR) .EQ. 0) RETURN

	symbol(1)  = 'getpar_3d'
	descrip(1) = 'Gets .PAR files'
	language(1) = 'FORTRAN'
	type(1) = 'INPUT' 
	EXTENT(1) = 'PAR'

	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'symbol', symbol, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'descrip', descrip, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'language', language, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'type', type, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'extent', extent, n)
	RETURN	
	END



