
	SUBROUTINE GETRUN_3D(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put

	INCLUDE 'canvas.inc'

	INTEGER unit, stat, num_p, num_e

	REAL e_points(ne+1),p_points(np+1),data(np,ne),error(np,ne)
	CHARACTER*80 input_file

* Check for GENIE version...
	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR).EQ.0) THEN
	  CALL MODULE_ERROR('GSASGEN','ERROR - Open GENIE version mismatch' 
     &  ,'Please re-compile this module')

	CALL MODULE_PUT_INT(pars_put,'stat',stat)
C	  RETURN
	ENDIF

	input_file=' '

	CALL MODULE_GET_STRING(pars_get, 'file', input_file)

	unit = 34
	stat = 0
		
	OPEN(unit,FILE=input_file,STATUS='OLD',readonly,shared,IOSTAT=stat)
	IF(stat.ne.0) THEN
		CALL module_error('GETRUN','Unable to open ' 
     &           //input_file,'Check file exists and is readable')
		RETURN
	ENDIF

* Read data from file

	READ(unit,*) num_p, num_e         
		
	CALL READ_RUN(num_e, num_p, e_points, p_points,
     &		   data, error, unit, stat, pars_put)	


	CALL MODULE_PUT_INT(pars_put,'stat',stat)

	RETURN
	END



	SUBROUTINE READ_RUN(num_e, num_p, e_points, p_points,
     &		         data, error, unit, stat, pars_put)

	IMPLICIT NONE
	INCLUDE 'genie_modules.inc'
	EXTERNAL pars_put
	INTEGER i, j, dims_array(2), unit, stat
	INTEGER num_e, num_p
	REAL e_points(num_e+1), p_points(num_p+1)
	REAL data(num_p,num_e), error(num_p,num_e)
	
* Read in the data set grids

	READ(unit,*)
	READ(unit,100) (p_points(i),i=1,num_p+1)	
	
	READ(unit,*) 
	READ(unit,100) (e_points(i),i=1,num_e+1)
	
	DO i = 1, num_p

* Read in as (num_p,num_e), which is the 'wrong' way
* round, but it will be reversed when passed back into GCL.

		READ(unit,*)
		READ(unit,100) (data(i,j),j=1,num_e)
		READ(unit,*)
		READ(unit,100) (error(i,j),j=1,num_e)		
	ENDDO	
	
* Return data to workspace

	dims_array(1) = num_e
	dims_array(2) = num_p
	                            	
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'error',error,
     &	dims_array,2)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'intensity',data,
     &	dims_array,2)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'energy',e_points,num_e+1)
*	CALL MODULE_PUT_REAL_ARRAY(pars_put,'phi',p_points,num_p+1)
	CALL MODULE_PUT_INT(pars_put,'num_e',num_e)
	CALL MODULE_PUT_INT(pars_put,'num_p',num_p)

	close(unit,iostat=stat)
	
100	FORMAT(8G10.4)

	RETURN
	END


	
	SUBROUTINE GMODULE_QUERY(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put

	INTEGER n
	PARAMETER (n=1)  
	CHARACTER*20 symbol(n),descrip(n),language(n),type(n),extent(n)
	
	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR) .EQ. 0) RETURN

	symbol(1)  = 'getrun_3d'
	descrip(1) = 'Gets .SPE files'
	language(1) = 'FORTRAN'
	type(1) = 'INPUT' 
	EXTENT(1) = 'SPE'

	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'symbol', symbol, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'descrip', descrip, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'language', language, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'type', type, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'extent', extent, n)
	RETURN	
	END


