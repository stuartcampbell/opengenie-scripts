	SUBROUTINE BIN1D_EXE(pars_get,pars_put)

	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
      INCLUDE 'canvas.inc'
	EXTERNAL pars_get, pars_put

	INTEGER num_p, num_e, size_x
      REAL selectx(ne,np), intensity(np,ne), error(np,ne)

	REAL int(np), err(np), count(np), x_values(np)

	WRITE(*,*) ' Entering BIN1D_EXE'

	CALL MODULE_GET_INT(pars_get,'num_e',num_e)             
	CALL MODULE_GET_INT(pars_get,'num_p',num_p)             
	CALL MODULE_GET_INT(pars_get,'size_x',size_x)             

	WRITE(*,*) 'SIZE_X = ',SIZE_X

	CALL BIN1D_SUB(num_e, num_p, size_x, int, err, count, 
     &		   x_values, selectx, intensity, error, 
     &		   pars_get, pars_put)
         
	   
	WRITE(*,*) ' Leaving BIN1D_EXE'
	               
	RETURN
	END


	SUBROUTINE BIN1D_SUB(num_e, num_p, size_x, int, err, count, 
     &		   x_values, selectx, intensity, error, 
     &		   pars_get, pars_put)
         
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put, bin1d

	INTEGER i, num_e, num_p, size_x, d, min_x, limit, flag
	REAL intensity(num_p,num_e), error(num_p,num_e), 
     &     selectx(num_e,num_p), new_p(3)
	REAL int(size_x), err(size_x), count(size_x), x_values(size_x)
      REAL tmpfunc, bin_size, bin1d      
	CHARACTER*80 info	


	WRITE(*,*) ' Entering BIN1D_SUB'

	DO i = 1, size_x
		int(i) = 0.0
		err(i) = 0.0
		count(i) = 0.0
		x_values(i) = 0.0
	ENDDO

	WRITE(*,*) 'SIZE_X = ',SIZE_X

	tmpfunc = 1.0
	d = 3
	limit = num_e*num_p

	CALL MODULE_GET_INT(pars_get,'min_x',min_x)

	CALL MODULE_GET_REAL_ARRAY(pars_get,'intensity',intensity,limit)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'error',error,limit)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'selectx',selectx,limit)

	CALL MODULE_GET_REAL_ARRAY(pars_get,'new_p',new_p,d)

	CALL MODULE_GET_REAL(pars_get,'del_x',bin_size)
	
	WRITE(*,*) ' Calling BIN1D'

	tmpfunc = BIN1D(intensity, error, num_e, num_p, min_x, size_x, 
     &		    new_p,bin_size,int,err,count,selectx,x_values,
     &		    info,flag)

	CALL MODULE_PUT_REAL_ARRAY(pars_put,'int',int,size_x)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'err',err,size_x)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'count',count,size_x)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'x',x_values,size_x)

	CALL MODULE_PUT_INT(pars_put,'flag',flag)
	CALL MODULE_PUT_STRING(pars_put,'info',info)


	WRITE(*,*) ' Leaving BIN1D_SUB'

	RETURN 
	END


