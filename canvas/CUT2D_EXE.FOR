	SUBROUTINE CUT2D_EXE(pars_get,pars_put)

	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
      INCLUDE 'canvas.inc'
	EXTERNAL pars_get, pars_put

	INTEGER num_p, num_e
      REAL vx(ne,np,5), vy(ne,np,5), selectx(ne,np)             
                  
	CALL MODULE_GET_INT(pars_get,'num_p',num_p)
	CALL MODULE_GET_INT(pars_get,'num_e',num_e)
      
	CALL CUT2D_SUB(num_e, num_p, selectx, vx, vy, pars_get, pars_put) 
                          
	RETURN
	END



	SUBROUTINE CUT2D_SUB(num_e, num_p, selectx, vx, vy, 
     &			   pars_get, pars_put)

	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put, cut2d

	
        INTEGER num_e, num_p, axes(3), d, min_x, size_x, 
     &	  flag, length, dims_array(2) 
	REAL tmpfunc, width, param, cut_bin, cut2d, q(2), r(2), new_p(3)
	REAL vx(num_e,num_p,5), vy(num_e,num_p,5), selectx(num_e,num_p)
 	CHARACTER*80 info                 

	d = 3
      tmpfunc = 0.0
	flag = 1

      q(1) = 0.0
	q(2) = 0.0
      r(1) = 0.0
	r(2) = 0.0

	length = 5*num_e*num_p

	CALL MODULE_GET_REAL_ARRAY(pars_get,'x',vx,length)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'y',vy,length)

      CALL MODULE_GET_INT_ARRAY(pars_get,'axes',axes,d)

	d = 2
	
      CALL MODULE_GET_REAL_ARRAY(pars_get,'x_value',q,d)
      CALL MODULE_GET_REAL_ARRAY(pars_get,'y_value',r,d)

	write(*,*) ' Origin of cut = ',q
	write(*,*) ' towards ...     ',r

      CALL MODULE_GET_INT(pars_get,'mode',param)

      CALL MODULE_GET_REAL(pars_get,'del_y',width)
      CALL MODULE_GET_REAL(pars_get,'del_x',cut_bin)
	 
      tmpfunc = CUT2D(axes,vx,vy,selectx,num_e,num_p,min_x, 
     &		    size_x,q,r,cut_bin,width,param,new_p,info,flag)

	WRITE(*,*) ' SIZE_X = ',SIZE_X

	CALL MODULE_PUT_STRING(pars_put,'info',info)
	CALL MODULE_PUT_INT(pars_put,'num_e',num_e)
	CALL MODULE_PUT_INT(pars_put,'num_p',num_p)
	CALL MODULE_PUT_INT(pars_put,'flag',flag)
	CALL MODULE_PUT_INT(pars_put,'min_x',min_x)
	CALL MODULE_PUT_INT(pars_put,'size_x',size_x)
      CALL MODULE_PUT_REAL(pars_put,'del_x',cut_bin)

	d = 3
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'new_p',new_p,d)

	write(*,*) 'new_p (cut2d_exe)= ',new_p

	dims_array(1) = num_e
	dims_array(2) = num_p

	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'selectx',selectx,
     &	  				dims_array,2)

	RETURN
	END
