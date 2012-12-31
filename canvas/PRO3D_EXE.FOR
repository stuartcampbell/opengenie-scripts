
	SUBROUTINE PRO3D_EXE(pars_get, pars_put)      
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
      INCLUDE 'canvas.inc'
	EXTERNAL pars_get, pars_put, pro3d_sub

	INTEGER num_p, num_e, total

      REAL vx(ne*np), vy(ne*np), vz(ne*np), temp(ne*np)
      REAL selectx(ne,np), selecty(ne,np)

	CALL MODULE_GET_INT(pars_get,'num_e',num_e)
	CALL MODULE_GET_INT(pars_get,'num_p',num_p)

	total = num_e*num_p

	CALL MODULE_GET_REAL_ARRAY(pars_get,'vx',vx,total)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'vy',vy,total)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'vy',vy,total)

	CALL PRO3D_SUB(num_e,num_p,selectx,selecty,
     &		   vx, vy, vz, pars_get, pars_put,temp)
                                                
	RETURN
	END
                 


	SUBROUTINE PRO3D_SUB(num_e, num_p, selectx, selecty, 
     &			   vx, vy, vz, pars_get, pars_put, temp)

	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
 	EXTERNAL pars_get, pars_put, pro3d

	INTEGER i, j, num_e, num_p, flag, axes(3), dims_array(2), d
	INTEGER size_x, size_y, min_x, min_y, mode
	INTEGER p1(3), p2(3), p3(3), p4(3), pnvector(3)

	REAL tmp, pro3d, temp(num_e*num_p)
      REAL selectx(num_e,num_p), selecty(num_e,num_p)
      REAL vx(num_e,num_p), vy(num_e,num_p), vz(num_e,num_p)

	CHARACTER*80 info

	tmp = 1.0
	flag = 0
	d = 3
	min_x = 0
	min_y = 0
	size_x = 0
	size_y = 0
            
	CALL MODULE_GET_STRING(pars_get,'info',info)
	CALL MODULE_GET_INT(pars_get,'flag',flag)
	CALL MODULE_GET_INT(pars_get,'mode',mode)
	CALL MODULE_GET_INT_ARRAY(pars_get,'axes',axes,d)
	CALL MODULE_GET_INT_ARRAY(pars_get,'p1',p1,d)
	CALL MODULE_GET_INT_ARRAY(pars_get,'p2',p2,d)
	CALL MODULE_GET_INT_ARRAY(pars_get,'p3',p3,d)
	CALL MODULE_GET_INT_ARRAY(pars_get,'p4',p4,d)
      
	tmp=PRO3D(p1, p2, p3, pnvector, p4, mode, axes,
     &	    vx, vy, vz, selectx, selecty, num_e, num_p, 
     &          min_x, min_y, size_x, size_y, info, flag)
	
	CALL MODULE_PUT_INT(pars_put,'flag',flag)
	CALL MODULE_PUT_INT_ARRAY(pars_put,'axes',axes,d)
	CALL MODULE_PUT_INT(pars_put,'num_e',num_e)
	CALL MODULE_PUT_INT(pars_put,'num_p',num_p)
	CALL MODULE_PUT_INT(pars_put,'min_x',min_x)
	CALL MODULE_PUT_INT(pars_put,'min_y',min_y)
	CALL MODULE_PUT_INT(pars_put,'size_x',size_x)
	CALL MODULE_PUT_INT(pars_put,'size_y',size_y)
	CALL MODULE_PUT_INT_ARRAY(pars_put,'pnvector',pnvector,d)
         
	dims_array(1) = num_e
	dims_array(2) = num_p
                  
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'selectx',selectx,dims_array,2)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'selecty',selecty,dims_array,2)

	CALL SWAP(selectx,temp,num_e,num_p)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'selectx',temp,num_e*num_p)
	CALL SWAP(selecty,temp,num_e,num_p)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'selecty',temp,num_e*num_p)
	

	RETURN
	END

	SUBROUTINE swap(in,temp,idim,jdim)

	REAL temp(jdim*idim), in(idim,jdim)
	INTEGER i, j, k, idim, jdim

	DO i = 1, idim
 	 DO j = 1, jdim
	   k = ((i-1)*idim) + j
	   temp(k) = in(i,j)
  	 ENDDO
	ENDDO

	RETURN
	END
