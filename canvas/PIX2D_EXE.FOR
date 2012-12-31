	SUBROUTINE PIX2D_EXE(pars_get,pars_put)

	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
      INCLUDE 'canvas.inc'
	EXTERNAL pars_get, pars_put

	INTEGER num_p, num_e
	REAL phi(np), width(np), theta(np), dist(np), energy(ne+1),
     &     x(ne,np,5), y(ne,np,5)             
	REAL c_x(np*ne,4)
	REAL c_y(np*ne,4)
                  


	CALL MODULE_GET_INT(pars_get,'num_p',num_p)
	CALL MODULE_GET_INT(pars_get,'num_e',num_e)
      
	CALL PIX2D_SUB(num_e, num_p, phi, width, theta,  
     &		   dist, energy, x, y, pars_get, pars_put,
     &		   c_x, c_y)
              

	RETURN
	END


	SUBROUTINE PIX2D_SUB(num_e, num_p, phi, width, theta, 
     &			   dist,energy,x,y,pars_get, 
     &			   pars_put,c_x,c_y)

	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put, pix2d

	INTEGER num_p, num_e, axes(3), d, flag, dims_array(2)
	INTEGER dims_array2(2), index, i,j,k
	REAL ei, psi, tmp, aspect, scale(3), pix2d
	REAL phi(num_p), width(num_p), theta(num_p), dist(num_p), 
     &     energy(num_e+1), x(num_e,num_p,5), y(num_e,num_p,5)
	REAL c_x(num_p*num_e,4), c_y(num_p*num_e,4)
	REAL centre_x(num_p*num_e), centre_y(num_p*num_e)

	CHARACTER*80 info                  

	d = 3
	flag = 0
	tmp=1.0

	CALL MODULE_GET_REAL_ARRAY(pars_get,'scale',scale,d)
      CALL MODULE_GET_STRING(pars_get,'info',info)
	CALL MODULE_GET_INT(pars_get,'flag',flag)
	CALL MODULE_GET_INT_ARRAY(pars_get,'axes',axes,d)

	CALL MODULE_GET_REAL(pars_get,'ei',ei)
	CALL MODULE_GET_REAL(pars_get,'psi',psi)
	CALL MODULE_GET_REAL(pars_get,'aspect',aspect)

	CALL MODULE_GET_REAL_ARRAY(pars_get,'phi',phi,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'width',width,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'theta',theta,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'dist',dist,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'energy',energy,num_e+1)

	print*,' aspect = ',aspect
	print*,' scale  = ',scale 

	tmp=PIX2D(phi,width,theta,dist,energy,ei,psi,scale,aspect,
     &          num_e,num_p,x,y,axes,info,flag)

	print*,' aspect = ',aspect
	print*,' scale  = ',scale 

	CALL MODULE_PUT_REAL(pars_put,'aspect',aspect)
	CALL MODULE_PUT_INT(pars_put,'flag',flag)
	CALL MODULE_PUT_INT(pars_put,'num_e',num_e)
	CALL MODULE_PUT_INT(pars_put,'num_p',num_p)
	CALL MODULE_PUT_INT_ARRAY(pars_put,'axes',axes,d)

	dims_array(1) = num_e
	dims_array(2) = num_p
	dims_array(3) = 5
                     	  
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'x',x,dims_array,3)	
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'y',y,dims_array,3)	
	
	do k = 1, 5
	 do i = 1, num_e
	  do j = 1, num_p
	   index = ((i-1)*num_p)+j
         IF (k.eq.1) THEN
	    centre_x(index) = x(i,j,k)
	    centre_y(index) = y(i,j,k)
	   ELSE
	    c_x(index,k-1) = x(i,j,k)
	    c_y(index,k-1) = y(i,j,k)
         ENDIF
	  enddo
	 enddo
	enddo

	dims_array2(1) = num_e*num_p
	dims_array2(2) = 4

	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'c_x',c_x,dims_array2,2)	
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'c_y',c_y,dims_array2,2)	

	CALL MODULE_PUT_REAL_ARRAY(pars_put,'centre_x',centre_x,num_e*num_p)	
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'centre_y',centre_y,num_e*num_p)	
	
	RETURN
	END
