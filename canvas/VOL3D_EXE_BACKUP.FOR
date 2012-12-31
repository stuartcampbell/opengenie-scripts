	SUBROUTINE VOL3D_EXE(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put

      INCLUDE 'canvas.inc'

	INTEGER num_p, num_e

	REAL phi(np), theta(np), dist(np), energy(ne+1),
     &     vx(ne,np), vy(ne,np), vz(ne,np), 
     &     vxt(np,ne), vyt(np,ne), vzt(np,ne) 

	CALL MODULE_GET_INT(pars_get,'num_e',num_e)
	CALL MODULE_GET_INT(pars_get,'num_p',num_p)

	CALL VOL3D_SUB(num_e,num_p,phi,theta,dist,energy,vx,vy,vz,
     &		   vxt,vyt,vzt,pars_get,pars_put)
                                                
	RETURN
	END
                 


	SUBROUTINE VOL3D_SUB(num_e,num_p,phi,theta,dist, energy, 
     &			   vx, vy, vz, vxt, vyt, vzt, pars_get, 
     &			   pars_put)

	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'

	INTEGER i, j, num_e, num_p, flag, axes(3), dims_array(2)
	REAL ei, psi, tmp, vol3d
	REAL phi(num_p), theta(num_p), dist(num_p), energy(num_e+1),
     &     vx(num_e,num_p), vy(num_e,num_p), vz(num_e,num_p),
     &     vxt(num_p,num_e), vyt(num_p,num_e), vzt(num_p,num_e)
	CHARACTER*80 info
 	EXTERNAL pars_get, pars_put, vol3d

	CALL MODULE_GET_STRING(pars_get,'info',info)
	CALL MODULE_GET_INT(pars_get,'flag',flag)
	CALL MODULE_GET_INT_ARRAY(pars_get,'axes',axes,3)

	CALL MODULE_GET_REAL(pars_get,'ei',ei)
	CALL MODULE_GET_REAL(pars_get,'psi',psi)

	CALL MODULE_GET_REAL_ARRAY(pars_get,'phi',phi,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'theta',theta,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'dist',dist,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'energy',energy,num_e+1)

	tmp = VOL3D(phi, theta, dist, energy, ei, psi, num_e,
     &		num_p, vx, vy, vz, axes, info, flag)
	
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'phi',phi,num_p)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'theta',theta,num_p)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'dist',dist,num_p)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'energy',energy,num_e+1)

	CALL MODULE_PUT_REAL(pars_put,'ei',ei)
	CALL MODULE_PUT_REAL(pars_put,'psi',psi)
	CALL MODULE_PUT_INT(pars_put,'flag',flag)
	CALL MODULE_PUT_INT_ARRAY(pars_put,'axes',axes,3)

	CALL MODULE_PUT_INT(pars_put,'num_e',num_e)
	CALL MODULE_PUT_INT(pars_put,'num_p',num_p)

	
	DO i = 1, num_p
	  DO j = i, num_e
		vxt(j,i) = vx(i,j)
		vyt(j,i) = vy(i,j)
		vzt(j,i) = vz(i,j)
	  ENDDO
	ENDDO

	dims_array(1) = num_p
	dims_array(2) = num_e
                  
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'vx',vxt,dims_array,2)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'vy',vyt,dims_array,2)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'vz',vzt,dims_array,2)

	RETURN
	END

