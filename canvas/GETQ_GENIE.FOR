	SUBROUTINE GETQ_GENIE(PARS_GET, PARS_PUT)

	IMPLICIT NONE
	INCLUDE 'genie_modules.inc'
	INCLUDE 'canvas_genie.inc'
	EXTERNAL PARS_GET, PARS_PUT

	INTEGER NUM_E, NUM_P
	REAL PHI(P_DIMEN), THETA(P_DIMEN), ENERGY(E_DIMEN+1), DIST(P_DIMEN)
	REAL QX(E_DIMEN, P_DIMEN), QY(E_DIMEN, P_DIMEN), QZ(E_DIMEN, P_DIMEN),
     &     ENG(E_DIMEN, P_DIMEN)  
	REAL EI, PSI, APHI, ATHETA, APSI

	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR) .EQ. 0) RETURN

	CALL MODULE_GET_INT(PARS_GET, 'NUM_E', NUM_E)		
	CALL MODULE_GET_INT(PARS_GET, 'NUM_P', NUM_P)		
	CALL MODULE_GET_REAL(PARS_GET, 'EI', EI)		
	CALL MODULE_GET_REAL(PARS_GET, 'PSI', PSI)		
	CALL MODULE_GET_REAL(PARS_GET, 'APHI', APHI)		
	CALL MODULE_GET_REAL(PARS_GET, 'ATHETA', ATHETA)		
	CALL MODULE_GET_REAL(PARS_GET, 'APSI', APSI)
	CALL MODULE_GET_REAL_ARRAY(PARS_GET, 'PHI', PHI, NUM_P)		
	CALL MODULE_GET_REAL_ARRAY(PARS_GET, 'DIST', DIST, NUM_P)		
	CALL MODULE_GET_REAL_ARRAY(PARS_GET, 'THETA', THETA, NUM_P)		
	CALL MODULE_GET_REAL_ARRAY(PARS_GET, 'ENERGY', ENERGY, NUM_E+1)		


	CALL CALCQ(PHI,THETA,DIST,ENERGY,EI,PSI,NUM_E,NUM_P,QX,QY,QZ,ENG,
     &	     APHI,ATHETA,APSI,PARS_PUT,PARS_GET)

                               
	RETURN
	END

C
C  CALCQ is a FORTRAN routine for LEONARDO
C  It calclates the Q vector for each detector and energy transfer
C
C

      SUBROUTINE CALCQ(phi, theta, dist, energy, ei, psi, 
     &                      num_e, num_p, qx, qy, qz, eng, 
     &                      aphi, atheta, apsi, pars_put, pars_get)

	IMPLICIT NONE
	INCLUDE 'genie_modules.inc'
	EXTERNAL pars_put, pars_get

      INTEGER num_e, num_p, i, j
      REAL phi(num_p), theta(num_p), dist(num_p), energy(num_e+1)
      REAL qx(num_e, num_p), qy(num_e, num_p), qz(num_e, num_p),
     &       eng(num_e, num_p)
      REAL aphi, atheta, apsi 
      REAL tx, ty, tz, nx, ny, sx, sy, sz
      REAL ei, psi, eng_t, rad

      rad = 3.141592654/180.0

      do i=1,num_e 
       eng_t = (energy(i) + energy(i+1))/2.0
        do  j=1,num_p

       tx = (-(sqrt((ei - eng_t)/2.072))*cos(theta(j)*rad)*sin(phi(j)*rad))
       ty = (-(sqrt((ei - eng_t)/2.072))*sin(theta(j)*rad)*sin(phi(j)*rad))
       tz = ((sqrt(ei/2.072))-(sqrt((ei - eng_t)/2.072))*cos(phi(j)*rad))

       sx = (tx*cos(psi*rad) - tz*sin(psi*rad))
       sy = ty
       sz = (tz*cos(psi*rad) + tx*sin(psi*rad))

       nx = (sx*cos(aphi) + sy*sin(aphi))*cos(atheta)-
     $              sz*sin(atheta)
       ny = sy*cos(aphi) - sx*sin(aphi)

       qx(i,j) = nx*cos(-apsi) + ny*sin(-apsi)
       qy(i,j) = ny*cos(-apsi) - nx*sin(-apsi)
       qz(i,j) = sz*cos(atheta) + (sx*cos(aphi) +
     $              sy*sin(aphi))*sin(atheta)

       eng(i,j) = eng_t

       end do 
      end do

	

	CALL MODULE_PUT_ND_ARRAY(

      return
      END



	SUBROUTINE GMODULE_QUERY(PARS_GET, PARS_PUT)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
      EXTERNAL pars_get, pars_put
	
	CHARACTER*20 symbol(1),descrip(1),language(1),type(1),extent(1)

	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR) .EQ. 0) RETURN

	SYMBOL(1) = ''
	DESCRIP(1) = ''
	LANGUAGE(1) = 'FORTRAN'
	TYPE(1) = 'TRANSFORM'
	EXTENT = '*'
                                                              
	CALL MODULE_PUT_STRING_ARRAY(PARS_PUT,'SYMBOL',SYMBOL,1)
	CALL MODULE_PUT_STRING_ARRAY(PARS_PUT,'DESCRIP',DESCRIP,1)
	CALL MODULE_PUT_STRING_ARRAY(PARS_PUT,'LANGUAGE',LANGUAGE,1)
	CALL MODULE_PUT_STRING_ARRAY(PARS_PUT,'TYPE',TYPE,1)
	CALL MODULE_PUT_STRING_ARRAY(PARS_PUT,'EXTENT',EXTENT,1)

	RETURN
	END

