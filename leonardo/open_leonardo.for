*******************************************************************
*                       OPEN_Leonardo v1.0
*
*                    S.I.Campbell & C.D.Frost
*                         ISIS Facility
*
* Developed from IDL version of leonardo written by C.D.Frost
*
*******************************************************************
* CONTAINS	:	GETRUN, GETPAR, CALCULATE 
*******************************************************************
        
*******************************************************************
*
* GETRUN :	reads in .SPE files produced by HOMER
*
*******************************************************************


	SUBROUTINE GETRUN(pars_get, pars_put)
	
	IMPLICIT NONE
	INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put
	
	INCLUDE 'leonardo.inc'
	REAL e_points(e_dimension+1), p_points(p_dimension+1)
	REAL data(p_dimension,e_dimension),
     &	     error(p_dimension,e_dimension)
	INTEGER unit, stat
	CHARACTER*80 input_file
	
	INTEGER num_p_points, num_e_points
	
	
	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR) .EQ. 0) RETURN
	
	input_file=' '

	CALL MODULE_GET_STRING(pars_get, 'file', input_file)

	CALL MODULE_INFORMATION('Executing GETRUN')
	
	unit = 34
	stat = 0
		
	OPEN(unit, FILE=input_file, STATUS='OLD', IOSTAT=stat)
	IF(stat.ne.0) THEN
		CALL module_error('GETRUN','Unable to open ' 
     &           //input_file,'Check file exists and is readable')
		RETURN
	ENDIF

* Read data from file

	READ(unit,*) num_p_points, num_e_points          
		
	CALL READ_DATA(num_e_points,num_p_points,e_points,p_points,
     &		data, error, unit,stat,pars_put)	
		                   
	RETURN
	END
	

	SUBROUTINE READ_DATA(num_e,num_p,e_points,p_points,
     &		data,error,unit,stat,pars_put)
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

	dims_array(1) = num_p
	dims_array(2) = num_e
	                            	
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'error',error,
     &	dims_array,2)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'values',data,
     &	dims_array,2)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'energy',e_points,num_e+1)
	CALL MODULE_PUT_REAL_ARRAY(pars_put,'phi',p_points,num_p+1)
	CALL MODULE_PUT_INT(pars_put,'num_e',num_e)
	CALL MODULE_PUT_INT(pars_put,'num_p',num_p)

	close(unit,iostat=stat)
	
C	CALL PGENV(0., 1.0 , 0., 1.0 , 1, -2)
C	CALL FREDDY(data,num_p,num_e,1.0,25.0)

100	FORMAT(8G10.4)

	RETURN
	END
	

*******************************************************************
*
* GETPAR :	reads in the parameter file that is associated with the 
*		data.	
*
*******************************************************************

	SUBROUTINE GETPAR(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put

	INCLUDE 'leonardo.inc'

	INTEGER i, j, unit, stat, num_p, number, dims_array(2)
	REAL det(p_dimension,5)
	CHARACTER*80 input_file

	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR).EQ.0) RETURN

	input_file=' '

	CALL MODULE_GET_STRING(pars_get, 'par_file', input_file)
	CALL MODULE_GET_INT(pars_get, 'num_p', num_p)

	CALL MODULE_INFORMATION('Executing GETPAR')
	
	unit = 34
	stat = 0
	
	OPEN(unit, FILE=input_file, STATUS='OLD', IOSTAT=stat)
	IF(stat.ne.0) THEN
		CALL module_error('GETPAR','Unable to open ' 
     &           //input_file,'Check file exists and is readable')
		RETURN
	ENDIF
	
	READ(unit,*) number
	
      IF(number.ne.num_p) THEN
		CALL module_error('GETPAR','PAR file does not match data',
     &	'Check file name is correct')
		RETURN
	ENDIF
	
	CALL READ_PAR(num_p,det,unit,pars_put)
	
	close(unit,iostat=stat)

	RETURN
	END
	
                            
	SUBROUTINE READ_PAR(n, det, unit, pars_put)
	IMPLICIT NONE
	INCLUDE 'genie_modules.inc'
	EXTERNAL pars_put
	INTEGER n, i, j, dims_array(2), unit
	REAL det(n,5)
	DO i = 1, n
		READ(unit,*) (det(i,j),j=1,5)
	ENDDO

	dims_array(1) = n
	dims_array(2) = 5

	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'det_array',
     &				det,dims_array,2)

	RETURN
	END
                                                                   
*******************************************************************
*
* CALCULATE : 
*
*******************************************************************

	SUBROUTINE CALCULATE(pars_get, pars_put)

	IMPLICIT NONE
	INCLUDE 'genie_modules.inc'
	INCLUDE 'leonardo.inc'
	INCLUDE 'constant.inc'
	EXTERNAL pars_get, pars_put
	
	INTEGER num_e, num_p
	REAL det_phi(p_dimension), det_theta(p_dimension),
     &	     det_width(p_dimension),det_height(p_dimension),
     &	     det_l2(p_dimension),energy(e_dimension+1)
	REAL epsi(5,e_dimension), phi(5,p_dimension), 
     &       theta(5,p_dimension), add_matrix(5,3)
	REAL qx(5,e_dimension,p_dimension),
     &	     qy(5,e_dimension,p_dimension),
     &       qz(5,e_dimension,p_dimension),
     &       eng(5,e_dimension,p_dimension)
	REAL c_x(p_dimension*e_dimension,4)
	REAL c_y(p_dimension*e_dimension,4)
	REAL centre_x(p_dimension,e_dimension)
	REAL centre_y(p_dimension,e_dimension)

C	integer*4 time, c1, c2, elapse
C	real*4 dtime, delta, tarray(2)
C 	external dtime, time
	
	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR).EQ.0) RETURN
	
C	delta = dtime ( tarray )
C	write(6,*) ' Start time',delta
	
C	c1 = time()
C	write(6,*) ' Start time',c1
	
	
	add_matrix(1,1) = 0.0 
	add_matrix(1,2) = 0.0
 	add_matrix(1,3) = 0.0
 	
	add_matrix(2,1) = 1.0 
	add_matrix(2,2) = 1.0
 	add_matrix(2,3) = 1.0
	
	add_matrix(3,1) =-1.0 
	add_matrix(3,2) = 1.0
 	add_matrix(3,3) = 1.0
	
	add_matrix(4,1) =-1.0 
	add_matrix(4,2) =-1.0
 	add_matrix(4,3) = 1.0

	add_matrix(5,1) = 1.0 
	add_matrix(5,2) =-1.0
 	add_matrix(5,3) = 1.0

	CALL MODULE_GET_INT(pars_get,'num_e',num_e)
	CALL MODULE_GET_INT(pars_get,'num_p',num_p)
                                               	
C        delta = dtime ( tarray )
C	write(6,*) ' Before calling CALC_SUB',delta    
        
C        c2 = time()
C       elapse = c2-c1
C        write(6,*) ' Elasped time so far ...',elapse
                                       
	CALL CALC_SUB(num_e,num_p,add_matrix,pars_get,pars_put,
     &	det_l2,det_phi,det_theta,det_height,det_width,energy,
     &	epsi,phi,theta,qx,qy,qz,eng,c_x,c_y,centre_x,centre_y)
C	,c1)

C	delta = dtime ( tarray )
C	write(6,*) ' Returning to GENIE',delta	
C        c2 = time()
C        elapse = c2-c1
C        write(6,*) ' Total time ...',elapse
		
	RETURN
	END


	SUBROUTINE CALC_SUB(num_e,num_p,add_matrix,pars_get,pars_put,
     &	det_l2,det_phi,det_theta,det_height,det_width,energy,
     &	epsi,phi,theta,qx,qy,qz,eng,c_x,c_y,centre_x,centre_y)
C	,c1)

	IMPLICIT NONE
	INCLUDE 'genie_modules.inc'
	INCLUDE 'leonardo.inc'
	INCLUDE 'constant.inc'
	EXTERNAL pars_get, pars_put
	
	
	INTEGER i, j, k, plot_mode, dims_array(2), index, nf
	REAL ei, psi, vtheta, apsi, vphi, factors(3)
	INTEGER num_e, num_p, total
	REAL energy1, energy2, delta_phi, delta_theta
	REAL det_phi(num_p), det_theta(num_p),
     &	     det_width(num_p), det_height(num_p),
     &	     det_l2(num_p), energy(num_e+1)
	REAL epsi(5,num_e), phi(5,num_p), 
     &       theta(5,num_p), add_matrix(5,3)
	REAL qx(5,num_e,num_p),
     &	     qy(5,num_e,num_p),
     &       qz(5,num_e,num_p),
     &       eng(5,num_e,num_p)
	REAL c_x(num_p*num_e,4)
	REAL c_y(num_p*num_e,4)
	REAL centre_x(num_p,num_e), centre_y(num_p,num_e)	
	
C	integer*4 time, c1, c2, elapse
C	real*4 dtime, delta, tarray(2)
C  	external dtime, time	
  		
	CALL MODULE_GET_REAL(pars_get,'ei',ei)
	CALL MODULE_GET_REAL(pars_get,'psi',psi)
	CALL MODULE_GET_REAL(pars_get,'vtheta',vtheta)
	CALL MODULE_GET_REAL(pars_get,'apsi',apsi)
	CALL MODULE_GET_REAL(pars_get,'vphi',vphi)
	CALL MODULE_GET_INT(pars_get,'plot_mode',plot_mode)
	nf = 3
	CALL MODULE_GET_REAL_ARRAY(pars_get,'factors',factors,nf)

	CALL MODULE_GET_REAL_ARRAY(pars_get,'energy',energy,
     &	num_e+1)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'det_phi',det_phi,
     &	num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'det_theta',det_theta,
     &	num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'det_height',
     &	det_height,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'det_width',
     &	det_width,num_p)
	CALL MODULE_GET_REAL_ARRAY(pars_get,'det_l2',det_l2,num_p)


	DO i = 1, num_e

	 energy1 = (energy(i)+energy(i+1))/2.0
	 energy2 = (energy(i)-energy(i+1))/2.0
	
	  epsi(1,i) = energy1+add_matrix(1,1)*energy2	
	  epsi(2,i) = energy1+add_matrix(2,1)*energy2
	  epsi(3,i) = energy1+add_matrix(3,1)*energy2
	  epsi(4,i) = energy1+add_matrix(4,1)*energy2
	  epsi(5,i) = energy1+add_matrix(5,1)*energy2
	
	ENDDO

	DO i = 1, num_p

	 delta_theta = (180.0/(2*pi))*asin(det_height(i)/det_l2(i))
	 delta_phi = (180.0/(2*pi))*asin(det_width(i)/det_l2(i))


	  phi(1,i) = det_phi(i)+add_matrix(1,2)*delta_phi	
	  phi(2,i) = det_phi(i)+add_matrix(2,2)*delta_phi	
	  phi(3,i) = det_phi(i)+add_matrix(3,2)*delta_phi	
	  phi(4,i) = det_phi(i)+add_matrix(4,2)*delta_phi	
	  phi(5,i) = det_phi(i)+add_matrix(5,2)*delta_phi
	  	
	  theta(1,i) = det_theta(i)+add_matrix(1,3)*delta_theta	
	  theta(2,i) = det_theta(i)+add_matrix(2,3)*delta_theta	
	  theta(3,i) = det_theta(i)+add_matrix(3,3)*delta_theta	
	  theta(4,i) = det_theta(i)+add_matrix(4,3)*delta_theta	
	  theta(5,i) = det_theta(i)+add_matrix(5,3)*delta_theta	
  	  
	ENDDO

	CALL CALCQ(qx,qy,qz,eng,epsi,phi,theta,ei,psi,vtheta,apsi,
     &			vphi,factors,num_e,num_p)

	IF(plot_mode.eq.9) THEN
 	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = eng(k,i,j)
 	    c_y(index,k-1) = phi(k,j)
 	    centre_x(j,i) = eng(1,i,j)
	    centre_y(j,i) = phi(1,j)
     	   ENDDO
	  ENDDO
	
C	delta = dtime ( tarray )
C	write(6,*) ' Time at loop ',k,' is',delta
C        c2 = time()
C        elapse = c2-c1
C        write(6,*)' Elasped time so far (IN CALCULATE) ...',elapse
	  
	 ENDDO
	ELSE IF (plot_mode.eq.8) THEN
 	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = 2.072*((qx(k,i,j)*qx(k,i,j))+(qy(k,i,j)*
     &				   qy(k,i,j))+(qz(k,i,j)*qz(k,i,j)))
 	    c_y(index,k-1) = eng(k,i,j)
	    centre_x(j,i) = 2.072*((qx(1,i,j)*qx(1,i,j))+(qy(1,i,j)*
     &				   qy(1,i,j))+(qz(1,i,j)*qz(1,i,j)))
	    centre_y(j,i) = eng(1,i,j)
	   ENDDO
	  ENDDO
	 ENDDO
	ELSE IF (plot_mode.eq.7) THEN
	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = sqrt((qx(k,i,j)*qx(k,i,j))+(qy(k,i,j)*
     &				   qy(k,i,j))+(qz(k,i,j)*qz(k,i,j)))
 	    c_y(index,k-1) = eng(k,i,j)
 	    centre_x(j,i) = sqrt((qx(1,i,j)*qx(1,i,j))+(qy(1,i,j)*
     &				   qy(1,i,j))+(qz(1,i,j)*qz(1,i,j)))
	    centre_y(j,i) = eng(1,i,j)
	   ENDDO
	  ENDDO
	 ENDDO	
	ELSE IF (plot_mode.eq.6) THEN
	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = qz(k,i,j)
 	    c_y(index,k-1) = eng(k,i,j)
	    centre_x(j,i) = qz(1,i,j)
	    centre_y(j,i) = eng(1,i,j)
	   ENDDO
	  ENDDO
	 ENDDO		
	ELSE IF (plot_mode.eq.5) THEN
	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = qy(k,i,j)
 	    c_y(index,k-1) = eng(k,i,j)
 	    centre_x(j,i) = qy(1,i,j)
	    centre_y(j,i) = eng(1,i,j)
	   ENDDO
	  ENDDO
	 ENDDO		
	ELSE IF (plot_mode.eq.4) THEN
	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = qx(k,i,j)
 	    c_y(index,k-1) = eng(k,i,j)
	    centre_x(j,i) = qx(1,i,j)
	    centre_y(j,i) = eng(1,i,j)
	   ENDDO
	  ENDDO
	 ENDDO		
	ELSE IF (plot_mode.eq.3) THEN
 	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = qx(k,i,j)
 	    c_y(index,k-1) = qy(k,i,j)
 	    centre_x(j,i) = qx(1,i,j)
	    centre_y(j,i) = qy(1,i,j)
	   ENDDO
	  ENDDO
	 ENDDO
	ELSE IF (plot_mode.eq.2) THEN
 	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = qz(k,i,j)
 	    c_y(index,k-1) = qy(k,i,j)
	    centre_x(j,i) = qz(1,i,j)
	    centre_y(j,i) = qy(1,i,j)
	   ENDDO
	  ENDDO
	 ENDDO
	ELSE IF (plot_mode.eq.1) THEN
 	 DO k = 2, 5
	  DO i = 1,num_e
	   DO j = 1,num_p
	    index = ((i-1)*num_p)+j
	    c_x(index,k-1) = qz(k,i,j)
 	    c_y(index,k-1) = qx(k,i,j)
 	    centre_x(j,i) = qz(1,i,j)
	    centre_y(j,i) = qx(1,i,j)
	   ENDDO
	  ENDDO
	 ENDDO
	ENDIF	

	total = num_e*num_p
	
	dims_array(1) = total
	dims_array(2) = 4

C	delta = dtime ( tarray )
C	write(6,*) ' Before putting all variables back',delta
C        c2 = time()
C        elapse = c2-c1
C        write(6,*)' Elasped time so far (IN CALCULATE) ...',elapse	
        
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'c_x',c_x,
     &			dims_array,2)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'c_y',c_y,
     &			dims_array,2)

	dims_array(1) = num_p
	dims_array(2) = num_e
	
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'centre_x',centre_x,
     &		dims_array,2)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'centre_y',centre_y,
     &		dims_array,2)
	
C	delta = dtime ( tarray )
C	write(6,*) ' After putting variables back',delta
C	       c2 = time()
C        elapse = c2-c1
C        write(6,*)' Elasped time so far (IN CALCULATE) ...',elapse
	
	RETURN
	END


*******************************************************************
*
* CALCQ : 
*
*******************************************************************

	SUBROUTINE CALCQ(qx,qy,qz,eng,energy,phi,theta,ei,psi,
     &			vtheta,apsi,vphi,factors,num_e,num_p)	
	
	IMPLICIT NONE
	INCLUDE 'constant.inc'

	INTEGER i, j, k, num_e, num_p
	REAL ei, psi, energy(5,num_e)
	REAL phi(5,num_p), theta(5,num_p)
	REAL tx, ty, tz, sx, sy, sz, nx, ny
	REAL vphi, vtheta, apsi, factors(3)
	REAL qx(5,num_e,num_p),
     &	     qy(5,num_e,num_p),
     &       qz(5,num_e,num_p),
     &	     eng(5,num_e,num_p)
	                                    	

	DO k = 1,5
C	DO k = 2,5
	
       	 DO i = 1, num_e
	  DO j = 1, num_p
	tx = -(SQRT((ei-energy(k,i))/2.072))*COS(theta(k,j)*dtr) 
     &		*sin(phi(k,j)*dtr)
	ty = -(SQRT((ei-energy(k,i))/2.072))*SIN(theta(k,j)*dtr)
     &	*sin(phi(k,j)*dtr)
	tz = (SQRT(ei/2.072))-(sqrt((ei-energy(k,i))/2.072))
     &	*cos(phi(k,j)*dtr)

* now transform the data to account for rotations of the crystal
* with respect to the neutron beam.
*
*OLD*	qx(j,i) = tx*COS(psi*dtr)+tz*SIN(psi*dtr)
*OLD*	qy(j,i) = ty
*OLD*	qz(j,i) = tz*COS(psi*dtr)-tx*SIN(psi*dtr)

	sx = (tx*cos(psi*dtr) - tz*sin(psi*dtr))
	sy = ty
	sz = (tz*cos(psi*dtr) + tx*sin(psi*dtr))

	nx = (sx*cos(vphi*dtr)+sy*sin(vphi*dtr))*cos(vtheta*dtr)-
     &	sz*sin(vtheta*dtr)
	ny = sy*cos(vphi*dtr)-sx*sin(vphi*dtr)

	qx(k,i,j) = nx*cos(-apsi*dtr) + ny*sin(-apsi*dtr)
	qy(k,i,j) = ny*cos(-apsi*dtr) - nx*sin(-apsi*dtr)
	qz(k,i,j) = sz*cos(vtheta*dtr) + (sx*cos(vphi*dtr) +
     &		sy*sin(vphi*dtr))*sin(vtheta*dtr)

	qx(k,i,j) = qx(k,i,j) / factors(1)
	qy(k,i,j) = qy(k,i,j) / factors(2)
	qz(k,i,j) = qz(k,i,j) / factors(3)

        eng(k,i,j) = energy(k,i)
	
	  ENDDO
	 ENDDO
	ENDDO
	                     	
	RETURN
	END

*******************************************************************

	SUBROUTINE GMODULE_QUERY(pars_get, pars_put)
	IMPLICIT NONE
      INCLUDE 'genie_modules.inc'
	EXTERNAL pars_get, pars_put

	INTEGER n
	PARAMETER (n=3)  
	CHARACTER*20 symbol(n),descrip(n),language(n),type(n),extent(n)
	
	IF(MODULE_VERSION_OK(GENIE_MAJOR, GENIE_MINOR) .EQ. 0) RETURN

	symbol(1)  = 'getrun'
	descrip(1) = 'Gets .SPE files'
	language(1) = 'FORTRAN'
	type(1) = 'INPUT' 
	EXTENT(1) = 'SPE'

	symbol(2)  = 'getpar'
	descrip(2) = 'Gets .PAR files'
	language(2) = 'FORTRAN'
	type(2) = 'INPUT' 
	EXTENT(2) = 'PAR'

	symbol(3)  = 'calculate'
	descrip(3) = 'Prepares data'
	language(3) = 'FORTRAN'
	type(3) = 'TRANSFORM' 
	EXTENT(3) = '*'

	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'symbol', symbol, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'descrip', descrip, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'language', language, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'type', type, n)
	CALL MODULE_PUT_STRING_ARRAY(pars_put, 'extent', extent, n)
	RETURN	
	END


