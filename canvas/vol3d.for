C-------------------------------------------------------------------------------
C
C ROUTINE: VOL3D
C
C WRITTEN: July-Sept 1998
C
C COMMENT:
C
C-------------------------------------------------------------------------------
C
C INPUT:
C
C Arrays which each have size (num_p)
C
C	phi		R spherical polar phi angle of detector pixel
C	theta		R spherical polar theta angle of detector pixel
C	dist		R sample to pixel distance
C
C Arrays which each have size (num_e+1)
C
C	energy		R  energy bin boundaries
C 
C Other variables
C
C	axes		I  defines what px and py are
C 	psi		R  crystal rotation angle psi
C	ei		R  incident neutron energy 
C       num_p		I  size of the input arrays 
C	num_e		I  size of the input arrays       
C
C-------------------------------------------------------------------------------
C
C OUTPUT:
C
C Arrays which have size (num_e, num_p) 	
C
C	vx		R  x pixel in reciprocal space coordinates
C	vy		R  y pixel in reciprocal space coordinates
C	vz		R  z pixel in reciprocal space coordinates
C
C Other variables
C
C       info		C  string containing error message
C	flag		I  flag to show sucessful execution of routine
C
C-------------------------------------------------------------------------------

      REAL FUNCTION VOL3D(phi, theta, dist, energy, ei, psi, 
     &                    num_e, num_p, vx, vy, vz, axes, info, flag)

      INTEGER num_e, num_p, i, j, axes(3), flag
      REAL phi(num_p), theta(num_p), dist(num_p), energy(num_e+1)
      REAL vx(num_e, num_p), vy(num_e, num_p), vz(num_e, num_p), v(4)
      REAL tx, ty, tz
      REAL ei, psi, rad, eng_t
      CHARACTER*50 info

C  Set the routine flag to OK.

      vol3d = 1.0
      flag  = 1

      rad = 3.141592654/180.0

      do i=1,num_e 
        eng_t = (energy(i) + energy(i+1))/2.0
        do  j=1,num_p

       tx = (-(sqrt((ei - eng_t)/2.072))*cos(theta(j)*rad)*sin(phi(j)*rad))
       ty = (-(sqrt((ei - eng_t)/2.072))*sin(theta(j)*rad)*sin(phi(j)*rad))
       tz = ((sqrt(ei/2.072))-(sqrt((ei - eng_t)/2.072))*cos(phi(j)*rad))

C  V is the four dimentional vector for each point

       v(1) = (tx*cos(psi*rad) - tz*sin(psi*rad))
       v(2) = ty
       v(3) = (tz*cos(psi*rad) + tx*sin(psi*rad))
       v(4) = energy(i) 

C  Q is the three dimentional vector selected from the 4-vector by axes

       vx(i,j) = v(axes(1))
       vy(i,j) = v(axes(2))
       vz(i,j) = v(axes(3))

       end do 
      end do

      return
      end
