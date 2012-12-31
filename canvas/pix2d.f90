!-------------------------------------------------------------------------------
!
! ROUTINE: PIX2D
!
! WRITTEN: July-Sept 1998
!
! COMMENT:
!
!-------------------------------------------------------------------------------
!
! INPUT:
!
! Arrays which each have size (num_p)
!
!	phi		R spherical polar phi angle of detector pixel
!	theta		R spherical polar theta angle of detector pixel
!	width		R width of detector pixel
!	dist		R sample to pixel distance
!
! Arrays which each have size (num_e+1)
!
!	energy		R  energy bin boundaries
! 
! Other variables
!
!	axes		I  defines what px and py are
! 	psi		R  crystal rotation angle psi
!	ei		R  incident neutron energy 
!       num_p		I  size of the input arrays 
!	num_e		I  size of the input arrays       
!       scale(3)        R scaling factors for vx, vy, vz
!
!-------------------------------------------------------------------------------
!
! OUTPUT:
!
! Arrays which have size (num_e, num_p, 5) 	
!
!	px		R  x pixel centre and vertices
!	py		R  y pixel centre and vertices
!
! Other variables
!
!       info		C  string containing error message
!	flag		I  flag to show sucessful execution of routine
!       aspect          R  aspect ratio for plot
!
!-------------------------------------------------------------------------------

      REAL FUNCTION PIX2D(phi, width, theta, dist, energy, ei, psi, scale, 
     &                    aspect, num_e, num_p, px, py, axes, info, flag)

      INTEGER num_e, num_p, i, j, k, axes(3), flag
      REAL phi(num_p), theta(num_p), dist(num_p), energy(num_e+1)
      REAL width(num_p), px(num_e, num_p, 5), py(num_e, num_p, 5), v(7)
      REAL tx, ty, tz, add_phi(5), add_eng(5)
      REAL ei, psi, rad, new_phi, new_eng, scale(3), aspect
      CHARACTER*50 info

!  Set the routine flag to OK.

      pix2d = 1.0
      flag  = 1

!  Create the arrays for the vertex calculation

      DATA (add_phi(i), i=1,5) / 0, 1, 1, -1,-1 /
      DATA (add_eng(i), i=1,5) / 0, 1,-1, -1, 1 /

      rad    = 3.141592654/180.0

!  CAlculate the vertices

      do k=1, 5
	
       do i=1,num_e 

         new_eng = ((energy(i) + energy(i+1))/2.0) + add_eng(k) *  
     &             ((energy(i+1) - energy(i))/2.0)  
         do  j=1,num_p

       new_phi = phi(j)*rad + add_phi(k)*atan((width(j)/2.0)/dist(j)) 
                                         
       tx = (-(sqrt((ei - new_eng)/2.072))*cos(theta(j)*rad)*sin(new_phi))
       ty = (-(sqrt((ei - new_eng)/2.072))*sin(theta(j)*rad)*sin(new_phi))
       tz = ((sqrt(ei/2.072))-(sqrt((ei - new_eng)/2.072))*cos(new_phi))


!  V is the four dimentional vector for each point

       v(1) = (tx*cos(psi*rad) - tz*sin(psi*rad))/scale(1)
       v(2) = ty/scale(2)
       v(3) = (tz*cos(psi*rad) + tx*sin(psi*rad))/scale(3)
       v(4) = new_eng
       v(5) = sqrt(v(1)**2 + v(2)**2 + v(3)**2) 
       v(6) = 2.072*(v(1)**2 + v(2)**2 + v(3)**2) 
	 v(7) = new_phi/rad

!  Q is the three dimentional vector selected from the 4-vector by axes

       px(i,j,k) = v(axes(1))
       py(i,j,k) = v(axes(2))

        end do 
       end do
      end do

      if ((axes(1) .ge. 4) .or. (axes(2) .ge. 4)) then
       aspect = 1.0
      else  
       aspect = scale(axes(1))/scale(axes(2))
      endif

      return
      END

