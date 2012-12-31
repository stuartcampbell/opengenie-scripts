C-------------------------------------------------------------------------------
C
C ROUTINE: PRO3D
C
C WRITTEN: July-Sept 1998
C
C COMMENT:
C
C-------------------------------------------------------------------------------
C
C INPUT:
C
C Arrays which each have size (num_e, num_p)
C
C 	vx		R	x axis recprocal space coordines of pixels
C 	vy		R	y axis recprocal space coordines of pixels
C 	vz		R	z axis recprocal space coordines of pixels
C 
C Other variables - those which define the plane from which the cut is made
C
C       x(3)		R  three vector of the x axis of the plane 	
C       y(3)		R  three vector of the y axis of the plane 	
C       p(3)		R  three vector to define the position of the plane 	
C       bin(3)		R  the size of the bins
C	axes		I  defines the what qx, qy and qz are 
C       num_p		I  size of the input arrays 
C	num_e		I  size of the input arrays       
C       param		I  determines the plane projection
C       p(3)		R  three vector defining the plane	
C
C-------------------------------------------------------------------------------
C
C OUTPUT:
C
C Arrays which have size (num_e, num_p) 	
C
C 	selectx		R  x index of bin for contributing pixel
C 	selecty		R  y index of bin for contributing pixel
C
C	1)	selectx, set to NULL_DATA for non-contributing pixels
C
C Other variables
C
C	min_x		I  minimum value of selectx
C	min_y		I  minimum value of selecty
C 	size_x		I  extend of selectx
C 	size_y		I  extend of selecty
C       info		C  string containing error message
C	flag		I  flag to show sucessful execution of routine
C       new_p(3)    	R  transformed vector defining the plane	
C
C-------------------------------------------------------------------------------

      REAL FUNCTION PRO3D(p, x, y, new_p, bin, param, axes, vx, vy, vz,
     &                    selectx, selecty, num_e, num_p, min_x, min_y,
     &                    size_x, size_y, info, flag)

      INTEGER errcode,num_e,num_p, i, j, axes(3), param, flag
      INTEGER size_x, size_y, max_x, min_x, max_y, min_y
      REAL selectx(num_e, num_p), selecty(num_e, num_p)
      REAL vx(num_e, num_p), vy(num_e, num_p), vz(num_e, num_p)
      REAL p(3), x(3), y(3), bin(3)
      REAL n(3), norm_x, norm_y, val, const, new_p(3), new_q(3)
      REAL testx_q, testx_e, testy_q, testy_e
      CHARACTER*50 info

C  Set the routines flag to OK.

      pro3d = 1.0
      flag  = 1

C  Check to see if x and y are orthogonal

      if ((x(1)*y(1)+x(2)*y(2)+x(3)*y(3)) .ne. 0.0) then
       info = 'ERROR: x NOT orthogonal to y'
       flag = -1
       return
      endif

C  Need to check whether calculations are possible.....

      testx_q =  x(1)*(axes(1) .le. 3)
     &         + x(2)*(axes(2) .le. 3)
     &         + x(3)*(axes(3) .le. 3)

      testx_e =  x(1)*(axes(1) .eq. 4)
     &         + x(2)*(axes(2) .eq. 4)
     &         + x(3)*(axes(3) .eq. 4)

      testy_q =  y(1)*(axes(1) .le. 3)
     &         + y(2)*(axes(2) .le. 3)
     &         + y(3)*(axes(3) .le. 3)

      testy_e =  y(1)*(axes(1) .eq. 4)
     &         + y(2)*(axes(2) .eq. 4)
     &         + y(3)*(axes(3) .eq. 4)

      if ((((testx_q .eq. 0.0) .or. (testx_e .eq. 0.0)) .eq. 0.0) .or. 
     &    (((testy_q .eq. 0.0) .or. (testy_e .eq. 0.0)) .eq. 0.0)) then
       info = 'ERROR: Incompatable axes for plane=0 cut'
       flag = -1
       return
      endif 

C  Calculate the normalising factors for the x and y vectors
  
      norm_x = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      norm_y = sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))

C  Calculate the cross-product of x and y to get the normal vector n

      write(6,*) norm_x, norm_y
      
      n(1) = (x(2)*y(3) - x(3)*y(2))/(norm_x*norm_y)
      n(2) = (x(3)*y(1) - x(1)*y(3))/(norm_x*norm_y)
      n(3) = (x(1)*y(2) - x(2)*y(1))/(norm_x*norm_y)

      write(6,*) n(1), n(2), n(3)

C  If we can cut in the x y param without projection then we should transform
C  p into the x y n coordinate system

      new_p(1) = (p(1)*x(1)+p(2)*x(2)+p(3)*x(3))/norm_x
      new_p(2) = (p(1)*y(1)+p(2)*y(2)+p(3)*y(3))/norm_y
      new_p(3) = (p(1)*n(1)+p(2)*n(2)+p(3)*n(3))

      write(6,*) new_p(1), new_p(2), new_p(3)

C  Calculate the constant value n.p = C

      const = n(1)*p(1)+n(2)*p(2)+n(3)*p(3)

C  const = new_p(3)!!!!!

C  We will also want the max and min for the binned x and y values - so 
C  get them now as we are performing loops 

      min_x =  1e4
      max_x = -1e4
      min_y =  1e4
      max_y = -1e4

      if (param .eq. 3.0) then
       info = 'INFO: The yz is the projected plane' 
       new_p(1) = p(2)
       new_p(2) = p(3)
       new_p(3) = const
       if (n(1) .ne. 0.0) then
        do i=1, num_e
         do j=1, num_p
          val = (const - n(2)*vy(i,j) - n(3)*vz(i,j))/n(1)
          if ((abs(vx(i,j)-val)) .le. (bin(3)/2.0)) then
           selectx(i,j) = nint((vy(i,j) - p(2))/bin(1))
           selecty(i,j) = nint((vz(i,j) - p(3))/bin(2))
           if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
           if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
           if (selecty(i,j) .lt. min_y) min_y = selecty(i,j)
           if (selecty(i,j) .gt. max_y) max_y = selecty(i,j)
          else 
           selectx(i,j) = -1e30
           selecty(i,j) = -1e30
          endif
         enddo
        enddo
       else
        info = 'ERROR: The yz plane is not projectable'
        flag = -1
        return
       endif
      endif

      if  (param .eq. 2) then
       info = 'INFO: The xz is the projected plane' 
       new_p(1) = p(1)
       new_p(2) = p(3)
       new_p(3) = const
       if (n(2) .ne. 0.0) then
        do i=1, num_e
         do j=1, num_p
          val = (const - n(1)*vx(i,j) - n(3)*vz(i,j))/n(2)
          if ((abs(vy(i,j)-val)) .le. (bin(3)/2.0)) then
           selectx(i,j) = nint((vx(i,j) - p(1))/bin(1))
           selecty(i,j) = nint((vz(i,j) - p(3))/bin(2))
           if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
           if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
           if (selecty(i,j) .lt. min_y) min_y = selecty(i,j)
           if (selecty(i,j) .gt. max_y) max_y = selecty(i,j)
          else 
           selectx(i,j) = -1e30
           selecty(i,j) = -1e30
          endif
         enddo
        enddo
       else
        info = 'ERROR: The xz plane is not projectable'
        flag = -1
        return
       endif
      endif


      if  (param .eq. 1) then
       info = 'INFO: The xy is the projected plane' 
       new_p(1) = p(1)
       new_p(2) = p(2)
       new_p(3) = const
       if (n(3) .ne. 0.0) then
        do i=1, num_e
         do j=1, num_p
          val = (const - n(2)*vy(i,j) - n(1)*vx(i,j))/n(3)
          if ((abs(vz(i,j)-val)) .le. (bin(3)/2.0)) then
           selectx(i,j) = nint((vx(i,j) - p(1))/bin(1))
           selecty(i,j) = nint((vy(i,j) - p(2))/bin(2))
           if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
           if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
           if (selecty(i,j) .lt. min_y) min_y = selecty(i,j)
           if (selecty(i,j) .gt. max_y) max_y = selecty(i,j)
         else 
           selectx(i,j) = -1e30
           selecty(i,j) = -1e30
          endif
         enddo
        enddo
       else
        info = 'ERROR: The xy plane is not projectable'
        flag = -1
        return
       endif
      endif

      if  (param .eq. 0) then
       info = 'INFO: The actual plane is the projected plane' 
       do i=1, num_e
        do j=1, num_p
	 new_q(3) = n(1)*vx(i,j)+n(2)*vy(i,j)+n(3)*vz(i,j)
         if ((abs(new_p(3)-new_q(3))) .le. (bin(3)/2.0)) then
          new_q(1) = (x(1)*vx(i,j)+x(2)*vy(i,j)+x(3)*vz(i,j))/norm_x
          new_q(2) = (y(1)*vx(i,j)+y(2)*vy(i,j)+y(3)*vz(i,j))/norm_y
          selectx(i,j) = nint((new_q(1) - new_p(1))/bin(1))
          selecty(i,j) = nint((new_q(2) - new_p(2))/bin(2))
           if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
           if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
           if (selecty(i,j) .lt. min_y) min_y = selecty(i,j)
           if (selecty(i,j) .gt. max_y) max_y = selecty(i,j)
         else 
          selectx(i,j) = -1e30
          selecty(i,j) = -1e30
         endif
        enddo
       enddo
      endif

      size_x = max_x - min_x + 1
      size_y = max_y - min_y + 1

      return
      END
       
	
