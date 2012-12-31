C-------------------------------------------------------------------------------
C
C ROUTINE: CUT3D
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
C	axes		I  defines the what vx, vy and vz are 
C       num_p		I  size of the input arrays 
C	num_e		I  size of the input arrays       
C       param		I  determines the plane projection
C
C Other variables - those which define the cut in the plane
C
C       q(3)		R  three vector of the start of the cut axis 	
C       r(3)		R  three vector to define direction of the cut axis 	
C       cut_bin		R  the size of the bins
C       width		R  the width of the cut
C
C-------------------------------------------------------------------------------
C
C OUTPUT:
C
C Arrays which have size (size_x, size_y) 	
C
C 	selectx		R  x index of bin for contributing pixel
C
C	1)	selectx, set to NULL_DATA for non-contributing pixels
C
C Other variables
C
C	min_x		I  minimum value of selectx
C 	size_x		I  extend of selectx
C       info		C  string containing error message
C	flag		I  flag to show sucessful execution of routine
C       new_p(3)    	R  transformed vector of the start of the cut axis 	
C
C-------------------------------------------------------------------------------

      REAL FUNCTION CUT3D(p, x, y, new_p, bin, param, axes, vx, vy, vz,
     &                        selectx, num_e, num_p, min_x, 
     &                        size_x, q, r, cut_bin, width, info, flag)


      INTEGER errcode, num_e, num_p, i, j, axes(3)
      INTEGER size_x, max_x, min_x , param, flag
      REAL vx(num_e, num_p), vy(num_e, num_p), vz(num_e, num_p)
      REAL selectx(num_e, num_p), selecty(num_e, num_p)
      REAL p(3), q(2), r(2), x(3), y(3), bin(3), width, cut_bin
      REAL new_x(3), new_y(3)
      REAL n(3), norm_x, norm_y, val, const, new_p(3), new_q(3), new_r(3)
      REAL testx_q, testy_q
      CHARACTER*50 info

C  Set routine flag to OK.

      cut3d = 1.0
      flag  = 1

C  Check to see if x and y are orthogonal

      if ((x(1)*y(1)+x(2)*y(2)+x(3)*y(3)) .ne. 0.0) then
       info = 'ERROR - x NOT orthogonal to y'
       flag = -1
       return
      endif

C  Need to check whether general cut is possible i.e a q vs q plane.....

      testx_q =  x(1)*(axes(1) .le. 3)
     &         + x(2)*(axes(2) .le. 3)
     &         + x(3)*(axes(3) .le. 3)

      testy_q =  y(1)*(axes(1) .le. 3)
     &         + y(2)*(axes(2) .le. 3)
     &         + y(3)*(axes(3) .le. 3)

      if ((q(1) .ne. r(1)) .and. (q(2) .ne. r(2))) then 
       if ((testx_q .ne. 0.0) .and. (testy_q .ne. 0.0)) then
        info = 'INFO: General Cut'
       else
        info = 'ERROR: Incompatable Axes for General Cut'
        flag = -1
        return
       endif 
      endif

      max_x = -1e4
      min_x =  1e4

C  First see if the cut is carried out without projection

      if  (param .eq. 0) then

      do i=1,3
        new_x(i) = (r(1)-q(1))*x(i) + (r(2)-q(2))*y(i)
        new_y(i) = (r(2)-q(2))*x(i) - (r(1)-q(1))*y(i)
        new_q(i) = q(1)*x(i) + q(2)*y(i) + p(i)
      enddo

      write(6,*) new_x(1), new_x(2), new_x(3)
      write(6,*) new_y(1), new_y(2), new_y(3)
      write(6,*) new_q(1), new_q(2), new_q(3)
       
C  Calculate the normalising factors for the x and y vectors
  
      norm_x = sqrt(new_x(1)*new_x(1)+new_x(2)*new_x(2)+new_x(3)*new_x(3))
      norm_y = sqrt(new_y(1)*new_y(1)+new_y(2)*new_y(2)+new_y(3)*new_y(3))

C  Calculate the cross-product of x and y to get the normal vector n

      n(1) = (new_x(2)*new_y(3) - new_x(3)*new_y(2))/(norm_x*norm_y)
      n(2) = (new_x(3)*new_y(1) - new_x(1)*new_y(3))/(norm_x*norm_y)
      n(3) = (new_x(1)*new_y(2) - new_x(2)*new_y(1))/(norm_x*norm_y)

C  If we can cut in the x y plane without projection then we should transform
C  p into the x y n coordinate system

      new_p(1) = (new_q(1)*new_x(1)+new_q(2)*new_x(2)+new_q(3)*new_x(3))/norm_x
      new_p(2) = (new_q(1)*new_y(1)+new_q(2)*new_y(2)+new_q(3)*new_y(3))/norm_y
      new_p(3) = (new_q(1)*n(1)+new_q(2)*n(2)+new_q(3)*n(3))

      info = 'INFO: The actual plane is the projected plane' 
      do j=1, num_p
       do i=1, num_e
	new_r(3) = n(1)*vx(i,j)+n(2)*vy(i,j)+n(3)*vz(i,j)
        new_r(2) =(new_y(1)*vx(i,j)+new_y(2)*vy(i,j)+new_y(3)*vz(i,j))/norm_y

        if (((abs(new_p(3)-new_r(3))) .le. (bin(3)/2.0)) .and. 
     &      ((abs(new_p(2)-new_r(2))) .le. (width/2.0))) then
         new_r(1) =(new_x(1)*vx(i,j)+new_x(2)*vy(i,j)+new_x(3)*vz(i,j))/norm_x
         selectx(i,j) = nint((new_r(1) - new_p(1))/cut_bin)
         if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
         if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
        else 
          selectx(i,j) = -1e30
        endif
        enddo
       enddo
      endif

      if  (param .ne. 0) then

C  Calculate the cross-product of x and y to get the normal vector n

      norm_x = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      norm_y = sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))

      n(1) = (x(2)*y(3) - x(3)*y(2))/(norm_x*norm_y)
      n(2) = (x(3)*y(1) - x(1)*y(3))/(norm_x*norm_y)
      n(3) = (x(1)*y(2) - x(2)*y(1))/(norm_x*norm_y)

      const = n(1)*p(1)+n(2)*p(2)+n(3)*p(3)

      new_x(1) = r(1) - q(1)
      new_x(2) = r(2) - q(2)
      new_y(1) = r(2) - q(2)
      new_y(2) = q(1) - r(1)

      write(6,*) new_x(1), new_x(2)
      write(6,*) new_y(1), new_y(2)

      norm_x = sqrt(new_x(1)*new_x(1)+new_x(2)*new_x(2))
      norm_y = sqrt(new_y(1)*new_y(1)+new_y(2)*new_y(2))

      new_p(1) = q(1)*new_x(1)+q(2)*new_x(2)
      new_p(2) = q(1)*new_y(1)+q(2)*new_y(2)
      new_p(3) = 0.0

      if  (param .eq. 3) then
       info = 'INFO: The yz plane is the cut plane'
       if (n(1) .ne. 0.0) then
        do j=1, num_p
         do i=1, num_e
          val = (const - n(2)*vy(i,j) - n(3)*vz(i,j))/n(1)
          new_r(2) = (new_y(1)*vy(i,j)+new_y(2)*vz(i,j))/norm_y
        if (((abs(vx(i,j)-val)) .le. (bin(3)/2.0)) .and. 
     &      ((abs(new_p(2)-new_r(2))) .le. (width/2.0))) then
          new_r(1) = (new_x(1)*vy(i,j)+new_x(2)*vz(i,j))/norm_x
          selectx(i,j) = nint((new_r(1) - new_p(1))/cut_bin)
           if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
           if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
          else 
           selectx(i,j) = -1e30
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
       info = 'INFO: The xz plane is the cut plane'
       if (n(2) .ne. 0.0) then
        do j=1, num_p
         do i=1, num_e
          val = (const - n(1)*vx(i,j) - n(3)*vz(i,j))/n(2)
          new_r(2) = (new_y(1)*vx(i,j)+new_y(2)*vz(i,j))/norm_y
        if (((abs(vy(i,j)-val)) .le. (bin(3)/2.0)) .and. 
     &      ((abs(new_p(2)-new_r(2))) .le. (width/2.0))) then
          new_r(1) = (new_x(1)*vx(i,j)+new_x(2)*vz(i,j))/norm_x
          selectx(i,j) = nint((new_r(1) - new_p(1))/cut_bin)
           if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
           if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
          else 
           selectx(i,j) = -1e30
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
       info = 'INFO: The xy plane is the cut plane'
       if (n(3) .ne. 0.0) then
        do j=1, num_p
         do i=1, num_e
          val = (const - n(1)*vy(i,j) - n(3)*vx(i,j))/n(3)
          new_r(2) = (new_y(1)*vy(i,j)+new_y(2)*vx(i,j))/norm_y
        if (((abs(vz(i,j)-val)) .le. (bin(3)/2.0)) .and. 
     &      ((abs(new_p(2)-new_r(2))) .le. (width/2.0))) then
          new_r(1) = (new_x(1)*vy(i,j)+new_x(2)*vx(i,j))/norm_x
          selectx(i,j) = nint((new_r(1) - new_p(1))/cut_bin)
           if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
           if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
          else 
           selectx(i,j) = -1e30
          endif
         enddo
        enddo
       else
        info = 'ERROR: The xy plane is not projectable'
        flag = -1
       endif
      endif

      endif

      write(6,*) min_x, max_x

      size_x = max_x - min_x + 1

      return
 
      end
