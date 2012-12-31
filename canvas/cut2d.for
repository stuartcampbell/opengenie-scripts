C-------------------------------------------------------------------------------
C
C ROUTINE: CUT2D
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
C
C	1)	here x refers to the cut axis
C	2)	selectx, set to NULL_DATA for non-contributing pixels
C 
C Other variables
C
C       q(3)		R  three vector of the start of the cut axis 	
C       r(3)		R  three vector to define direction of the cut axis 	
C       cut_bin		R  the size of the bins
C       width		R  the width of the cut
C       num_p		I  size of the input arrays 
C	num_e		I  size of the input arrays       
C       param		I  determines the plane projection
C
C-------------------------------------------------------------------------------
C
C OUTPUT:
C
C Arrays which have size (size_x, size_y) 	
C
C 	selectx		R  x index of bin for contributing pixel
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

      REAL FUNCTION CUT2D(axes, vx, vy, selectx, num_e, num_p, min_x, size_x, 
     &                    q, r, cut_bin, width, param, new_p, info, flag)

      INTEGER errcode, num_e, num_p, i, j, axes(3)
      INTEGER size_x, max_x, min_x, param, flag
      REAL selectx(num_e, num_p)
      REAL vx(num_e, num_p,5), vy(num_e, num_p,5)
      REAL p(3), q(2), r(2), width, cut_bin
      REAL new_x(2), new_y(2)
      REAL n(3), norm_x, norm_y, val, const, new_p(3), new_r(2)
      REAL testx_q, testy_q
      CHARACTER*50 info

C  Set routine flags to OK.

      cut2d  = 1.0
      flag   = 1

C  Check to see if a general cut is valid and set INFO and FLAG


      if ((q(1) .ne. r(1)) .and. (q(2) .ne. r(2))) then 
       if ((axes(1) .le. 3) .and. (axes(2) .le. 3)) then
        info = 'INFO: General Cut'
       else
        info = 'ERROR: Incompatable Axes for General Cut'
        flag = -1
       endif
      endif

      max_x = -1e4
      min_x =  1e4

      new_x(1) = r(1)-q(1)
      new_x(2) = r(2)-q(2)
      new_y(1) = r(2)-q(2)
      new_y(2) = q(1)-r(1)
 
C  Calculate the normalising factors for the x and y vectors
  
      norm_x = sqrt(new_x(1)*new_x(1)+new_x(2)*new_x(2))
      norm_y = sqrt(new_y(1)*new_y(1)+new_y(2)*new_y(2))


C  Param now detemines which axis to measure cut along

      new_p(1) = (q(1)*new_x(1) + q(2)*new_x(2))/norm_x
      new_p(2) = (q(1)*new_y(1) + q(2)*new_y(2))/norm_y
      new_p(3) = 0.0

      if (param .eq. 0.0) then 
      info = 'INFO: No projection of bins'
      do i=1, num_e
       do j=1, num_p
	new_r(1) = (new_x(1)*vx(i,j,1)+new_x(2)*vy(i,j,1))/norm_x
        new_r(2) = (new_y(1)*vx(i,j,1)+new_y(2)*vy(i,j,1))/norm_y
        if (abs(new_p(2)-new_r(2)) .le. (width/2.0)) then
         selectx(i,j) = nint((new_r(1) - new_p(1))/cut_bin)
         if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
         if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
        else 
          selectx(i,j) = -1e30
        endif
       enddo
      enddo
      endif

      if (param .eq. 1.0) then 
      info = 'INFO: Projection of bins onto x-axis'
      do i=1, num_e
       do j=1, num_p
	new_r(1) = (new_x(1)*vx(i,j,1)+new_x(2)*vy(i,j,1))/norm_x
        new_r(2) = (new_y(1)*vx(i,j,1)+new_y(2)*vy(i,j,1))/norm_y
        if (abs(new_p(2)-new_r(2)) .le. (width/2.0)) then
         selectx(i,j) = nint((vx(i,j,1) - q(1))/cut_bin)
         if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
         if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
        else 
          selectx(i,j) = -1e30
        endif
       enddo
      enddo
      new_p(1) = q(1)
      new_p(2) = q(2)
      endif

      if (param .eq. 2.0) then 
      info = 'INFO: Projection of bins onto y-axis'
      do i=1, num_e
       do j=1, num_p
	new_r(1) = (new_x(1)*vx(i,j,1)+new_x(2)*vy(i,j,1))/norm_x
        new_r(2) = (new_y(1)*vx(i,j,1)+new_y(2)*vy(i,j,1))/norm_y
        if (abs(new_p(2)-new_r(2)) .le. (width/2.0)) then
         selectx(i,j) = nint((vy(i,j,1) - q(2))/cut_bin)
         if (selectx(i,j) .lt. min_x) min_x = selectx(i,j)
         if (selectx(i,j) .gt. max_x) max_x = selectx(i,j)
        else 
          selectx(i,j) = -1e30
        endif
       enddo
      enddo

      new_p(1) = q(2)
      new_p(2) = q(1)

      endif

      size_x = max_x - min_x + 1

	write(*,*) ' size_x (in cut2d) = ',size_x

      return
      end
