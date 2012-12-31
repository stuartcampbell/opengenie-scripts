C-------------------------------------------------------------------------------
C
C ROUTINE: BIN1D
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
C 	intensity	R	intensity in pixels
C 	error		R	error on intensity in pixels
C 	selectx		R	x index of bin for contributing pixel
C
C	1)	here x refers to the cut axis
C	2)	selectx, set to NULL_DATA for non-contributing pixels
C 
C Other variables
C
C	min_x		I  minimum value of selectx
C 	size_x		I  extend of selectx
C       bin_size	R  the size of the bins
C       p(3)		R  three vector of the start of the cut axis 	
C       num_p		I  size of the input arrays 
C	num_e		I  size of the input arrays       
C
C-------------------------------------------------------------------------------
C
C OUTPUT:
C
C Arrays which have size (size_x) 	
C
C	int		R  intensity in bins
C	err		R  errors on intensity in bins
C       count		R  counts how many pixels per bin
C       x_values	R  the x_vales of the bins
C
C Other variables
C
C       info		C  string containing error message
C	flag		I  flag to show sucessful execution of routine
C
C-------------------------------------------------------------------------------

      REAL FUNCTION BIN1D(intensity, error, num_e, num_p,  
     &                    min_x, size_x, p, bin_size,
     &                    int, err, count, selectx, x_values, info, flag)

      INTEGER num_e, num_p, i, j, min_x, size_x, m, flag
      REAL selectx(num_e, num_p)
      REAL intensity(num_p, num_e), error(num_p, num_e), bin_size
      REAL int(size_x), err(size_x), count(size_x)
      REAL x_values(size_x), p(3)
      CHARACTER*50 info

C  Set the routine flag to OK.

      bin1d = 1.0
      flag  = 1

	WRITE(*,*) ' INSIDE BIN1D...'

C  Perform the summation loop

      do j=1,num_p
       do i=1,num_e            
        if ((selectx(i,j) .ne. -1e30) .and. (intensity(j,i) .ne. -1e30)) then
          m = selectx(i,j) - min_x + 1
          int(m) = int(m) + intensity(j,i) 
          err(m) = err(m) + error(j,i)*error(j,i)
          count(m) = count(m) + 1.0
        endif
       end do
      end do

	WRITE(*,*) ' INSIDE BIN1D...1ST LOOP'

C  Perform the normalisation loop      

	WRITE(*,*) ' SIZE_X = ',SIZE_X

      do i=1, size_x
       if (count(i) .ne. 0.0) then
        int(i) = int(i)/count(i)
        err(i) = SQRT(err(i))/count(i)
       else
        int(i) = -1e30
        err(i) = 0.0
       endif 
       x_values(i)   = (i+min_x-1)*bin_size + p(1)
      enddo

      return
      end

