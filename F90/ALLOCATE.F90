
 SUBROUTINE test_alloc(pars_get, pars_put)
 
 INCLUDE 'genie_modules_f90.inc'
 
 REAL, ALLOCATABLE 	:: a(:)
 INTEGER 			:: n

  CALL MODULE_GET_INT(pars_get, 'n', n)

  ALLOCATE(a(n))
  
  DO i = 1, n
   a(i) = i
  ENDDO 

  CALL MODULE_PUT_REAL_ARRAY(pars_put,'a', a, n)

 RETURN
 END
