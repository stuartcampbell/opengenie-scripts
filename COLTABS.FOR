	SUBROUTINE COLOURS(pars_get,pars_put)

	INCLUDE 'genie_modules.inc'
	INCLUDE 'coltabs.inc'

      INTEGER n, dims_array(2)

	EXTERNAL pars_get, pars_put

	n = 2                
       
	dims_array(1) = 3
	dims_array(2) = 256

	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'bgyrw',bgyrw,dims_array,n)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'heat',heat,dims_array,n)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'spectrum',spectrum,dims_array,n)
	CALL MODULE_PUT_ND_REAL_ARRAY(pars_put,'serp',serp,dims_array,n)
                                                  
	RETURN
	END
