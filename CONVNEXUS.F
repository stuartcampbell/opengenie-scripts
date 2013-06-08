***********************************************************************************
*
*	CONVNEXUS.F
*
*	C. Moreton-Smith January 2000	
*
*	Example program to create a NeXus Level 0 file from any ISIS raw file.
*
*	This program reads some data from a standard ISIS raw data file and	
*	converts it by rewriting back into a NeXus level 0 file.
*
***********************************************************************************

	PROGRAM convert

	IMPLICIT NONE

	CALL GF_activate_session('NEXUS', 0, 1)
	CALL GF_select_source('hrp00273.raw')
	CALL GF_select_destination('hrp00273.nx0', 'NEXUS')
	CALL GF_get('h_FILE', 'FILE', 0)
	CALL GF_put('h_FILE', 'file', 0, 'new',
     +              'Level 0 NeXus conversion of an ISIS raw file')
	CALL GF_deactivate_session()
	
	END
