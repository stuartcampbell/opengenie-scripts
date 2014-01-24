$ SET NOVERIFY
$ ON CONTROL_Y THEN GOTO CTRL_EXIT
$!
$!***************************************************************************
$!	ISISBIN.COM	March 3, 1994
$!	Joyce A. Goldstone
$!
$! 	Automate the binning of TEST BEAM POSITION 2 data -- prototype ENGIN
$!	This version for two detector banks.  
$!	
$!***************************************************************************
$!
$!	Command: @ISISBIN P1 P2 P3
$!
$!	where	IF P1 = FILE, then input is read from a file named by
$!		the parameter P2.  File specification must sufficient to
$!		find the file from the current default directory.  For
$!		example, one can activate an input file on the same node and
$!		disk by giving the directory and name, e.g. [jag.junk]run.input.
$!		It is recommended that this is the minimum infornation
$!		supplied.
$!
$!		IF P1 = anything else, then it is assumed that input is to
$!		done interactively and the user is asked for the required
$!		information: directory and run numbers.
$!
$!		IF P3 present and = 1 indicates that the bin grouping should be
$!		changed to [0.001]
$!
$!	Examples:	@ISISBIN
$!			    Totally interactive -- information requested
$!			@ISISBIN FILE RUN.INPUT
$!			    Totally automatic -- information from a file
$!			    named RUN.INPUT
$!			@ISISBIN FILE RUN.INPUT 1
$!			    Totally automatic -- information from a file
$!			    named RUN.INPUT, bin grouping will be changed
$!		
$!***************************************************************************
$!
$!	Input file structure:  This command file needs two pieces of 
$!	information -- the directory for the output files and the list of
$!	run numbers to be binned.  The input file should contain this 
$!	on consecutive lines.  A blank line indicates the end of the run
$!	number list.
$!	
$!	RUN.INPUT example
$!		[JAG.JUNK]
$!		[jag.ccsl]teb_24mar94.iparm
$!		4580 1 2 
$!		4581 1 2 
$!		4582 1 2
$!
$!		<EOF>
$!
$!	This can also be used as the input file for ISISREF by modifying the
$!	the run number lines and adding lines below the run number list.
$!	See ISISREF documentation for further information.
$!
$!***************************************************************************
$!
$!	The data are placed into two histogram in blocks 1 and 2 of the common
$!	.HIS file.  The convention is that the left detector is in block 1 and
$!	right detector is in block 2.  Left and right are determined by looking
$!	at the target from the sample position.  Thus the left detector is on
$!	the HRPD side.
$!
$!***************************************************************************
$!
$ SPACE = " "				! Separator for input from a file or
$!					   interactive for the run numbers
$ SEP = "/"				! List separator for the runlist
$ ibin = 2				! Flag for rebin 1=yes, 2=no
$!
$ RUNLIST = ""
$ IF p1 .eqs. "FILE"			! Check to see if file input
$   THEN
$	IF (P3 .eq. 1) THEN ibin=1	! Set flag to rebin
$	IFILE = F$PARSE(P2)
$	OPEN/READ/ERROR=NO_FILE INP 'IFILE'
$	READ INP adir     
$	READ INP iparm			! Get the instrument parameter file
$	FIRST = "TRUE"
$ READRUN:				! Now read the individual run numbers
$	READ INP LINE
$	IF (LINE .eqs. "")		! Blank line ends input
$	   THEN 
$	   CLOSE INP
$	   GOTO START
$	ENDIF
$   	RUNNUM =F$ELEMENT(0," ",LINE)	! Construct the run list for binning
$	IF FIRST .EQS. "TRUE"
$	   THEN
$	   RUNLIST = RUNLIST + RUNNUM
$	   FIRST = "FALSE"
$	ELSE
$	   RUNLIST = RUNLIST +SEP +RUNNUM
$	ENDIF
$	GOTO READRUN
$!
$ ELSE	 			! If NOT FILE, request the directory,runlist
$     write sys$output "Enter directory to store binned data file(s). "
$     inquire adir "(Ex: [jag.mit])"
$     write sys$output ""
$     write sys$output "Enter the run number(s) (slash separated) of the data to be binned"
$     inquire runlist "(Ex: 1218/1219/1220) "
$     write sys$output ""
$     write sys$output "Do you want to change the bin grouping to [0.001]? "
$     inquire wbin "(1 =YES, anything else, including <CR> means no) "
$     write sys$output ""
$ ENDIF
$!
$     IF (wbin .eq. 1 ) THEN ibin = 1		! Set rebin flag to 1
$!
$!     Write a temporary batch file saving user entered data
$!
$START:
$ open/write temp 'adir'bin.batch
$ write temp "$ SET NOVERIFY"
$ write temp "$ SET DEFAULT ''adir'"
$ write temp "$ DELETE bin.log;*"
$ write temp "$ GENIE"
$ write temp "@utilityn:focus25.com ''ibin'"
$!
$ count=0			! Parse the run number list for each run number
$ runcount:
$    runnum =f$element ( count,SEP,runlist)
$    IF runnum .eqs. SEP THEN GOTO runwrite
$    count = count +1 
$    GOTO runcount
$!
$ runwrite:
$ element = 0
$ get_run:
$    runnum =f$element ( element,SEP,runlist)
$    IF runnum .eqs. SEP THEN GOTO done
$    write temp "''runnum'"
$    element = element + 1
$    IF element .eq. count 
$	THEN write temp "1"			! 1 flags the end of run numbers
$       ELSE write temp "2"
$    ENDIF
$    GOTO get_run
$!
$ done:
$ write temp "exit"
$ write temp "$EOJ"
$ close temp
$!
$!	Submit the batch job to do the focussing and write out the file
$!
$ write sys$output "The binned file(s) will be saved as ''adir'run#.his"
$ write sys$output ""
$ submit/queue=PEARL$BATCH 'adir'bin.batch/delete
$EXIT                                                  
$! 
$CTRL_EXIT:
$   close temp
$   close inp
$EXIT
$!
$NO_FILE:
$    write sys$output "Input file ''IFILE' not found"
