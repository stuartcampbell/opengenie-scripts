$!
$!	Can't use SET NOVERIFY HERE -- GENIE does not understand this command
$!
$!***************************************************************************
$!
$!	BIN.COM		March 3, 1994
$!	J. A. Goldstone
$!
$!	Bins <run#>.RAW files for the ENGIN instrument and
$!	writes the binned file to <run#>.HIS file for later
$!	use.
$!
$!***************************************************************************
$!
$!	This command file is driven by ISISBIN.COM, which creates a 
$!	temporary file named BIN.BATCH which invokes GENIE and then calls
$!	this command file to do the actual binning (focussing) of the data
$!	from each detector into a final spectrum.
$!
$!	This command file can also be invoked interactively from within 
$!	GENIE.  It will prompt for the desired input.  Binned data files
$!	will be written to the current default directory.
$!
$!	Binning sections stolen from ENGINFOCL.COM and ENGINFOCRT.COM.  
$!	It must be updated with the appropriate no sample and vanadium 
$!	spectra as necessary.  The left bank is binned first and will 
$!	be in block 1 of the .HIS file and thr right bank will be in
$!	block 2 of the .HIS file.
$!
$!***************************************************************************
$!
$ mpath=2
$ icount=1
$ ibin = p1
$ enter:
$ INQUIRE mrun  ">> Enter next run number "
$ INQUIRE mpath ">> Enter 1 if last run number, 2 if not"
$ mrun'icount'=mrun
$ icount=icount+1
$ GOTO (start,enter)mpath
!
$ start:
$ icount=icount-1
!
$DO K=1,icount
$ nwork=1
$ NRUN=mrun'k'
!
> ASS 'nrun'
>@utilityn:focusl_10_98 nwork nrun
$ GOTO (change1,nochange1) ibin
$ CHANGE1:
>   REBIN w'nwork' 1000 [0.001] 19000
$ NOCHANGE1:
>write/min w'nwork' 'nrun'l.his
>write/close
>write/min w'nwork' 'nrun'.his
!$nwork=2
> ASS 'nrun'
>@utilityn:focusr_10_98.com nwork nrun
$ GOTO (change2,nochange2) ibin
$ CHANGE2:
>   REBIN w'nwork' 1000 [0.001] 19000
$ NOCHANGE2:
>write w'nwork'
>write/close
>write/min w'nwork' 'nrun'r.his
>write/close
$ENDDO
