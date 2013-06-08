!!	Genie Routine to read in IN5 data from ILL   as CROSSX output
!!	Instrument abbr-name is IN5 and ext is ILL
$ inquire nrun "Run number"
$ nwr=1
$ ngrps=8
$ do n=1,ngrps
$ v1=1
> lo w'n' inx'nrun'.dat i_i:lo_in10_'n'		!!load data
$ goto (wr_1,wr_2) nwr
$ wr_1:
> @i_p:wr_int INX 'nrun' 'n' IQW  	!!Write Genie file
$ nwr=2
$ goto wr_3
$ wr_2:
> wr w'n'
$ wr_3:
$ end do
> wr/close
