        SUBROUTINE FUNCTION_G2(pars_get, pars_put)
        CHARACTER*8 fname,file*12
        INCLUDE 'FUNCTCOM.CMN'
	EXTERNAL pars_get, pars_put

        call FUNCTION_IN(pars_get)
        
!	write(6,1001)
!        accept 1002,kr,fname
!        file=fname//'.DAT'

	call module_get_string(pars_get, 'infile', file)

        open(unit=11,name=file,type='new',form='formatted')
40      rewind 11
        write(11,1009)in_title
        write(11,1010)
        write(6,1003)
        accept *, emin,emax
        npts=0
        do 10 i=1,lpt
        if(xin(i).lt.emin.or.xin(i).gt.emax)goto 10
        npts=npts+1
        write(11,1000)xin(i),yin(i),ein(i)
        if(npts.gt.3999)goto 20
10      continue
        write(6,1008)npts
        goto 30
20      write(6,1004)
        goto 40
30      close(unit=11)
        call FUNCTION_OUT(pars_put)
        STOP 'b2a> file closed'
1000    FORMAT(1P3E13.6)
1001    FORMAT(' b2a> filename ? ',$)
1002    FORMAT(q,a)
1003    FORMAT(' input, min mev, (,) max mev>')
1004    FORMAT('NO. OF POINTS READ WOULD EXCEED 3000')
1008    FORMAT('NO. OF POINTS,',I4)
1009    FORMAT(' ',A)
1010    FORMAT('M  (=MEV)',5X,'COUNTS',7X,'ERROR')
        END