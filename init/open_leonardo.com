$ wso:==write sys$output
$ purge/nolog leonardo.inc
$ define/nolog/user genie_gcl_init user$disk:[sic01.opengenie.init]open_leonardo.gcl
!$ opengenie
$ newgenie 
!$ opengenie "user$disk:[sic01.opengenie.init]open_leonardo.gcl"
       