$ define/user genie_gcl_init het$disk:[hetmgr.axp.opengenie]canvas_init.gcl
$ define/user genie_smalltalk_image het$disk:[hetmgr.axp.opengenie]canvas.im
$ @axplib$disk:[opengenie_test]genie_setup -l
  