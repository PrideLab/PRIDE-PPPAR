type lsqcfg
!
!! start & stop time
  integer*4 jd0,jd0_ses,jd1
  real*8 dintv,sod0,sod0_ses,sod1,seslen
!
!! number of satellite & PRNs
  integer*4 nprn
  character*3 prn(MAXSAT)
!
!! number of stations
  integer*4 nsit
!
!! ztd model
  character*15 ztdmod,htgmod
!
!! common file
  character*20 flnorb,flnsck,flnrck,flnpos,flnneq,flnerp
  character*20 flnztd,flnamb,flnres,flncon,flnhtg,flnvmf,flnfcb
  character*20 flnrck_G,flnrck_R,flnrck_E,flnrck_C,flnrck_3,flnrck_J
!
!! keep bias
  logical*1 lrmbias
!
!! ion model
  logical*1 lioh
  character*20 flntec
!
!!
  character*20 sys
  integer*4 sysnum
end type
