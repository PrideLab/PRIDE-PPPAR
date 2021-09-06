type lsqcfg
!
!! start & stop time
  integer*4 jd0,jd_beg,jd1,jd_end
  real*8 dintv,sod0,sod_beg,sod1,sod_end,seslen
!
!! number of satellite & PRNs
  integer*4 nprn
  character*3 prn(MAXSAT)
!
!! ztd model
  character*15 ztdmod,htgmod
!
!! common file
  character*20 flnorb,flnsck,flnrck,flnpos,flnneq,flnerp
  character*20 flnztd,flnamb,flnres,flncon,flnhtg,flnvmf,flnfcb
  character*20 flnatt
  character*60 flnorb_real,flnsck_real,flnfcb_real,flnatt_real,flnerp_real,flnatx_real
!
!! read ATT
  logical*1 attuse
!
!! keep bias
  logical*1 lrmbias
!
!! ion model
  logical*1 lioh
  character*20 flntec
!
!! tide displacement
  logical*1 otluse
  character*60 tide
!
!! fixing ambiguity
  logical*1 fcbuse,pcowl
  character*3 fcbprn(MAXSAT)
  integer*4 fcbnprn,nconG,nconE,nconC2,nconC3,nconJ
!
!!
  character*20 sys
  integer*4 sysnum
!
!!
  character*15 edit
!
!!
  integer*4 maxdel,minsav
  real*8 minsec_common,wl_maxdev,wl_maxsig,wl_alpha,nl_maxdev,nl_maxsig,nl_alpha,cutoff,chisq,ratio
end type
