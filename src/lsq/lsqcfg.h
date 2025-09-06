type lsqcfg
!
!! start & stop time
  integer*4 jd0, jd_beg
  integer*4 jd1, jd_end
  real*8 dintv, seslen,dwnd
  real*8 sod0, sod_beg
  real*8 sod1, sod_end
!
!! number of satellite & PRNs
  integer*4    nprn
  character*3  prn(MAXSAT)
  character*10 sys, edit
!
!! parameter model
  character*15 rckmod, isbsys
  character*15 ztdmod
  character*15 htgmod
!
!! common file
  character*20 flnorb, flnsck, flnerp, flnatt, flnfcb, flnlat
  character*20 flncon, flnneq, flnvmf
  character*20 flnamb, flnpos, flnrck
  character*20 flnres, flnhtg, flnztd
  character*60 flnorb_real, flnsck_real, flnerp_real
  character*60 flnatt_real, flnfcb_real, flnatx_real, flnlat_real
!
!! read ATT/OSB/LEO
  logical*1 attuse
  logical*1 fcbuse
  logical*1 latuse
!
!! remove ambiguity parameters
  logical*1 lrmbias
!
!! remove observations without bias corrections 
  logical*1 del_nobia
!
!! ion model
  logical*1 lioh
  character*20 flntec
!
!! mhm model
  logical*1 mhmuse
  character*20 flnmhm
!
!! tide displacement
  logical*1 otluse
  character*60 tide
!
!! fixing ambiguity
  logical*1 pcowl, isini
  character*3 fcbprn(MAXSAT)
  integer*4 fcbnprn
  integer*4 nconG, nconE, nconC2, nconC3, nconJ
!
!! options
  integer*4 maxdel, minsav
  real*8 minsec_common, cutoff, chisq, ratio
  real*8 wl_maxdev, wl_maxsig, wl_alpha
  real*8 nl_maxdev, nl_maxsig, nl_alpha
  character*15 amb_at_dbd
end type
