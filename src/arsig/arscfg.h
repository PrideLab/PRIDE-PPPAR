!
!! AmBiguity Fractional cycle biases
type arscfg
  logical*1 lsearch, lverbose, lambsvm
  integer*4 jd0,jd1,maxdel,minsav
  real*8 sod0,sod1,dintv
  real*8 minsec_common,wl_maxdev,wl_maxsig,wl_alpha,nl_maxdev,nl_maxsig,nl_alpha,cutoff,chisq,ratio
  character*3 fcbprn(MAXSAT)
  integer*4 fcbnprn
!! common files
  character*20 flnamb,flnneq,flncon,flnwlf
  real*8 features(8)
!! # of obervations
  real*8 nobs
!! satellites
  integer*4 nprn
  character*3 prn(MAXSAT)
!! statistics
  integer*4 ntot_G,nwlfx_G,nwnfx_G
  integer*4 ntot_E,nwlfx_E,nwnfx_E
  integer*4 ntot_C,nwlfx_C,nwnfx_C
  integer*4 ntot_3,nwlfx_3,nwnfx_3
  integer*4 ntot_J,nwlfx_J,nwnfx_J
  integer*4 ntot_ind,nwl_ind,nwn_ind
end type

