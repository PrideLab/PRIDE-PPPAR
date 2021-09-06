!
!! AmBiguity Fractional cycle biases
type arscfg
  logical*1 lsearch
  integer*4 jd0,jd1,maxdel,minsav
  real*8 sod0,sod1,dintv
  real*8 minsec_common,wl_maxdev,wl_maxsig,wl_alpha,nl_maxdev,nl_maxsig,nl_alpha,cutoff,chisq,ratio
  character*3 fcbprn(MAXSAT)
  integer*4 fcbnprn
!! common files
  character*20 flnamb,flnneq,flncon,flnwlf
!! satellites
  integer*4 nprn
  character*3 prn(MAXSAT)
end type
