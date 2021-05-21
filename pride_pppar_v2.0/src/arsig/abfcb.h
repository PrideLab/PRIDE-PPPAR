!
!! AmBiguity Fractional cycle biases
type abfcb
  logical*1 lsearch
  integer*4 jd0,jd1,maxdel,minsav
  real*8 sod0,sod1,dintv
  real*8 minsec_common,wl_maxdev,wl_maxsig,wl_alpha,nl_maxdev,nl_maxsig,nl_alpha,cutoff,chisq,ratio
!! common files
  character*20 flnamb,flnneq,flnfcb,flncon,flnwlf
!! stations & satellites
  integer*4 nsit,nprn,npp,iss(2,MAXPP),wcad(MAXPP),ncad(MAXPP)
  character*3 prn(MAXSAT)
  real*8 fwl(MAXPP),vwl(MAXPP),swl(MAXPP)
  integer*4 isc(MAXPP)
  real*8 fnl(MAXPP),vnl(MAXPP),snl(MAXPP)
end type
