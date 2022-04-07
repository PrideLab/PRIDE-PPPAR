!
!! resdig configure
!
type rescfg
  character*20 obstyp,flnrhd
!
!! start & stop time
  integer*4 jd0
  real*8 sod0,dintv
!
!! number of satellite & PRNs
  integer*4 nprn
  character*3 prn(MAXSAT)
!
!! number of stations
  character*4 snam
!
!! screening control
  integer*4 nsht
  real*8 jump
!
!! common file
  integer*4 lfnres,lfnstt,lfnrhd
end type
