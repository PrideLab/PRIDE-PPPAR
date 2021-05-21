!
!! resdig configure
!
type rescfg
  character*10 obstyp
!
!! start & stop time
  integer*4 jd0
  real*8 sod0,dintv,seslen
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
  real*8 jump,xres
!
!! common file
  logical*1 lupd
  integer*4 lfnres,lfnstt,lfnrhd
end type
