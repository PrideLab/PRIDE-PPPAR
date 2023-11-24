!
!! antenna atx
!
type antatx
  character*20 antnam
  character*20 antnum 
  integer*4    nfreq
  character*3  frq(MAXFRQ, MAXSYS)
  real*8       zen1, zen2, dzen, dazi
  real*8       neu(3, MAXFRQ, MAXSYS)
  real*8       pcv(50, 0:200, MAXFRQ, MAXSYS)
  character*5  sys_multi
end type
