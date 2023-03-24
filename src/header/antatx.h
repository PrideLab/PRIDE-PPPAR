!
!! antenna atx
!
type antatx
  character*20 antnam
  character*20  antnum 
  integer*4    nfreq
  character*3  frq(10,5)
  real*8 zen1,zen2,dzen,dazi
  real*8 neu(3,10,5)
  real*8 pcv(50,0:200,10,5)
  character*5 sys_multi
  character*5 sys_multi2
end type
