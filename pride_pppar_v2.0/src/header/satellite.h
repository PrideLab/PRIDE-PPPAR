type satellite
  integer*4 iptatx
  character*3 prn
  character*1 sys
  character*20 typ
!! antenna
  real*8 xyz(3,2)
!! satellite attitude
  real*8 xscf(3),yscf(3),zscf(3)
!! clock correction
  real*8 sclock
end type
