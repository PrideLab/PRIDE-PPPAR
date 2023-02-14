type satellite
  integer*4 iptatx
  character*3 prn
  character*20 typ
  real*8 var
!! antenna
  real*8 xyz(3,2), pcc(2)
!! satellite attitude
  real*8 xscf(3), yscf(3), zscf(3)
  real*8 nadir
end type
