!
!! Binary orbit header
!
type orbhdr
!
!! satelilte 
  character*10 :: sattyp = ''   ! satellite type  
  integer*4 nprn                ! # of satellites
  character*3 prn(MAXSAT)       ! satellite PRN
!
!! system tag
  character*80 :: iers = ''     ! IERS Conventions
!
!! time tag
  integer*4 jd0,jd1             !
  real*8 sod0,sod1              ! start & end time
!
!! sp3 interval (second)
  real*8 dintv
end type
!
!! orbit block
!
type sp3block
    integer*4 jd
    real*8    sod
    real*8    x(6,MAXSAT)  ! x, y, z, vx, vy, vz
    logical*1 flag(MAXSAT) ! lost sat : false ; have sat :true; 
end type

