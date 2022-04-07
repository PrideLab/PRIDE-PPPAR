type station
  character*4 name
  character*2 skd
  integer*4 iunit,lfnjmp,ikin,iion,imet,iptatx
  real*8 x(6),dx0(3),geod(3)
!
!! related file names
  character*256 obsfil,metfil
  character*20 kinfil,ionfil
!
!! receiver & antenna
  real*8 enu0(3),enu(3,2),rot_l2f(3,3)
  real*8 enu_G(3,2),enu_R(3,2),enu_E(3,2),enu_C(3,2),enu_J(3,2)
  character*20 rectyp,anttyp
!
!! clock correction
  real*8 dclk0,rclock_G,rclock_R,rclock_E,rclock_C,rclock_3,rclock_J
!
!! observation info
  real*8 sigr,sigp,cutoff
!
!! meteorology info
!  integer*4 nmet
!  character*2 mete(10)
!
!! meteorology
  character*3 map
  integer*4 nvm,jdv(13)
  real*8 ztdcor,dztd0,qztd,dhtg0,qhtg,p0,t0,hr0,undu,sodv(13),vm1(4,13)
!
!! phase wind-up
  logical*1 first(maxsat)
  real*8 prephi(maxsat)
!
!! ocean load
  real*8 olc(11,6),rlat,rlon
end type 
