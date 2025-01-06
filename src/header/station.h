type station
  character*4 name
  character*2 skd
  integer*4 iunit, lfnjmp, ikin, iion, imet, iptatx
  integer*4 pospd
  real*8 x(6), dx0(3), rx0, geod(3)
!
!! LEO science reference frame
  real*8 xsrf(3), ysrf(3), zsrf(3)
!
!! related file names
  character*256 obsfil, metfil
  character*20  kinfil, ionfil
!
!! receiver & antenna
  real*8 enu0(3), enu(3,2), rot_l2f(3,3)
  real*8 enu_sys(3,2,MAXSYS)
  character*20 rectyp, anttyp
!
!! mhm model
  real*8 mhm(360, 90, 2, 5)
!
!! clock correction
  real*8 dclk0, rclock_G, rclock_R, rclock_E, rclock_C, rclock_3, rclock_J, qrck
!
!! observation info
  real*8 sigr, sigp, cutoff
!
!! meteorology info
!  integer*4   nmet
!  character*2 mete(10)
!
!! meteorology
  character*3 map
  integer*4   nvm, jdv(maxday*4+1)
  real*8 ztdcor, dztd0, qztd, dhtg0, qhtg, p0, t0, hr0, undu
  real*8 sodv(maxday*4+1), vm1(4,maxday*4+1)
!
!! phase wind-up
  logical*1 first(maxsat)
  real*8   prephi(maxsat)
!
!! ocean load
  real*8 olc(11,6), rlat, rlon
  character*20 otlfil
end type 
