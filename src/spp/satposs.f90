
! Satellite positions and clocks --------------------------------------------
! compute satellite positions, velocities and clocks
! args   : gtime_t teph     I   time to select ephemeris (gpst)
!          obsd_t *obs      I   observation data
!          integer*4 n      I   number of observation data
!          nav_t  *nav      I   navigation data
!          integer*4 ephopt I   ephemeris option (EPHOPT_???)
!          real*8 *rs       O   satellite positions and velocities (ecef)
!          real*8 *dts      O   satellite clocks
!          real*8 *var      O   sat position and clock error variances (m^2)
!          integer*4 *svh   O   sat health flag (-1:correction not available)
! return : none
! notes  : rs ((0:2)+i*6)= obs(i) sat position {x,y,z} (m)
!          rs ((3:5)+i*6)= obs(i) sat velocity {vx,vy,vz} (m/s)
!          dts((0:1)+i*2)= obs(i) sat clock {bias,drift} (s|s/s)
!          var(i)        = obs(i) sat position and clock error variance (m^2)
!          svh(i)        = obs(i) sat health flag
!          if no navigation data, set 0 to rs(_), dts(_), var(_) and svh(_)
!          satellite position and clock are values at signal transmission time
!          satellite position is referenced to antenna phase center
!          satellite clock does not include code bias correction (tgd or bgd)
!          any pseudorange and broadcast ephemeris are always needed to get
!          signal transmission time
!----------------------------------------------------------------------------
subroutine satposs(teph, obs, n, nav, ephopt, rs, dts, var, svh)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: teph
type(obsd_t), intent(in) :: obs(n)
integer*4, intent(in) :: n, ephopt
type(nav_t), intent(in) :: nav
real*8, intent(out) :: rs(n,6), dts(n,2), var(n)
integer*4, intent(out) :: svh(n)
type(gtime_t) time(2*MAXOBS),timeadd
real*8 dt,pr,SQR
integer*4 i,j,ulimit,info
external :: timeadd,SQR
time=gtime_t(0,0.d0)
ulimit=min(n,2*MAXOBS)
do i=1,ulimit
    rs(i,:)=0.d0; dts(i,:)=0.d0
    var(i)=0.d0; svh(i)=0
    
    ! search any psuedorange 
    pr=0.d0
    do j=1,NFREQ
        pr=obs(i)%P(j)
        if(dabs(pr-0.d0)>=1d-20) exit
    enddo
    if (j>NFREQ) cycle
    
    ! transmission time by satellite clock 
    time(i)=timeadd(obs(i)%time,-pr/CLIGHT)
        
    ! satellite clock bias by broadcast ephemeris 
    call ephclk(time(i),teph,obs(i)%sat,nav,dt,info)
    if (info==0) cycle
    
    time(i)=timeadd(time(i),-dt)
    
    ! satellite position and clock at transmission time 
    call satpos(time(i),teph,obs(i)%sat,ephopt,nav,rs(i,:),dts(i,:),var(i),svh(i),info)
    if(info==0) cycle
    
    ! if no precise clock available, use broadcast clock instead 
    if (dabs(dts(i,1)-0.d0)<=1d-20)then
        call ephclk(time(i),teph,obs(i)%sat,nav,dts(i,1),info)
        if (info==0) cycle
        dts(i,2)=0.d0
        var(1)=SQR(STD_BRDCCLK)
    endif        
enddo
end subroutine

! satellite position and clock ----------------------------------------------
! compute satellite position, velocity and clock
! args   : gtime_t time     I   time (gpst)
!          gtime_t teph     I   time to select ephemeris (gpst)
!          integer*4 sat    I   satellite number
!          nav_t  *nav      I   navigation data
!          integer*4 ephopt I   ephemeris option (EPHOPT_???)
!          real*8 *rs       O   sat position and velocity (ecef)
!                               {x,y,z,vx,vy,vz} (m|m/s)
!          real*8 *dts      O   sat clock {bias,drift} (s|s/s)
!          real*8 *var      O   sat position and clock error variance (m^2)
!          integer*4 *svh   O   sat health flag (-1:correction not available)
! return : status (1:ok,0:error)
! notes  : satellite position is referenced to antenna phase center
!          satellite clock does not include code bias correction (tgd or bgd)
!----------------------------------------------------------------------------
subroutine satpos(time, teph, sat, ephopt, nav, rs, dts, var, svh, stat)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time, teph
integer*4, intent(in) :: sat, ephopt
type(nav_t), intent(in) :: nav
real*8, intent(out) :: rs(6), dts(2), var
integer*4, intent(out) :: svh, stat

svh=0
select case(ephopt)
case(EPHOPT_BRDC)
    call ephpos(time,teph,sat,nav,-1,rs,dts,var,svh,stat); return
end select
svh=-1
stat=0
end subroutine
