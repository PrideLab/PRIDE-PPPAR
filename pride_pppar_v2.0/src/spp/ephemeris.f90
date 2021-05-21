!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
! constants and macros ------------------------------------------------------

module ephe_f90_
use rtkcmn_f90_
implicit none
real*8, parameter :: RE_GLO = 6378136.d0        ! radius of earth (m)            ref (2) 
real*8, parameter :: MU_GPS = 3.9860050d14      ! gravitational constant         ref (1) 
real*8, parameter :: MU_GLO = 3.9860044d14      ! gravitational constant         ref (2) 
real*8, parameter :: MU_GAL = 3.986004418d14    ! earth gravitational constant   ref (7) 
real*8, parameter :: MU_CMP = 3.986004418d14    ! earth gravitational constant   ref (9) 
real*8, parameter :: J2_GLO = 1.0826257d-3      ! 2nd zonal harmonic of geopot   ref (2) 

real*8, parameter :: OMGE_GLO = 7.292115d-5      ! earth angular velocity (rad/s) ref (2) 
real*8, parameter :: OMGE_GAL = 7.2921151467d-5  ! earth angular velocity (rad/s) ref (7) 
real*8, parameter :: OMGE_CMP = 7.292115d-5      ! earth angular velocity (rad/s) ref (9) 
real*8, parameter :: SIN_5 = -0.0871557427476582d0 ! dsin(-5.d0 deg) 
real*8, parameter :: COS_5 =  0.9961946980917456d0 ! dcos(-5.d0 deg) 

real*8, parameter :: ERREPH_GLO  = 5.d0          ! error of glonass ephemeris (m) 
real*8, parameter :: TSTEP       = 60.d0         ! integration step glonass ephemeris (s) 
real*8, parameter :: RTOL_KEPLER = 1d-13         ! relative tolerance for Kepler equation 

real*8, parameter :: DEFURASSR  = 0.15d0         ! default accurary of ssr corr (m) 
real*8, parameter :: MAXECORSSR = 10.d0          ! max orbit correction of ssr (m) 
real*8, parameter :: MAXCCORSSR = (1d-6*CLIGHT)  ! max clock correction of ssr (m) 
real*8, parameter :: MAXAGESSR  = 90.d0          ! max age of ssr orbit and clock (s) 
real*8, parameter :: MAXAGESSR_HRCLK = 10.d0     ! max age of ssr high-rate clock (s) 
real*8, parameter :: STD_BRDCCLK = 30.d0         ! error of broadcast clock (m) 
real*8, parameter :: STD_GAL_NAPA = 500.d0       ! error of galileo ephemeris for NAPA (m) 

integer*4, parameter :: MAX_ITER_KEPLER = 30     ! max number of iteration of Kelpler 

! ephemeris selections ------------------------------------------------------
integer*4, parameter :: eph_sel(6)=(/0,0,1,0,0,0/) ! GPS,GLO,GAL,QZS,BDS,SBS 

contains
!! variance by ura ephemeris (ref (1) 20.3.3.3.1.1) -------------------------
real*8 function var_uraeph(sys, ura)
implicit none
integer*4, intent(in) :: sys, ura
if(sys==SYS_GAL)then    ! galileo sisa (ref (7) 5.1.11) 
    if (ura<= 49) var_uraeph=SQR(ura*0.01d0)
    if (ura<= 74) var_uraeph=SQR(0.5+(ura-50)*0.02d0)
    if (ura<= 99) var_uraeph=SQR(1.d0+(ura-75)*0.04d0)
    if (ura<=125) var_uraeph=SQR(2.d0+(ura-100)*0.16d0)
    if (ura> 125) var_uraeph=SQR(STD_GAL_NAPA)
else    ! gps ura (ref (1) 20.3.3.3.1.1) 
    var_uraeph=rcond(ura<0.or.14<ura,SQR(6144.d0),SQR(ura_value(ura+1)))
endif
end function

! broadcast ephemeris to satellite clock bias -------------------------------
! compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
! args   : gtime_t time     I   time by satellite clock (gpst)
!          eph_t *eph       I   broadcast ephemeris
! return : satellite clock bias (s) without relativeity correction
! notes  : see ref (1),(7),(8)
!          satellite clock does not include relativity correction and tdg
!----------------------------------------------------------------------------
real*8 function eph2clk(time, eph)
implicit none
type(gtime_t),intent(in) :: time
type(eph_t), intent(in) :: eph
real*8 t
integer*4 i
t=timediff(time,eph%toc)
do i=1,2
    t=t-(eph%f0+eph%f1*t+eph%f2*t*t)
enddo
eph2clk=eph%f0+eph%f1*t+eph%f2*t*t
end function

! broadcast ephemeris to satellite position and clock bias --------------------
! compute satellite position and clock bias with broadcast ephemeris (gps,
! galileo, qzss)
! args   : gtime_t time     I   time (gpst)
!          eph_t *eph       I   broadcast ephemeris
!          real*8 *rs       O   satellite position (ecef) {x,y,z} (m)
!          real*8 *dts      O   satellite clock bias (s)
!          real*8 *var      O   satellite position and clock variance (m^2)
! return : none
! notes  : see ref (1),(7),(8)
!          satellite clock includes relativity correction without code bias
!          (tgd or bgd)
!-----------------------------------------------------------------------------
subroutine eph2pos(time, eph, rs, dts, var)
implicit none
type(gtime_t), intent(in) :: time
type(eph_t), intent(in) :: eph
real*8, intent(out) :: rs(3),dts,var
real*8 tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge1
real*8 xg,yg,zg,sino1,coso1
integer*4 n,sys,prn
if(eph%A<=0.d0)then
    rs=0.d0; dts=0.d0; var=0.d0; return
endif
tk=timediff(time, eph%toe)

call satsys(eph%sat,prn,sys)
select case(sys)
case (SYS_GAL)
    mu=MU_GAL; omge1=OMGE_GAL
case (SYS_CMP)
    mu=MU_CMP; omge1=OMGE_CMP
case default
    mu=MU_GPS; omge1=OMGE
end select
M=eph%M0+(dsqrt(mu/(eph%A*eph%A*eph%A))+eph%deln)*tk
n=0;E=M;Ek=0.d0

do while(dabs(E-Ek)>RTOL_KEPLER .and. n<MAX_ITER_KEPLER)
    Ek=E; E=E-(E-eph%e*dsin(E)-M)/(1.d0-eph%e*dcos(E))
    n=n+1
enddo

if(n>=MAX_ITER_KEPLER) return
sinE=dsin(E); cosE=dcos(E)

u=datan2(dsqrt(1.d0-eph%e*eph%e)*sinE,cosE-eph%e)+eph%omg
r=eph%A*(1.d0-eph%e*cosE)
i=eph%i0+eph%idot*tk
sin2u=dsin(2.d0*u); cos2u=dcos(2.d0*u)
u=u+eph%cus*sin2u+eph%cuc*cos2u
r=r+eph%crs*sin2u+eph%crc*cos2u
i=i+eph%cis*sin2u+eph%cic*cos2u
x=r*dcos(u); y=r*dsin(u); cosi=dcos(i)

if (sys==SYS_CMP .and. (eph%flag==2 .or. (eph%flag==0 .and. prn<=5)))then  ! to be verified
	O=eph%OMG0+eph%OMGd*tk-omge1*eph%toes
	sinO=dsin(O); cosO=dcos(O)
	xg=x*cosO-y*cosi*sinO
	yg=x*sinO+y*cosi*cosO
	zg=y*dsin(i)
	sino1=dsin(omge1*tk); coso1=dcos(omge1*tk)
	rs(1)= xg*coso1+yg*sino1*COS_5+zg*sino1*SIN_5
	rs(2)=-xg*sino1+yg*coso1*COS_5+zg*coso1*SIN_5
	rs(3)=-yg*SIN_5+zg*COS_5
else
	O=eph%OMG0+(eph%OMGd-omge1)*tk-omge1*eph%toes
    sinO=dsin(O); cosO=dcos(O)
    rs(1)=x*cosO-y*cosi*sinO
    rs(2)=x*sinO+y*cosi*cosO
    rs(3)=y*dsin(i)
endif

tk=timediff(time,eph%toc)
dts=eph%f0+eph%f1*tk+eph%f2*tk*tk
    
! relativity correction 
dts=dts-2.d0*dsqrt(mu*eph%A)*eph%e*sinE/SQR(CLIGHT)
    
! position and clock error variance 
var=var_uraeph(sys,eph%sva)
end subroutine

!! glonass orbit differential equations --------------------------------------
subroutine deq(x, xdot, acc)
implicit none
real*8, intent(in) :: x(6),acc(3)
real*8, intent(out) :: xdot(6)
real*8 a,b,c,r2,r3,omg2
r2=dot(x,x,3)
r3=r2*dsqrt(r2)
omg2=SQR(OMGE_GLO)
if(r2<=0.d0)then
    xdot=0.d0; return
endif
a=1.5d0*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3  ! 3/2*J2*mu*Ae^2/r^5 
b=5.d0*x(3)*x(3)/r2                      ! 5*z^2/r^2 
c=-MU_GLO/r3-a*(1.d0-b)                  ! -mu/r^3-a(1-b) 
xdot(1)=x(4); xdot(2)=x(5); xdot(3)=x(6)
xdot(4)=(c+omg2)*x(1)+2.d0*OMGE_GLO*x(5)+acc(1)
xdot(5)=(c+omg2)*x(2)-2.d0*OMGE_GLO*x(4)+acc(2)
xdot(6)=(c-2.d0*a)*x(3)+acc(3)
end subroutine

!! glonass position and velocity by numerical integration --------------------
subroutine glorbit(t, x, acc)
implicit none
real*8, intent(in) :: t,acc(3)
real*8, intent(out) :: x(6)
real*8 k1(6),k2(6),k3(6),k4(6),w(6)
integer*4 i

call deq(x,k1,acc); w=x+k1*t/2.d0
call deq(w,k2,acc); w=x+k2*t/2.d0
call deq(w,k3,acc); w=x+k3*t
call deq(w,k4,acc)
x=x+(k1+2.d0*k2+2.d0*k3+k4)*t/6.d0
end subroutine

! glonass ephemeris to satellite clock bias -----------------------------------
! compute satellite clock bias with glonass ephemeris
! args   : gtime_t time     I   time by satellite clock (gpst)
!          geph_t *geph     I   glonass ephemeris
! return : satellite clock bias (s)
! notes  : see ref (2)
!-----------------------------------------------------------------------------
real*8 function geph2clk(time, geph)
implicit none
type(gtime_t), intent(in) :: time
type(geph_t), intent(in) :: geph
real*8 t
integer*4 i
t=timediff(time,geph%toe)
do i=1,2
    t=t-(-geph%taun+geph%gamn*t)
enddo
geph2clk=-geph%taun+geph%gamn*t
end function

! glonass ephemeris to satellite position and clock bias ----------------------
! compute satellite position and clock bias with glonass ephemeris
! args   : gtime_t time     I   time (gpst)
!          geph_t *geph     I   glonass ephemeris
!          real*8 *rs       O   satellite position {x,y,z} (ecef) (m)
!          real*8 *dts      O   satellite clock bias (s)
!          real*8 *var      O   satellite position and clock variance (m^2)
! return : none
! notes  : see ref (2)
!-----------------------------------------------------------------------------
subroutine geph2pos(time, geph, rs, dts, var)
implicit none
type(gtime_t), intent(in) :: time
type(geph_t), intent(in) :: geph
real*8, intent(out) :: rs(3),dts,var
real*8 t,tt,x(6)
integer*4 i

t=timediff(time,geph%toe)
dts=-geph%taun+geph%gamn*t
do i=1,3
    x(i  )=geph%pos(i)
    x(i+3)=geph%vel(i)
enddo
tt=rcond(t<0.d0,-TSTEP,TSTEP)
do while(dabs(t)>1d-9)
    if (dabs(t)<TSTEP) tt=t
    call glorbit(tt,x,geph%acc)
    t=t-tt
enddo
rs=x(1:3)
var=SQR(ERREPH_GLO)
end subroutine

!! select ephememeris --------------------------------------------------------
function seleph(time, sat, iode, nav)
implicit none
type(eph_t), pointer :: seleph
type(gtime_t), intent(in) :: time
integer*4, intent(in) :: sat,iode
type(nav_t), intent(in) :: nav
real*8 t,tmax,tmin
integer*4 :: i,j
integer*4 prn, sys, sel
j=-1; sel=0
call satsys(sat,prn,sys)
select case(sys)
case(SYS_GPS)
    tmax=MAXDTOE+1.d0    ; sel=eph_sel(1)
case(SYS_GAL)
    tmax=MAXDTOE_GAL     ; sel=eph_sel(3)
case(SYS_QZS)
    tmax=MAXDTOE_QZS+1.d0; sel=eph_sel(4)
case(SYS_CMP)
    tmax=MAXDTOE_CMP+1.d0; sel=eph_sel(5)
case default
    tmax=MAXDTOE+1.d0
end select
tmin=tmax+1.d0

do i=1,nav%n
    if (nav%eph(i)%sat/=sat) cycle
    if (iode>=0 .and. nav%eph(i)%iode/=iode) cycle
    if (sys==SYS_GAL .and. sel/=0)then
        if (sel==1 .and. and(nav%eph(i)%code,lshift(1,9))==0) cycle  ! I/NAV 
        if (sel==2 .and. and(nav%eph(i)%code,lshift(1,8))==0) cycle  ! F/NAV 
    endif
    t=dabs(timediff(nav%eph(i)%toe,time))
    if (t>tmax) cycle
    if (iode>=0)then
        seleph=>nav%eph(i); return
    endif
    if (t<=tmin)then  ! toe closest to time 
        j=i; tmin=t
    endif
enddo
if (iode>=0 .or. j<0)then
    nullify(seleph); return
endif
seleph=>nav%eph(j)
end function

! select glonass ephememeris ------------------------------------------------
function selgeph(time, sat, iode, nav)
implicit none
type(geph_t), pointer :: selgeph
type(gtime_t), intent(in) :: time
integer*4, intent(in) :: sat,iode
type(nav_t), intent(in) :: nav
real*8 :: t,tmax,tmin
integer*4 :: i,j
j=-1
tmax=MAXDTOE_GLO; tmin=tmax+1.d0
do i=1,nav%ng
    if (nav%geph(i)%sat/=sat) cycle
    if (iode>=0 .and. nav%geph(i)%iode/=iode) cycle
    t=dabs(timediff(nav%geph(i)%toe,time))
    if (t>tmax) cycle
    if (iode>=0)then
        selgeph=>nav%geph(i); return
    endif
    if (t<=tmin)then  ! toe closest to time 
        j=i; tmin=t
    endif
enddo
if (iode>=0 .or. j<0)then
    nullify(selgeph); return
endif
selgeph=>nav%geph(j)
end function

! satellite clock with broadcast ephemeris ----------------------------------
subroutine ephclk(time, teph, sat, nav, dts, stat)
implicit none
type(gtime_t), intent(in) :: time,teph
integer*4, intent(in) :: sat
type(nav_t), intent(in) :: nav
real*8, intent(out) :: dts
integer*4, intent(out) :: stat
type(eph_t), pointer :: eph
type(geph_t), pointer :: geph
integer*4 prn,sys
call satsys(sat,prn,sys)

if (sys==SYS_GPS .or. sys==SYS_GAL .or. sys==SYS_QZS .or. sys==SYS_CMP)then
    eph=>seleph(teph,sat,-1,nav)
    if (.not.associated(eph))then
        stat=0; return
    endif
    dts=eph2clk(time,eph)
elseif(sys==SYS_GLO)then
    geph=>selgeph(teph,sat,-1,nav)
    if (.not.associated(geph))then
        stat=0; return
    endif
    dts=geph2clk(time,geph)
else
    stat=0; return
endif
stat=1
end subroutine

! satellite position and clock by broadcast ephemeris -----------------------
subroutine ephpos(time_in, teph, sat, nav, iode, rs, dts, var, svh, stat)
implicit none
type(gtime_t), intent(in) :: time_in, teph
integer*4, intent(in) :: sat, iode
type(nav_t), intent(in) :: nav
real*8, intent(out) :: rs(6), dts(2), var
integer*4, intent(out) :: svh, stat
type(gtime_t) time
type(eph_t), pointer :: eph
type(geph_t), pointer :: geph
real*8 :: rst(3),dtst(1),tt
integer*4 i,sys,prn
tt=1d-3
call satsys(sat,prn,sys)
time=time_in
svh=-1
if(sys==SYS_GPS .or. sys==SYS_GAL .or. sys==SYS_QZS .or. sys==SYS_CMP)then
    eph=>seleph(teph,sat,iode,nav)
    if (.not.associated(eph))then
        stat=0; return
    endif
    call eph2pos(time,eph,rs(1:3),dts(1),var)
    time=timeadd(time,tt)
    call eph2pos(time,eph,rst(1:3),dtst(1),var)
    svh=eph%svh
elseif(sys==SYS_GLO)then
    geph=>selgeph(teph,sat,iode,nav)
    if (.not.associated(geph))then
        stat=0; return
    endif
    call geph2pos(time,geph,rs(1:3),dts(1),var)
    time=timeadd(time,tt)
    call geph2pos(time,geph,rst(1:3),dtst(1),var)
    svh=geph%svh
else
    stat=0; return
endif
! satellite velocity and clock drift by differential approx 
rs(4:6)=(rst(1:3)-rs(1:3))/tt
dts(2)=(dtst(1)-dts(1))/tt
stat=1
end subroutine

! satellite position and clock ------------------------------------------------
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
!-----------------------------------------------------------------------------
subroutine satpos(time, teph, sat, ephopt, nav, rs, dts, var, svh, stat)
implicit none
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

! satellite positions and clocks ----------------------------------------------
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
!-----------------------------------------------------------------------------
subroutine satposs(teph, obs, n, nav, ephopt, rs, dts, var, svh)
implicit none
type(gtime_t), intent(in) :: teph
type(obsd_t), intent(in) :: obs(n)
integer*4, intent(in) :: n, ephopt
type(nav_t), intent(in) :: nav
real*8, intent(out) :: rs(n,6), dts(n,2), var(n)
integer*4, intent(out) :: svh(n)
type(gtime_t) time(2*MAXOBS)
real*8 dt,pr
integer*4 i,j,ulimit,info

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
end module ephe_f90_