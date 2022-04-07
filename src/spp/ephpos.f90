
! Calculate satellite position ----------------------------------------------
! satellite position and clock by broadcast ephemeris -----------------------
subroutine ephpos(time_in, teph, sat, nav, iode, rs, dts, var, svh, stat)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time_in, teph
integer*4, intent(in) :: sat, iode
type(nav_t), intent(in) :: nav
real*8, intent(out) :: rs(6), dts(2), var
integer*4, intent(out) :: svh, stat
type(gtime_t) time,timeadd
type(eph_t), pointer :: eph  !,seleph
type(geph_t), pointer :: geph  !,selgeph
real*8 :: rst(3),dtst(1),tt,var_uraeph
external :: timeadd,var_uraeph  !,seleph,selgeph

interface
  function seleph(time, sat, iode, nav)
  implicit none
  include 'file_para.h'
  type(gtime_t), intent(in) :: time
  integer*4, intent(in) :: sat,iode
  type(nav_t), intent(in) :: nav
  type(eph_t), pointer :: seleph
  end function

  function selgeph(time, sat, iode, nav)
  implicit none
  include 'file_para.h'
  type(gtime_t), intent(in) :: time
  integer*4, intent(in) :: sat,iode
  type(nav_t), intent(in) :: nav
  type(geph_t), pointer :: selgeph
  end function
end interface

integer*4 i,sys,prn
tt=1d-3
call satsys(sat,prn,sys)
time=time_in
svh=-1
if(sys==SYS_GPS .or. sys==SYS_GAL .or. sys==SYS_QZS .or. sys==SYS_CMP)then
    !eph=seleph(teph,sat,iode,nav)
    eph=>seleph(teph,sat,iode,nav)
    if (.not.associated(eph))then
        stat=0; return
    endif
    call eph2pos(time,eph,rs(1:3),dts(1),var)
    time=timeadd(time,tt)
    call eph2pos(time,eph,rst(1:3),dtst(1),var)
    svh=eph%svh
elseif(sys==SYS_GLO)then
    !geph=selgeph(teph,sat,iode,nav)
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

! broadcast ephemeris to satellite position and clock bias ------------------
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
!----------------------------------------------------------------------------
subroutine eph2pos(time, eph, rs, dts, var)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time
type(eph_t), intent(in) :: eph
real*8, intent(out) :: rs(3),dts,var
real*8 tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge1
real*8 xg,yg,zg,sino1,coso1,timediff,SQR,var_uraeph
external :: timediff,SQR,var_uraeph
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

! glonass ephemeris to satellite position and clock bias --------------------
! compute satellite position and clock bias with glonass ephemeris
! args   : gtime_t time     I   time (gpst)
!          geph_t *geph     I   glonass ephemeris
!          real*8 *rs       O   satellite position {x,y,z} (ecef) (m)
!          real*8 *dts      O   satellite clock bias (s)
!          real*8 *var      O   satellite position and clock variance (m^2)
! return : none
! notes  : see ref (2)
!----------------------------------------------------------------------------
subroutine geph2pos(time, geph, rs, dts, var)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time
type(geph_t), intent(in) :: geph
real*8, intent(out) :: rs(3),dts,var
real*8 t,tt,x(6),timediff,rcond,SQR
external :: timediff,rcond,SQR
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

!! glonass position and velocity by numerical integration -------------------
subroutine glorbit(t, x, acc)
implicit none
include 'file_para.h'
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

!! glonass orbit differential equations -------------------------------------
subroutine deq(x, xdot, acc)
implicit none
include 'file_para.h'
real*8, intent(in) :: x(6),acc(3)
real*8, intent(out) :: xdot(6)
real*8 a,b,c,r2,r3,omg2,dot,SQR
external :: dot,SQR
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
