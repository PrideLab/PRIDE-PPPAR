!
!! glsbrd.f90
!!
!!    Copyright (C) 2023 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program. If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Jianghui Geng, Songfeng Yang
!! refer to RTKLIB
!! 
!!
!!
!! Calculate Glonass position, velocity and clocks from broadcast ephemeris
!
subroutine glsbrd(cmode,neph,ephem,iprn,jeph,jd,sod,xsat,v,dtsat,dtref)
implicit none
include '../header/brdeph.h'

character*3 cmode
integer*4 neph
type(brdeph) ephem(1:*)
integer*4 jd,jeph
character*3 iprn
real*8 sod,xsat(1:*),v,dtsat,dtref
!
!! local
integer*4 i
real*8 dt,tt,x(6)
!
!! function used
real*8 timdif
!
!! find out the nearest ephem.
dtref=0.d0
if(jeph.eq.0) then
  dtref=3600.d0    ! seconds
  do i=1,neph
    if(ephem(i)%svn.eq.iprn) then
      dt=timdif(jd,sod,ephem(i)%jd,ephem(i)%sod)
      if(dabs(dt).le.dtref) then
        dtref=dabs(dt) ! in seconds
        jeph =i
      endif
    endif
  enddo
endif
dtref=dtref/3600.d0 ! seconds to hours
if(jeph.eq.0) return
!
!! dt 
dt=timdif(jd,sod,ephem(jeph)%jd,ephem(jeph)%sod) ! seconds
!
!! satellite clock correction
if(cmode(3:3).eq.'y') dtsat=-ephem(jeph)%a0+ephem(jeph)%a1*dt
if(cmode(1:2).eq.'nn') return
!
!! orbit integration
x(1:3)=ephem(jeph)%pos(1:3)
x(4:6)=ephem(jeph)%vel(1:3)
if(dt.lt.0.d0) then
  tt=-30.d0   ! integration step in seconds
else
  tt= 30.d0
endif
do while(dabs(dt).gt.0.d0)
  if(dabs(dt).lt.30.d0) tt=dt
  call glorbit(tt,x,ephem(jeph)%acc)
  dt=dt-tt
enddo
!
!! transformed to WGS84
call pz2wgs84(jd,sod,x,xsat) ! position and velocity in meter

v=0.d0 !ysf
return
end
!
!! glonass orbit differential equations
subroutine deq(x,xdot,acc)
implicit none
include '../header/const.h'

real*8 x(1:*),xdot(1:*),acc(1:*)
!
!! local
real*8 a,b,c,r2,r3,omg2,gme,wearth,j2
!
!! function used
real*8 dot

data gme,wearth,j2/3.986005d14,7.2921151467d-5,1.0826257d-3/
save gme,wearth,j2

r2  =dot(3,x,x)
r3  =r2*dsqrt(r2)
if(r2.le.0.d0) then
  xdot(1:6)=0.d0
  return
endif
omg2=wearth**2
a=1.5d0*j2*gme*erad**2*1.d6/r2/r3 ! 3/2*J2*mu*Ae^2/r^5
b=5.d0*x(3)**2/r2                 ! 5*z^2/r^2
c=-gme/r3-a*(1.d0-b)              ! -mu/r^3-a(1-b)
xdot(1)=x(4)
xdot(2)=x(5)
xdot(3)=x(6)
xdot(4)=(c+omg2)*x(1)+2.d0*wearth*x(5)+acc(1)
xdot(5)=(c+omg2)*x(2)-2.d0*wearth*x(4)+acc(2)
xdot(6)=(c-2.d0*a)*x(3)+acc(3)

return
end
!
!! glonass position and velocity by numerical integration
subroutine glorbit(t,x,acc)
implicit none

real*8 t,x(1:*),acc(1:*)
!
!! local
real*8 k1(6),k2(6),k3(6),k4(6),w(6)

call deq(x,k1,acc)
w(1:6)=x(1:6)+k1(1:6)*t/2.d0
call deq(w,k2,acc)
w(1:6)=x(1:6)+k2(1:6)*t/2.d0
call deq(w,k3,acc)
w(1:6)=x(1:6)+k3(1:6)*t
call deq(w,k4,acc)
x(1:6)=x(1:6)+(k1(1:6)+2.d0*k2(1:6)+2.d0*k3(1:6)+k4(1:6))*t/6.d0

return
end
