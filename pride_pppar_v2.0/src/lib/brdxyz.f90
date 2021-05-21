!
!! brdxyz.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License (version 3) as
!!    published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose   : compute position and velocity of GNSS satellite using broadcast ephemeris
!!
!! parameters:
!!     cmode  == string with 3 characters for pos,vel and clock, respectively.
!!             for example,   yyy -- all required, nny -- only clock
!!     lstick == use the last one (jeph) if it is not zero, to avoid discontinuity by two brdeph
!!     neph, ephem == total number of brdeph in 'ephem'
!!     iprn == PRN of the satellite its pos and vel to be computed
!!     jd, sod ==  epoch time, (#week, sow) or (mjd,sod) GPS time
!!     xsat, dtsat == position , velocity and clock correction
!!     v ==
!!     dtref == t - toe/toc
!!
!
subroutine brdxyz(cmode, lstick, neph, ephem, iprn, jeph, jd, sod, xsat, v, dtsat, dtref)
  implicit none
  include '../header/brdeph.h'

  character*3 cmode
  logical*1 lstick
  integer*4 neph
  type(brdeph) ephem(1:*)
  integer*4 jd, jeph
  character*3 iprn
  real*8 sod, xsat(1:*), dtsat, dtref, v
!
!! local variables
  integer*4 i, j, ieph
  real*8 gme, wearth
  real*8 pi, dt, a, xn, xm, ex, e, v0, vs, vc, phi, ccc, sss, du, dr, di, r, u, xi, xx, yy, &
    xnode, term, xpdot, ypdot, asc, xinc, xp, yp, asctrm,d_ex,d_v,d_u,d_r,d_i,d_omega
  real*8 xg,yg,zg,cos5,sin5
!
!! Angular velocity of the earth  and  Earth's gravitational constant
  data gme, wearth,cos5,sin5/3.986005d14, 7.2921151467d-5,0.9961946980917456,-0.0871557427476582/
  pi = datan(1.d0)*4.d0
  if(iprn(1:1).eq."G".or.iprn(1:1).eq.'R') then
    gme=3.986005d14
  elseif(iprn(1:1).eq."E") then
    gme=3.986004418d14
  elseif(iprn(1:1).eq."C") then
    gme=3.986004418d14
  elseif(iprn(1:1).eq."J") then
    gme=3.986005d14
  endif
!
!! find out the nearest ephem.
  dtref=0.d0
  if (jeph .eq. 0) then
    dtref = 12.d0    ! hours
    do ieph = 1, neph
      if (ephem(ieph)%svn .eq. iprn) then
        dt = (jd - ephem(ieph)%jd)*24.d0 + (sod - ephem(ieph)%sod)/3600.d0
        if (dabs(dt) .le. dtref) then
          dtref = dabs(dt)
          jeph = ieph
          if (dtref .le. 1.d0) exit
        endif
      endif
    enddo
  endif
  if (jeph .eq. 0) return
!
!! dt
  i = jeph
  dt = (jd - ephem(jeph)%jd)*86400.d0 + sod - ephem(jeph)%sod
!
!! satellite clock correction
  if (cmode(3:3) .eq. 'y') dtsat = ephem(i)%a0 + (ephem(i)%a1 + ephem(i)%a2*dt)*dt

  if (cmode(1:2) .eq. 'nn') return
!
!! compute sat. coordinate in wgs-84
  a = ephem(i)%roota**2
  xn = dsqrt(gme/a/a/a)
  xn = xn + ephem(i)%dn
!
!! iterate to solve EX
  xm = ephem(i)%m0 + xn*dt
  ex = xm
  e = ephem(i)%e
  do j = 1, 12
    ex = xm + e*dsin(ex)
  enddo
!
!! determination of V
  v0 = 1.d0 - e*dcos(ex)
  vs = dsqrt(1.d0 - e*e)*dsin(ex)/v0
  vc = (dcos(ex) - e)/v0
  v = datan2(vs, vc)
  phi = v + ephem(i)%omega

  ccc = dcos(2*phi)
  sss = dsin(2*phi)
  du = ephem(i)%cuc*ccc + ephem(i)%cus*sss
  dr = ephem(i)%crc*ccc + ephem(i)%crs*sss
  di = ephem(i)%cic*ccc + ephem(i)%cis*sss
  r = a*(1 - e*dcos(ex)) + dr
  u = phi + du

  xi = ephem(i)%i0 + di + ephem(i)%idot*dt
  xx = r*dcos(u)
  yy = r*dsin(u)
  if(iprn.eq."C01".or.iprn.eq."C02".or.iprn.eq."C03".or.iprn.eq."C04".or.iprn.eq."C05") then !!for 1-5 beidou GEO satellite
    xnode=ephem(i)%omega0+ephem(i)%omegadot*dt-wearth*ephem(i)%toe
    xg=xx*dcos(xnode)-yy*dcos(xi)*dsin(xnode)
    yg=xx*dsin(xnode)+yy*dcos(xi)*dcos(xnode)
    zg=yy*dsin(xi)
    xsat(1)=xg*dcos(wearth*dt)+yg*dsin(wearth*dt)*cos5+zg*dsin(wearth*dt)*sin5
    xsat(2)=-xg*dsin(wearth*dt)+yg*dcos(wearth*dt)*cos5+zg*dcos(wearth*dt)*sin5
    xsat(3)=-yg*sin5+zg*cos5
  else
    xnode = ephem(i)%omega0 + (ephem(i)%omegadot - wearth)*dt
    xnode = xnode - wearth*ephem(i)%toe
    xsat(1) = xx*dcos(xnode) - yy*dcos(xi)*dsin(xnode)
    xsat(2) = xx*dsin(xnode) + yy*dcos(xi)*dcos(xnode)
    xsat(3) = yy*dsin(xi)
  endif

  if (cmode(2:2) .eq. 'n') return
!
!! velocity
  d_ex=xn/(1.d0-e*dcos(ex))
  d_v=dsin(ex)*d_ex*(1.d0+e*dcos(v))/(1.d0-e*dcos(ex))/dsin(v)
  d_u=d_v+2.d0*(ephem(i)%cus*ccc-ephem(i)%cuc*sss)*d_v 
  d_r=a*e*dsin(ex)*d_ex+2.d0*(ephem(i)%crs*ccc-ephem(i)%crc*sss)*d_v
  d_i=ephem(i)%idot+2.d0*(ephem(i)%cis*ccc-ephem(i)%cic*sss)*d_v

  d_omega=(ephem(i)%omegadot-wearth)

  xpdot=d_r*dcos(u)-r*dsin(u)*d_u
  ypdot=d_r*dsin(u)+r*dcos(u)*d_u

  xsat(4)=xpdot*dcos(xnode)-ypdot*dcos(xi)*dsin(xnode)+yy*dsin(xi)*d_i*dsin(xnode)-xsat(2)*d_omega
  xsat(5)=xpdot*dsin(xnode)+ypdot*dcos(xi)*dcos(xnode)-yy*dsin(xi)*d_i*dcos(xnode)+xsat(1)*d_omega
  xsat(6)=ypdot*dsin(xi)+yy*dcos(xi)*d_i
  return

end
