!
!! xyzblh.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!!! purpose  :  computation of ellipsoidal coordinates rb,rl,rh
!!             given the cartesian coordinates x,y,z
!!
!! parameter:    x,y,z:    cartesian coord. of the point in metres
!!               scale:    scale to meter
!!            dx,dy,dz:    translation components from the origin of
!!                         the cart. coord. system (x,y,z) to the
!!                         center of the ref.ellipsoid (in metres)
!!                 a,b:    semi-major and semi-minor axes of the
!!                         ref.ellipsoid in metres
!!                  rb:    ell. latitude (arc)
!!                  rl:    ell. longitude (pos. east of greenwich)
!!                  rh:    ell. height (m)
!!
!
subroutine xyzblh(x, scale, a0, b0, dx, dy, dz, geod)
  implicit none
  real*8 x(1:*), scale, a0, b0, dx, dy, dz, geod(1:*)
!
!! local variable
  real*8 a, b, xp, yp, zp, e2, s, zps, n, hp, rbi, rbd, rhd, pi
!
!! default ellipsoid WGS84
  if (a0 .eq. 0.d0) then
    a = 6378137.d0
    b = 298.257223563d0
  else
    a = a0
    b = b0
  endif
  pi = 4.d0*datan(1.d0)
!
!! to meters
  xp = x(1)*scale + dx
  yp = x(2)*scale + dy
  zp = x(3)*scale + dz
!
!! semi axis or finv
  if (b .le. 6000000.d0) b = a - a/b
  e2 = (a*a - b*b)/(a*a)

  s = dsqrt(xp*xp + yp*yp)
  geod(2) = datan2(yp, xp)
  if (geod(2) .lt. 0.d0) geod(2) = 2.d0*pi + geod(2)
  zps = zp/s
  geod(3) = dsqrt(xp*xp + yp*yp + zp*zp) - a
  geod(1) = datan(zps/(1.d0 - e2*a/(a + geod(3))))
!
10 n = a/dsqrt(1.d0 - e2*dsin(geod(1))**2)
  hp = geod(3)
  rbi = geod(1)
  geod(3) = s/dcos(geod(1)) - n
  geod(1) = datan(zps/(1.d0 - e2*n/(n + geod(3))))
  rbd = dabs(rbi - geod(1))
  rhd = dabs(hp - geod(3))
  if (rbd*n .gt. 1.d-4 .or. rhd .gt. 1.d-4) goto 10

  return
end
