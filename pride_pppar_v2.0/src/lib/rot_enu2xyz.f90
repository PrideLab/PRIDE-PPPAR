!
!! rot_enu2xyz.f90
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
!! Contributor: Maorong Ge
!! 
!!
!!
subroutine rot_enu2xyz(lat, lon, rotmat)
! purpose  : rotation matrix from local station system to earth-fixed system.
!            Local system is east-north-radial/up  right-hand system.
!
! parameter:  lat lon -- latitude and longitode in radians
!             rotmat  -- rotation matrix
!

  implicit none
  real*8 lat, lon, rotmat(3, 3)

  real*8 coslat, sinlat, coslon, sinlon, pi

  pi = 4.d0*datan(1.d0)

! pi/2.d0-lat clockwise rotating around east to Z
  coslat = dcos(lat - pi/2.d0)
  sinlat = dsin(lat - pi/2.d0)

! pi/2.d0+lon clockwise rotating around up(Z)
  coslon = dcos(-pi/2.d0 - lon)
  sinlon = dsin(-pi/2.d0 - lon)

  rotmat(1, 1) = coslon
  rotmat(1, 2) = sinlon*coslat
  rotmat(1, 3) = sinlon*sinlat
  rotmat(2, 1) = -sinlon
  rotmat(2, 2) = coslon*coslat
  rotmat(2, 3) = coslon*sinlat
  rotmat(3, 1) = 0.d0
  rotmat(3, 2) = -sinlat
  rotmat(3, 3) = coslat

  return
end
