!
!! phase_windup.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
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
!! purpose  : phase wind-up correction ( the receiver and satellite antenna orientation dependent
!!            phase corrections). See Wu J.T., et al., Manuscripta Geogetica (1993) 18, pp91-98
!!
!! parameter:
!!           first -- first call for this satellite-station pair
!!           rot_f2j -- rotation matrix from earth-fixed to inertial (J2000)
!!           rot_l2f -- rotation matrix from station system (enu right-hand system) to earth-fixed
!!           xbf,ybf,zbf -- unit vectors of spacecraft-fixed system (rotation matrix from body-fixed
!!                          to inertial)
!!           xrec2sat -- vector from rec. to satellite  for k-direction
!!           dphi0 -- initial dphi
!!           dphi  --  phase correction
!!
!
subroutine phase_windup(first, rot_f2j, rot_l2f, xbf, ybf, zbf, xrec2sat, dphi0, dphi)
  implicit none
  logical*1 first
  real*8 rot_f2j(3, 3), rot_l2f(3, 3), xbf(1:*), ybf(1:*), zbf(1:*), xrec2sat(1:*)
  real*8 dphi0, dphi
!
!! local
  logical*1 lfirst
  integer*4 j, n
  real*8 x_l(3), y_l(3), x_f(3), y_f(3), x_j(3), y_j(3), dummy(3), rlength, kusi
  real*8 k(3), d_r(3), d_s(3), twopi
!
!! function called
  real*8 dot

  save lfirst, twopi
!
!! enu (right-hand system)
  data lfirst/.true./, x_l/0.0, 1.0, 0.0/, y_l/-1.0, 0.0, 0.0/
!
!! initial value set to zero
  if (lfirst) then
    lfirst = .false.
    twopi = 8.d0*datan(1.d0)
  endif
!
!! unit vector of the signal transmitting direction
  call unit_vector(3, xrec2sat, k, rlength)
  do j = 1, 3
    k(j) = -k(j)
  enddo
!
!! equivalent antenna dipole for both receiver and satellite antenna
!!    D = x_j - k (k . x_j) - k x y_j
!! k is the unit vector of signal transmitting dirrection
!
!! transfer the unit vector of the antenna dipole unit in local/antenna system to
!! inertial system.
!! For ground receiver local   => earth fixed      => inertial
!! For LEO    receiver antenna => spacecraft fixed => inertial
  call matmpy(rot_l2f, x_l, x_f, 3, 3, 1)
  call matmpy(rot_l2f, y_l, y_f, 3, 3, 1)
  call matmpy(rot_f2j, x_f, x_j, 3, 3, 1)
  call matmpy(rot_f2j, y_f, y_j, 3, 3, 1)
!
!! D = x_j - k(k . x_j) + k x y_j
  rlength = dot(3, k, x_j)
  call cross(k, y_j, dummy)
  do j = 1, 3
    d_r(j) = x_j(j) - rlength*k(j) + dummy(j)
  enddo
  call unit_vector(3, d_r, d_r, rlength)
!
!! the same for the satellite antenna
  rlength = dot(3, k, xbf)
  call cross(k, ybf, dummy)
  do j = 1, 3
    d_s(j) = xbf(j) - rlength*k(j) - dummy(j)
  enddo
  call unit_vector(3, d_s, d_s, rlength)
!
!! kusi  k . (D_s x D_r)
  call cross(d_s, d_r, dummy)
  kusi = dot(3, k, dummy)
  dphi = dot(3, d_s, d_r)
  if (dabs(dphi) .gt. 1.d0) dphi = nint(dphi)*1.d0
  dphi = dsign(1.d0, kusi)*dacos(dphi)/twopi
!      if(kusi.lt.0.d0) dphi=-dphi
!
!! interger part
  if (first) then
    first = .false.
    n = 0
  else
    n = nint(dphi0 - dphi)
  endif
  dphi = n + dphi
!
!! save for the next epoch
  dphi0 = dphi

  return
end
