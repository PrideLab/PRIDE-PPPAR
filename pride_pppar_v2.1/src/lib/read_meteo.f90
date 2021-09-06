!
!! read_meteo.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : read meterology information
!! parameter:
!!    input : jd,sod -- requested time
!!            map  -- map options
!!            geod -- ellipsoid coordinates
!!    output: pres   -- pressure, mbar
!!            temp   -- temperature, C
!!            humd   -- humidity (0.0-1.0)
!
subroutine read_meteo(jd, sod, imet, map, geod, pres, temp, humd, undu)
  implicit none

  integer*4 jd, imet
  real*8 sod, geod(3), pres, temp, humd, undu
  character*(*) map
!
!! local
  real*8 V(10, 10), W(10, 10)
!
!! VMF1 and VMF3 do not need meteorological information
  if (map(1:3) .eq. 'VM1' .or. map(1:3) .eq. 'VM3') return
!
!! initialization
  V = 0.d0
  W = 0.d0
!
!! compute 9x9 spherical harmonics and geoid undulation
  call spherical_harmonics(geod(1), geod(2), V, W, undu)
!
!! meteorology items
  if (imet .eq. 0) then
    pres = 1013.25d0
    temp = 20.d0
    humd = 0.5d0
  else if (imet .ne. 0) then

  endif

  if (geod(3) .lt. 0.d0 .or. geod(3) .gt. 9.d0) return
!! get subsequent meteorology info at station
  if (map(1:3) .eq. 'GMF') then
    call global_meteo(jd*1.d0, geod(3)*1.d3, V, W, undu, pres, temp, humd)
  else
    call mete_sea2site(.true., pres, temp, humd, geod(3)*1.d3 - undu, pres, temp, humd)
  endif

  return
end
