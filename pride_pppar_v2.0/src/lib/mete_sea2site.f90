!!
!! mete_sea2site.f90
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
!
!! purpose  : conversion meteo-data on see level to station height.
!!
!! parameter:
!!    input:  ts,ps,hs --- temperature, pressure, humidity(%) on sea level
!!              height --- distance from the ellipesoid to the station (m)
!!
!!   output:     t,p,h --- temperature,pressure,humidity(%) at height
!!
!!
subroutine mete_sea2site(lsea, ps, ts, hs, height, p, t, h)
  implicit none

  logical*1 lsea
  real*8 ts, ps, hs, height, t, p, h

  if (height .gt. -1.d4 .and. height .lt. 2.d4) then
    if (lsea) then
      p = ps*(1.d0 - 2.26d-5*height)**5.225d0
      t = ts - height*0.0065d0
      h = hs*exp(-6.396d-4*height)
    else
      p = ps/(1.d0 - 2.26d-5*height)**5.225d0
      t = ts + height*0.0065d0
      h = hs/exp(-6.396d-4*height)
    endif
  else
    p = ps
    t = ts
    h = hs
  endif

  return
end
