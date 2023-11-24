!
!! gpsweek.f90
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
!! Contributor: Maorong Ge
!! 
!!
!!
!! purpose   : get GPS week number and day of week given day month and year
!!             or doy and  year
!!
!! parameters: day month year -- if month==0, day = doy
!!             week, wd       -- GPS week and day of week
!!
!
subroutine gpsweek(day, month, year, week, wd)
  implicit none
  integer*4 year, month, day, week, wd
  real*8 sod, sow
!
!! local
  integer*4 mjd
!
!! function called
  integer*4 modified_julday

  if (year .ne. 0) then
    if (month .ne. 0) then
      mjd = modified_julday(day, month, year)
    else
      mjd = modified_julday(1, 1, year) + day - 1
    endif
  else
    mjd = day
  endif

  week = (mjd - 44244)/7
  wd = mjd - 44244 - week*7
  return
end
