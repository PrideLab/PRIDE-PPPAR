!
!! yeardoy2monthday.f90
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
!!! purpose   : get month and day given year and day of year
!! parameters: iyear,idoy -- year and day of year
!!             imonth,iday -- month and day in the yar
!!
!!
subroutine yeardoy2monthday(iyear, idoy, imonth, iday)
  implicit none
  integer*4 iyear, idoy, imonth, iday
!
!! local
  integer*4 days_in_month(12), id

  data days_in_month/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

  days_in_month(2) = 28
  if (mod(iyear, 4) .eq. 0 .and. (mod(iyear, 100) .ne. 0 .or. mod(iyear, 400) .eq. 0)) days_in_month(2) = 29
  id = idoy
  do imonth = 1, 12
    id = id - days_in_month(imonth)
    if (id .gt. 0) cycle
    iday = id + days_in_month(imonth)
    exit
  enddo
  return
end
