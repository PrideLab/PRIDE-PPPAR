!
!! day_of_year.f90
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
!!    purpose   : given year,month and day, return day of the year
!!
!
integer*4 function day_of_year(iday, imonth, iyear)
  integer*4 iday, imonth, iyear
!
!! function called
  integer*4 modified_julday

  day_of_year = modified_julday(iday, imonth, iyear) - &
                modified_julday(1, 1, iyear) + 1
  return
end
