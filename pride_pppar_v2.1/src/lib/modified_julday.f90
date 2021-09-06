!
!! modified_julday.f90
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
!!!
!! purpose  : return modified julian day given year, month and day
!!
!! parameter: iyear -- year 4-digits
!!            imonth -- 1 to 12 or
!!                      0   indicates iday is day of year
!!            iday   -- 1 to 31 or
!!                      day of year, if imonth is zero
!!
integer*4 function modified_julday(iday, imonth, iyear)
  implicit none

  integer*4 imonth, iday, iyear
!
!! local
  integer*4 iyr, doy_of_month(12)
  data doy_of_month/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
!
!! check the input data
  if ((iyear .lt. 0 .or. imonth .lt. 0 .or. iday .lt. 0 .or. imonth .gt. 12 .or. iday .gt. 366) &
      .or. (imonth .ne. 0 .and. iday .gt. 31)) then
    write (*, '(a)') '***ERROR(modified_julday): incorrect date(year,month,day): ', &
      iyear, imonth, iday
    call exit(1)
  endif

  iyr = iyear
  if (imonth .le. 2) iyr = iyr - 1
  modified_julday = 365*iyear - 678941 + iyr/4 - iyr/100 + iyr/400 + iday
  if (imonth .ne. 0) modified_julday = modified_julday + doy_of_month(imonth)

  return
end
