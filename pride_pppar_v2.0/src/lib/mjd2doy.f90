!
!! mjd2doy.f90
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
!!!
!! purpose  : get year and day of year given modified julian day
!!            this algorithm is a little bit slower but validates for long time
!!
!! parameter: jd  -- modified julian day
!!            iyear, idoy -- year and day of year
!!
!!

subroutine mjd2doy(jd, iyear, idoy)
  implicit none
  integer*4 jd, iyear, idoy
!
!! local
  integer*4 modified_julday

  iyear = (jd + 678940)/365
  idoy = jd - modified_julday(1, 1, iyear)
  do while (idoy .le. 0)
    iyear = iyear - 1
    idoy = jd - modified_julday(1, 1, iyear) + 1
  enddo

  return
end
