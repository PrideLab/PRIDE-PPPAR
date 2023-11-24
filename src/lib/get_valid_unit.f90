!
!! get_valid_unit.f90
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
!! purpose  : get a valid file unit for open a new file.
!!
!! parameter: iunit0 -- start unit
!!             get_valid_unit -- output unit
!!
!!
!! note     : should check the max. unit number of the system
!!
!
integer*4 function get_valid_unit(iunit0)
  implicit none
  integer*4 iunit0
!
!! local
  integer*4 iunit
  character*256 string

  get_valid_unit = 0
  iunit = iunit0

  if (iunit0 .le. 10) iunit = 10
  do while (get_valid_unit .eq. 0)
    string = ""
    inquire (unit=iunit, name=string)
    if (len_trim(string) .ne. 0) then
      iunit = iunit + 1
    else
      get_valid_unit = iunit
    endif
  enddo

  return
end
