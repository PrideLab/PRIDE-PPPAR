!
!! pointer_string.f90
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
!! purpose  : return the FIRST position of an element appears in a string
!!            array.
!!
!! parameter: n -- number of elements in the array
!!            string_array --  string array
!!            string    -- the element
!!            point_int -- position of the element in the array.
!!                         zere is it is not in the array
!!
!!
integer*4 function pointer_string(n, string_array, string)
  implicit none
  integer*4 n
  character*(*) string_array(1:*), string
!
!! local
  integer*4 i, l

  pointer_string = 0
!
!! if input string is empty
  l = len_trim(string)
  if (l .eq. 0) return

  do i = 1, n
    if (index(string_array(i), string(1:l)) .eq. 1) then
      pointer_string = i
      exit
    endif
  enddo
  return
end
