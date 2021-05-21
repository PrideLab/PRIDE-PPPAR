!
!! lower_string.f90
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
!!
!! purpose  : return the lower case of the input string
!!
!! parameter: string -- input string
!!            lower_string -- return string
!!
!!
character*(*) function lower_string(string)
  implicit none
  character*(*) string
!
!! local
  integer*4 i, length, iadd
!
!! ascii code difference between upper case and lower case
  iadd = ichar('a') - ichar('A')

  length = min(len(string), len(lower_string))
  do i = 1, length
    if (lge(string(i:i), 'A') .and. lle(string(i:i), 'Z')) then
      lower_string(i:i) = char(ichar(string(i:i)) + iadd)
    else
      lower_string(i:i) = string(i:i)
    endif
  enddo
!
!! put SPACE for the rest
  do i = length + 1, len(lower_string)
    lower_string(i:i) = ' '
  enddo

  return
end
