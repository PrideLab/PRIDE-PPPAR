!
!! left_justify_string.f90
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
!! purpose  : left just a string, remove space at the beginning of a string
!!
!! parameter: string -- input and ouput string
!!
!!
subroutine left_justify_string(string)
  implicit none
  character*(*) string

!
!! local
  integer*4 i, j, n
  n = len(string)
  do i = 1, n
    if (string(i:i) .ne. ' ') exit
  enddo
  if (i .le. n) string = string(i:n)
!
!! empty the rest
  j = n - i + 1
  do i = j + 1, n
    string(i:i) = ' '
  enddo

  return
end

