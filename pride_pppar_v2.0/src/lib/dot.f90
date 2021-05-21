!
!! dot.f90
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
!! purpose  : dot product of twe vectors
!!
!! parameter: n -- dimmesion of vectors
!!             v1, v2 -- input vecotrs
!!             dot -- dot product of v1 and v2
!!
!

real*8 function dot(n, v1, v2)
  implicit none
  integer*4 n
  real*8 v1(1:*), v2(1:*)
!
!! local
  integer*4 i

  dot = 0.d0
  do i = 1, n
    dot = dot + v1(i)*v2(i)
  enddo
  return
end
