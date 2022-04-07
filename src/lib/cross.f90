!
!! cross.f90
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
!! purpose  :  cross product of two 3-dimmension vectors
!!
!! parameter:  v1,v2 -- 2 input 3-d vector
!!             vout  -- cross product of v1 and v2
!!
!
subroutine cross(v1, v2, vout)
  implicit none
  real*8 v1(1:*), v2(1:*), vout(1:*)
!
!! local
  real*8 tmp(3)

  tmp(1) = v1(2)*v2(3) - v1(3)*v2(2)
  tmp(2) = v1(3)*v2(1) - v1(1)*v2(3)
  tmp(3) = v1(1)*v2(2) - v1(2)*v2(1)

  vout(1) = tmp(1)
  vout(2) = tmp(2)
  vout(3) = tmp(3)

  return
end
