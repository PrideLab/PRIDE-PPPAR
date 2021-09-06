!
!! slerp.f90
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
!! Contributor: Shuyin Mao
!
subroutine slerp(qs,qe,t,q)
implicit none
!
!! local
integer*4 i
real*8 qs(4),qe(4),q(4),t
real*8 sina,cosa,k0,k1,a

cosa = qs(1)*qe(1) + qs(2)*qe(2) +qs(3)*qe(3) + qs(4)*qe(4);
! If the dot product is negative, the quaternions have opposite 
! handed-ness and slerp won't take
! the shorter path. Fix by reversing one quaternion.
if(cosa.lt.0) then 
  qe(1:4) = -1.d0*qe(1:4)
  cosa = -1.d0*cosa
endif
! If the inputs are too close for comfort, linearly interpolate
if(cosa.gt.0.9995d0) then
  k0 = 1.0d0-t;
  k1 = t;
else 
  sina = dsqrt(1.d0-cosa*cosa)
  a = datan2(sina,cosa)
  k0=dsin((1.d0-t)*a)/sina
  k1=dsin(t*a)/sina;
endif
do i=1,4
  q(i) = qs(i)*k0 + qe(i)*k1
enddo
return
end
