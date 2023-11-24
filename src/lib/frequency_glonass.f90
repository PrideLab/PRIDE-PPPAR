!
!! frequency_glonass.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
!! 
!!
!
subroutine frequency_glonass(FREQ1_R,FREQ2_R)
implicit none
include '../header/const.h'
real*8 :: FREQ1_R(-50:50),FREQ2_R(-50:50)
integer*4 ks
real*8 g_R,lambda1_R(-50:50),lambda2_R(-50:50),lambdaw_R(-50:50),c1_R,c2_R
do ks=-50,50
  FREQ1_R(ks)=1.602d9+0.0005625d9*ks
  FREQ2_R(ks)=FREQ1_R(ks)*7.d0/9.d0
end do
end
