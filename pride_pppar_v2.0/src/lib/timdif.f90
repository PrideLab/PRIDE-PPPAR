!
!! timdif.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose   : get time difference
!!
!! parameter :
!!    input  : jd2,sod2 -- time to be differenced
!!             jd1,sod1 -- time difference
!!
!!
!
real*8 function timdif(jd2, sod2, jd1, sod1)
  implicit none

  integer*4 jd1, jd2
  real*8 sod1, sod2
!
!! local
  real*8 spd

  data spd/86400.d0/

  timdif = spd*(jd2 - jd1) + sod2 - sod1

  return
end
