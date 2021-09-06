!
!! sod2hms.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
!! 
!!
!
!! purpose   : seconds of day to hour minute and second
!!
!
subroutine sod2hms(sod, ih, im, sec)
  implicit none
  integer*4 ih, im
  real*8 sod, sec

  ih = int(sod/3600.d0)
  im = int((sod - ih*3600.d0)/60.d0)
  sec = sod - ih*3600.d0 - im*60.d0

  return
end
