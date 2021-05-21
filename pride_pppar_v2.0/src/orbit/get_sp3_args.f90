!
!! get_sp3orb_args.f90
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
!! Contributor: Maorong Ge, Jianghui Geng
!! 
!!
!
subroutine get_sp3orb_args(sp3fil1, sp3fil2, sp3fil3, mersp3fil)
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'

  character*256 sp3fil1      ! last day
  character*256 sp3fil2      ! processing day
  character*256 sp3fil3      ! next day
  character*256 mersp3fil

!
!! functions called
  integer*4 get_valid_unit
!
  sp3fil1 = ' '
  sp3fil2 = ' '
  sp3fil3 = ' '
  mersp3fil = ' '
  call getarg(1, sp3fil1)        ! last day
  call getarg(2, sp3fil2)        ! processing day
  call getarg(3, sp3fil3)        ! next day
  call getarg(4, mersp3fil)
  return
end
