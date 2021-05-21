!
!! find_flag.f90
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
!!!
!! purpose  : find char in flg and return it in iepo
!!            reverse search is also allowed
!! parameter:
!!    input : istart,istop -- beginning and ending tag
!!            flg -- flag array
!!            char -- searching target
!!    output: iepo -- epoch number
!!
subroutine find_flag(istart, istop, flg, char, iepo)
  implicit none
  integer*4 flg(1:*), iepo, istart, istop
  character*(*) char
! local
  integer*4 i, istep
! function called
  logical*1 istrue

  istep = 1
  if (istop .le. istart) istep = -1

  iepo = -1
  do i = istart, istop, istep
    if (istrue(flg(i), char)) then
      iepo = i
      return
    endif
  enddo

  return
end
