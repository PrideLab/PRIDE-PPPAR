!
!! find_ambd.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : find death epoch of an ambiguity
!! parameter:
!!    input : nepo -- total number of flag
!!            flag -- flag array
!!            iepo -- starting epoch of amb
!!    output: find_flag -- death epoch of amb
!
integer*4 function find_ambd(nepo, flag, iepo)
  implicit none

  integer*4 nepo, flag(1:*), iepo
!
!! local
  integer*4 i

  if (iepo .eq. nepo) then
    find_ambd = iepo
    return
  endif
  call find_flag(iepo + 1, nepo, flag, 'amb', i)
  if (i .lt. 0) i = nepo + 1
  call find_flag(i - 1, iepo, flag, 'ok', find_ambd)

  return
end
