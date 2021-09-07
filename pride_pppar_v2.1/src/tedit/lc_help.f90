!
!! lc_help.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
subroutine lc_help(nepo, flagall)
  integer*4 nepo, flagall(1:*)

  integer*4 set_flag
  logical*1 istrue

  integer*4 iepo

  do iepo = 1, nepo
    if (istrue(flagall(iepo), 'ok')) then
      if (.not. istrue(flagall(iepo), 'bigsd') .and. &
          .not. istrue(flagall(iepo), 'gap')) then
        if (.not. istrue(flagall(iepo), 'lgjump') .or. &
            .not. istrue(flagall(iepo), 'lwjump')) then
          flagall(iepo) = set_flag(flagall(iepo), 'lggood')
          flagall(iepo) = set_flag(flagall(iepo), 'lwgood')
        endif
      endif
    endif
  enddo

  return
end
