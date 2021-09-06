!
!! remove_short.f90
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
!!
!!    remove short piece and mark large gap
!
subroutine remove_short(keep_end, nepo, ti, flg, len_short, len_gap, &
                        interval, flag_shrt, removed)
  implicit none

  integer*4 nepo, flg(1:*), len_short, flag_shrt, len_gap
  real*8 ti(1:*), interval
  logical*1 removed, keep_end

! local
  logical*1 found, remove_it, istrue
  integer*4 i, j, k, iepo, ngood, set_flag, ilast

! mark large gap and first epoch (for check_sd)
  ilast = 0
  do iepo = 1, nepo
    if (istrue(flg(iepo), 'ok')) then
      if ((ilast .ne. 0 .and. ti(iepo) - ti(ilast) .gt. len_gap) .or. ilast .eq. 0) then
        flg(iepo) = set_flag(flg(iepo), 'gap')
      endif
      ilast = iepo
    endif
  enddo

! last amb, remove it
  if (ilast.ne.0 .and. istrue(flg(ilast), 'AMB')) then
    flg(ilast) = set_flag(flg(ilast), 'shrt')
  endif

! remove short piece
  j = 1
  removed = .false.
  found = .true.
  do while (found)
!!
!! find arc corresponding to an ambiguity
    call find_flag(j, nepo, flg, 'AMB', i)
    if (i .gt. 0) then
      call find_flag(i + 1, nepo, flg, 'AMB', j)
      if (j .lt. 0) then
        found = .false.
        j = nepo + 1
      endif
      call find_flag(j - 1, i, flg, 'OK', k)
      remove_it = .false.
      ngood = 1
      do iepo = i + 1, k
        if (istrue(flg(iepo), 'ok')) ngood = ngood + 1
      enddo
!!
!! check arc length
!! keep_end : keep short ending arc
      if ((ti(k) - ti(i + 1) .lt. len_short .or. ngood*interval .le. len_short/2) &
          .and. ((k .lt. nepo - 1 .and. i .gt. 2) .or. .not. keep_end)) then
        remove_it = .true.
      endif
      if (remove_it) then
        do iepo = i, k
          if (istrue(flg(iepo), 'ok')) then
            flg(iepo) = set_flag(flg(iepo), 'shrt')
          endif
        enddo
        removed = .true.
      endif
      j = k + 1
    else
      found = .false.
    endif
  enddo

  return
end
