!
!! check_ionosphere.f90
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
!! Check ionosphere observations: LG
!!
!! input : pg = obs(iepo,iprn,1)    geometry-free
!
subroutine check_ionosphere(nepo, flagall, ti, pg, res, sig, interval, lglimit, lgrmslimit)
  implicit none
  integer*4 nepo, flagall(1:*)
  real*8 ti(1:*), pg(1:*), lglimit, lgrmslimit, res(1:*), sig(1:*)
  logical*1 use_brdeph, strict, set

  logical*1 istrue
  integer*4 set_flag
  real*8 v(nepo), rms, interval, a0
  integer*4 iepo, flg(nepo), nwithin, npt0, npt1, npt2, ndgr, niter, jepo, ipt(nepo), i0, i1, kobs, ierr

  npt0 = 3
  nwithin = 8
  ndgr = 2
  niter = 2

  do iepo = 1, nepo
    res(iepo) = 0.d0
    sig(iepo) = 999.d0
    flg(iepo) = 2
    if (istrue(flagall(iepo), 'ok')) then
      flg(iepo) = 1
      if (.not. istrue(flagall(iepo), 'lgjump')) flg(iepo) = 0
    endif
  enddo

  do iepo = 1, nepo
    if (flg(iepo) .ne. 2) then
      set = .false.
! identify usable pg before iepo
! jepo = last epoch
      jepo = iepo - 1
      npt1 = 0
      i0 = iepo
      do while (npt1 .lt. npt0 .and. iepo - jepo .le. nwithin .and. jepo .ge. 1)
        if (flg(jepo) .lt. 2) then
          npt1 = npt1 + 1
          i0 = jepo
        endif
        jepo = jepo - 1
      enddo

      if (npt1 .eq. 0) then
        set = .true.
        goto 100
      endif
! identify usable pg after iepo
      jepo = iepo
      npt2 = 0
      i1 = iepo
      do while (npt2 .lt. npt0 .and. jepo - iepo .lt. nwithin .and. jepo .le. nepo)
        if (flg(jepo) .lt. 2) then
          npt2 = npt2 + 1
          i1 = jepo
        endif
        jepo = jepo + 1
      enddo
! polydf to identify carrier phase jumps
      kobs = i1 - i0 + 1
!! purpose  : Check jump using polynomial fitting
      call check_for_jump('LGC', -1, kobs, ti(i0), pg(i0), flg(i0), ndgr, niter, &
                          0.3d0, 0.3d0, a0, rms, v(i0), ipt(i0), ierr, interval)
      if (ierr .eq. 0) then
        res(iepo) = v(iepo)
        sig(iepo) = rms
        if (lgrmslimit .eq. 0.d0) then
          if (dabs(v(iepo)) .gt. lglimit) set = .true.
        else
          if (rms .gt. lgrmslimit .or. dabs(v(iepo)) .gt. max(3*rms, lglimit)) set = .true.
        endif
      else
        set = .true. ! if polydf failed, then ...
      endif
100   continue
      if (set) flagall(iepo) = set_flag(flagall(iepo), 'lgjump')
    endif
  enddo

  return
end
