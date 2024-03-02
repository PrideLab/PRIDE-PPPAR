!
!! tedit.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin, Jing Zeng
!!
!!
!!
!! Turboedit
!
program tedit
  implicit none
  include '../header/const.h'
  include '../header/absbia.h'
  include '../header/brdeph.h'
  include 'data_flag.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
!
  integer*4     jd0
  integer*4     nsat, msat
  integer*4     nobs(MAXSAT)
  integer*4, pointer :: flag(:, :)
  real*8, pointer :: ti(:), ts(:)
  real*8, pointer :: x(:), y(:), z(:)
  integer*4     kday, nday
  integer*4     kepo, jepo
  integer*4     mepo, nepo
  integer*4     ierr
  type(absbia)  bias(MAXSAT, MAXTYP)
  character*3   prn0(MAXSAT)
!! broadcast ephemerides
  integer*4     neph
  type(brdeph), pointer :: ephem(:)
!! logical parameters
  logical*1     again
  logical*1     use_brdeph
  logical*1     check_pc, check_lc
  logical*1     turbo_edit, lm_edit, ltighter
  logical*1     debug_sd, debug_tb
  logical*1     keep_end
  logical*1     lexist
! control parameters
  character*1   trunc_dbd
  character*256 flnrnx, flneph, flnrhd
  character*20  string
  character*4   stanam
  integer*4     tstart(5), length_gap, length_short
  real*8        sstart, interval, cutoff_elevation
  real*8        pclimit, lclimit, lglimit, lgrmslimit
  real*8        max_mean_namb, min_percent, min_mean_nprn, session_length
! local
  integer*4     lfnsd
  integer*4     nused
  integer*4     isat, iepo
  integer*4, pointer :: tmpflg(:, :)
  real*8, pointer :: obs(:, :, :)
  real*8, pointer :: tti(:), tts(:)
  real*8, pointer :: xt(:), yt(:), zt(:)
  real*4, pointer :: vel(:, :)
  real*8        t_first_in_rinex, t_last_in_rinex
  real*8        fjd0, fjd1, overlap_sec, tmp_session
! for one or single difference
  real*8, pointer :: valpg(:), sigpg(:), respg(:)
! function used
  logical*1     istrue
  integer*4     get_valid_unit
  integer*4     modified_julday

  call prn_matbld(prn0)

!
!! initialize
  idxfrq = 0
  debug_tb = .true.
  debug_sd = .false.
!
!! get input from command line
  call get_control_parameter(flnrnx, flneph, flnrhd, &
                             check_lc, turbo_edit, lm_edit, ltighter, use_brdeph, check_pc, keep_end, trunc_dbd, &
                             tstart, sstart, session_length, length_gap, length_short, &
                             cutoff_elevation, max_mean_namb, min_percent, min_mean_nprn, &
                             interval, lclimit, pclimit, lglimit, lgrmslimit, &
                             stanam)
!
!! calculate time span
  fjd0 = modified_julday(tstart(3), tstart(2), tstart(1)) + tstart(4)/24.d0 + tstart(5)/144.d1 + sstart/864.d2
  fjd1 = fjd0 + session_length/864.d2

  nday = ceiling(fjd1 + 1d-9) - floor(fjd0)
  nepo = nint(session_length/interval) + 1
!
!! allocate global arrays
  allocate (ti(1:nepo), stat=ierr)
  allocate (ts(1:nepo), stat=ierr)
  allocate (flag(1:nepo, 1:MAXSAT), stat=ierr)
!
!! determine the first session
  if (nday .le. 3) then
  !
  !! while nday <= 3, all data will be processed in a single session
  !!   (e.g. 59488.5 - 59500.5)
  !
    kday = nday
    overlap_sec = 0.d0
    tmp_session = session_length
  else
  !
  !! the first session is ended at 24 h in the first day
  !! an overlap of 0.5 h is appended to avoid saving short piece while save_end=true
  !
    kday = 1
    overlap_sec = 36.d2
    tmp_session = 864.d2*(ceiling(fjd0 + 1.d-9) - fjd0) + overlap_sec
  end if

  ! round to the nearest millisecond
  tmp_session = nint(tmp_session*1.d3)/1.d3

!
!! divide a long-period processing into several sessions
  kepo = 1
  jepo = 1
  do while (kday .le. nday)

    ! number of epoch for this session
    mepo = nint(tmp_session/interval) + 1

    ! allocate temporary arrays
    allocate (xt(1:mepo),    stat=ierr)
    allocate (yt(1:mepo),    stat=ierr)
    allocate (zt(1:mepo),    stat=ierr)
    allocate (tmpflg(1:mepo, 1:MAXSAT),   stat=ierr)
    allocate (obs(1:mepo, 1:MAXSAT, 1:7), stat=ierr)
    allocate (vel(1:mepo, 1:MAXSAT),      stat=ierr)
    allocate (tti(1:mepo),   stat=ierr)
    allocate (tts(1:mepo),   stat=ierr)
    allocate (valpg(1:mepo), stat=ierr)
    allocate (sigpg(1:mepo), stat=ierr)
    allocate (respg(1:mepo), stat=ierr)

    ! read a priori xyz
    call get_xyz(use_brdeph, tstart, sstart, tmp_session, mepo, interval, xt, yt, zt)

    ! read broadcast ephemerids
    if (use_brdeph) then
      neph = 3*(MAXSAT*24 + MAXSAT_E*120)
      allocate (ephem(neph))
      call rdrnxn(flneph, fjd0 + sstart/864.d2, fjd0 + (sstart + tmp_session)/864.d2, neph, ephem)
    end if

    ! read rinex observation file
    call read_rinex_file(flnrnx, tstart, sstart, interval, &
                         check_pc, pclimit, &
                         cutoff_elevation, use_brdeph, neph, ephem, lm_edit, ltighter, &
                         stanam, xt, yt, zt, t_first_in_rinex, t_last_in_rinex, vel, &
                         mepo, nsat, jd0, nobs, tmpflg, tti, tts, obs, bias)

    ! remove short piece and flag gap
    do isat = 1, nsat
      again = .true.
      do while (again)
        ! ********************************************************************** !
        !               remove short piece and mark large gap                    !
        ! ********************************************************************** !
        if (prn0(isat)(1:1) .eq. 'R') then
          call remove_short(keep_end,       'y', mepo, tti, tmpflg(1, isat), length_short, length_gap, interval, flag_shrt, again)
        else
          call remove_short(keep_end, trunc_dbd, mepo, tti, tmpflg(1, isat), length_short, length_gap, interval, flag_shrt, again)
        end if
      end do
    end do

    ! SD LC checking when necessary
    lfnsd = -1
    if (debug_sd) lfnsd = 0
    if (check_lc) then
      ! ********************************************************************** !
      !                  check single-difference LC                            !
      ! ********************************************************************** !
      call check_sd(lfnsd, neph, ephem, 2, 3, 20, 8, xt, yt, zt, interval, lclimit, &
                    mepo, nsat, jd0, nobs, tmpflg, tti, tts, obs)
    end if
    !
    ! widelane checking etc.
    do isat = 1, nsat
      nused = 0
      do iepo = 1, mepo
        obs(iepo, isat, 4) = 0.d0
        obs(iepo, isat, 5) = 0.d0
        if (istrue(tmpflg(iepo, isat), 'ok')) then
          ! geometry-free
          valpg(iepo) = obs(iepo, isat, 1)
          nused = nused + 1
        end if
      end do
      if (nused .ne. 0 .and. turbo_edit) then
  !! ****************************************************************** !!
  !!                       Used MW obs to check epochs
  !!           Nov. 1, 2007. "limit" should not be too small
  !!           Dec. 6, 2022. improved with FBMWA and STPIR algorithm
  !! ****************************************************************** !!
        if (.not. lm_edit) then
          call edit_widelane(mepo, tti, obs(1, isat, 2), tmpflg(1, isat), 3.0d0)
        else
          call lfbmwa(mepo, tti, obs(1, isat, 2), obs(1, isat, 7), tmpflg(1, isat), 1.0d0)
          call edit_widelane(mepo, tti, obs(1, isat, 2), tmpflg(1, isat), 3.0d0)
          call mstpir(mepo, tmpflg(1, isat), tti, valpg)
        end if
        if (check_lc) then
  !! ****************************************************************** !!
  !!                Check ionosphere observations: LG
  !!           Nov. 1, 2007. "limit" should not be too small
  !! ****************************************************************** !!
          call check_ionosphere(mepo, tmpflg(1, isat), tti, valpg, respg, sigpg, interval, lglimit, lgrmslimit)
          call lc_help(mepo, tmpflg(1, isat))
        end if
      end if
      if (nused .ne. 0) then
        again = .true.
        do while (again)
          if (prn0(isat)(1:1) .eq. 'R') then
            call remove_short(keep_end,       'y', mepo, tti, tmpflg(1, isat), length_short, length_gap, interval, flag_shrt, again)
          else
            call remove_short(keep_end, trunc_dbd, mepo, tti, tmpflg(1, isat), length_short, length_gap, interval, flag_shrt, again)
          end if
        end do
        if (debug_tb) then
          do iepo = 1, mepo
            if (istrue(tmpflg(iepo, isat), 'nodata')) cycle
            string = 'ok'
            if (istrue(tmpflg(iepo, isat), 'no4')) then
              string = 'No 4'
            else if (istrue(tmpflg(iepo, isat), 'lowele')) then
              string = 'Low Elevation'
            else if (istrue(tmpflg(iepo, isat), 'shrt')) then
              string = 'Short Piece'
            else if (istrue(tmpflg(iepo, isat), 'lwbad')) then
              string = 'Bad Widelane'
            else if (istrue(tmpflg(iepo, isat), 'lgbad')) then
              string = 'Bad Ionosphere'
            else if (istrue(tmpflg(iepo, isat), 'lccheck')) then
              string = 'Can not check LC'
            else if (istrue(tmpflg(iepo, isat), 'pcbad')) then
              string = 'Bad Range   '
            else if (istrue(tmpflg(iepo, isat), 'pc1ms')) then
              string = 'Bad Range 1ms'
            else if (istrue(tmpflg(iepo, isat), 'amb')) then
              string = 'AMB'
            end if
          end do
        end if
      end if
    end do

    if (nday .gt. 3 .and. kday .ne. 1) then
      jepo = jepo - nint(length_short/interval)
    end if

    ! save time tags
    ti(jepo:(kepo + mepo - 1)) = tti((jepo - kepo + 1):mepo)
    ts(jepo:(kepo + mepo - 1)) = tts((jepo - kepo + 1):mepo)

    ! save flags
    flag(jepo:(kepo + mepo - 1), :) = tmpflg((jepo - kepo + 1):mepo, :)

    ! deallocate temporary arrays
    deallocate (tmpflg)
    deallocate (obs)
    deallocate (vel)
    deallocate (tti)
    deallocate (tts)
    deallocate (valpg)
    deallocate (sigpg)
    deallocate (respg)
    deallocate (ephem)

    ! accumulate sstart
    sstart = 864.d2*(kday - 1 + ceiling(fjd0 + 1.d-9) - fjd0)
    jepo = kepo + mepo 
    kepo = nint(sstart/interval) + 1

    ! come to next session
    kday = kday + 1
    if (kday .eq. nday) then
      tmp_session = 864.d2*(fjd1 - floor(fjd1))
    else
      tmp_session = 864.d2 + overlap_sec
    end if

    ! round to the nearest millisecond
    sstart = dnint(sstart*1.d3)/1.d3
    tmp_session = dnint(tmp_session*1.d3)/1.d3

  end do

!
!! write rhd file for further processing
  if (flnrhd(1:1) .ne. ' ') then
    call write_diag_rpt(flnrhd, nepo, MAXSAT, stanam, jd0, ts, session_length, flag, interval, trunc_dbd)
  end if
!
!! deallocate global arrays
  deallocate (ti)
  deallocate (ts)
  deallocate (flag)
end program
