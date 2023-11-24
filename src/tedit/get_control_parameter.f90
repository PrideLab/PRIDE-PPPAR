!
!! get_control_parameter.f90
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
!
subroutine get_control_parameter(flnrnx, flneph, flnrhd, &
                                 check_lc, turbo_edit, lm_edit, use_brdeph, check_pc, keep_end, trunc_dbd, &
                                 tstart, sstart, session_length, length_gap, length_short, &
                                 cutoff_elevation, max_mean_namb, min_percent, min_mean_nprn, &
                                 interval, lclimit, pclimit, lglimit, lgrmslimit, &
                                 stanam)
  implicit none
  include '../header/const.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  character*256 flnrnx, flneph, flnrhd
  logical*1     check_lc, turbo_edit, lm_edit, use_brdeph, check_pc, keep_end
  character*1   trunc_dbd
  integer*4     tstart(5)
  real*8        sstart, session_length
  integer*4     length_gap, length_short
  real*8        cutoff_elevation
  real*8        max_mean_namb, min_percent, min_mean_nprn
  real*8        interval, lclimit, pclimit, lglimit, lgrmslimit
  character*4   stanam
! local
  character*256 arg(50, 20), arg0, message, chr*1
  integer*4     i0, i, j, k, ipar, nval(20), npar, narg, ioerr
  logical*1     lexist, display_help
! function called
  character*1   lower_string

  flneph = ' '
  flnrhd = ' '

!
!! set default parameters
  check_lc = .true.
  check_pc = .false.
  turbo_edit = .true.
  lm_edit = .false.
  use_brdeph = .true.
  keep_end = .true.
  trunc_dbd = 'n'
  tstart = 0
  sstart = 0.d0
  length_gap = 600
  length_short = 1800
  cutoff_elevation = 0.d0
  max_mean_namb = 3.d0
  min_percent = 80.d0
  min_mean_nprn = 5.d0
  interval = 0.d0
  lclimit = 1.d0
  pclimit = 250.d0
  lglimit = -10.d0
  lgrmslimit = -10.d0
!
!! parse arguements
  message = ' '
  display_help = .false.
  narg = iargc()
  if (narg .eq. 0) display_help = .true.
  npar = 0
! first argument is always command name. But the first one
! is start at 0 or 1 is compiler dependant.
  call getarg(0, arg0)
  if (arg0(1:1) .eq. ' ') then
    i0 = 2
    narg = narg + 1
  else
    i0 = 1
  end if
  i = i0
  do while (i .le. narg .and. .not. display_help)
    call getarg(i, arg0)
    if (arg0(1:5) .eq. '-help') then
      display_help = .true.
    else if (i .eq. i0) then
      flnrnx = arg0
      inquire (file=arg0, exist=lexist)
      display_help = .not. lexist
      message = 'Error: rinex_obs_file '
    else
      chr = lower_string(arg0(2:2))
      if (arg0(1:1) .eq. '-' .and. ichar(chr) .ge. ichar('a') .and. &
          ichar(chr) .le. ichar('z')) then
        npar = npar + 1
        nval(npar) = 1
      else
        nval(npar) = nval(npar) + 1
      end if
      arg(npar, nval(npar)) = arg0
    end if
    i = i + 1
  end do

!
!! read arguments
  ipar = 1
  do while (ipar .le. npar .and. .not. display_help)
    message = 'Error: '//arg(ipar, 1)
    if (arg(ipar, 1)(1:5) .eq. '-rnxn') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        flneph = arg(ipar, 2)
        inquire (file=flneph, exist=lexist)
        display_help = .not. lexist
      end if
    else if (arg(ipar, 1)(1:4) .eq. '-rhd') then
      display_help = nval(ipar) .lt. 1
      if (.not. display_help) flnrhd = arg(ipar, nval(ipar))
      if (.not. display_help) stanam = flnrhd(13:16)
    else if (arg(ipar, 1)(1:5) .eq. '-time') then
      display_help = nval(ipar) .ne. 7
      do j = 1, nval(ipar) - 2
        read (arg(ipar, j + 1), *, iostat=ioerr) tstart(j)
        if (ioerr .ne. 0) display_help = .true.
      end do
      read (arg(ipar, nval(ipar)), *, iostat=ioerr) sstart
      if (ioerr .ne. 0) display_help = .true.
    else if (arg(ipar, 1)(1:4) .eq. '-int') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) interval
        if (ioerr .ne. 0) display_help = .true.
      end if
    else if (arg(ipar, 1)(1:4) .eq. '-len') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) session_length
        if (ioerr .ne. 0) display_help = .true.
      end if
    else if (arg(ipar, 1)(1:6) .eq. '-short') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) length_short
        if (ioerr .ne. 0) display_help = .true.
      end if
    else if (arg(ipar, 1)(1:10) .eq. '-trunc_dbd') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) trunc_dbd
        if (ioerr .ne. 0) display_help = .true.
      end if
    else if (arg(ipar, 1)(1:9) .eq. '-lc_check') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        if (arg(ipar, 2)(1:1) .eq. 'y') then
          check_lc = .true.
        else if (arg(ipar, 2)(1:1) .eq. 'n') then
          check_lc = .false.
        else if (arg(ipar, 2)(1:1) .eq. 'o') then
          check_lc = .true.
          turbo_edit = .false.
        else if (arg(ipar, 2)(1:1) .eq. 'l') then
          check_lc = .false.
          lm_edit = .true.
        else
          display_help = .true.
        end if
      end if
    else if (arg(ipar, 1)(1:9) .eq. '-pc_check') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        check_pc = .true.
        read (arg(ipar, 2), *, iostat=ioerr) pclimit
        display_help = ioerr .ne. 0
      end if
    else if (arg(ipar, 1)(1:4) .eq. '-xyz') then
      display_help = nval(ipar) .ne. 4 .and. nval(ipar) .ne. 2
    else if (arg(ipar, 1)(1:5) .eq. '-elev') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) cutoff_elevation
        display_help = ioerr .ne. 0
      end if
    else if (arg(ipar, 1)(1:5) .eq. '-freq') then
      display_help = (nval(ipar) .lt. 2 .or. nval(ipar) .gt. MAXSYS + 1)
      do j = 2, nval(ipar)
        if (len_trim(arg(ipar, j)) .ne. 3) then
          display_help = .true.
        else
          i0 = index(GNSS_PRIO, arg(ipar, j)(1:1))
          if (i0 .le. 0) display_help = .true.
          do k = 1, 2
            read (arg(ipar, j)(k+1:k+1), '(i1)', iostat=ioerr) idxfrq(i0, k)
            if (ioerr .ne. 0) display_help = .true.
          end do
        end if
      end do
    else
      display_help = .true.
    end if
    ipar = ipar + 1
  end do
!
!! set default frequency
  i0 = index(GNSS_PRIO, 'G')
  if (any(idxfrq(i0, :) .eq. 0)) idxfrq(i0, :) = [1, 2]
  i0 = index(GNSS_PRIO, 'R')
  if (any(idxfrq(i0, :) .eq. 0)) idxfrq(i0, :) = [1, 2]
  i0 = index(GNSS_PRIO, 'E')
  if (any(idxfrq(i0, :) .eq. 0)) idxfrq(i0, :) = [1, 5]
  i0 = index(GNSS_PRIO, 'C')
  if (any(idxfrq(i0, :) .eq. 0)) idxfrq(i0, :) = [2, 6]
  i0 = index(GNSS_PRIO, 'J')
  if (any(idxfrq(i0, :) .eq. 0)) idxfrq(i0, :) = [1, 2]

!
!! check consistence
  if (.not. display_help) then
    if (flneph(1:1) .eq. ' ') then
      if (check_lc) then
        message = '***ERROR(get_control_parameter): LC check need broadcast ephem. file. &
&                  Use -rnxn or -no_lc_check'
        display_help = .true.
      else if (cutoff_elevation .ne. 0) then
        message = '***ERROR(get_control_parameter): elevation edit need broadcast ephem. &
&                  file. Use -rnxn or set elev. to zero'
        display_help = .true.
      else
        use_brdeph = .false.
      end if
    end if
    i = len_trim(flnrnx)
    if (flnrhd(1:1) .eq. '-') flnrhd = flnrnx(1:i)//'.rhd'
    if (lglimit .eq. -10.d0) then
      if (check_lc) then
        lglimit = 1.d0
        lgrmslimit = 1.8d0
      else
        lglimit = 1.d0
        lgrmslimit = 0.3d0
      end if
    end if
  end if

!
!! display help
100 continue
  if (display_help) then
    call getarg(0, arg0)
    if (arg0(1:1) .eq. ' ') call getarg(1, arg0)
    write (*, '(/2a/)') trim(arg0), ' Version 1.0, Wuhan University, Aug. 2007 '
    write (*, '(a,2a)') 'usage: ', trim(arg0), ' rinex_obs_file [option] / '
    write (*, '(2x,a)') '-rnxn filename '
    write (*, '(2x,a)') ' broadcast ephem. file. If -check_lc is active or -elev is on, this file is required'
    write (*, '(2x,a)') '-rhd [filename]'
    write (*, '(2x,a)') ' output rhd file. Default is no. Default name is {rinx_obs_file}.rhd'
    write (*, '(2x,a)') '-time year month day hour minut second '
    write (*, '(2x,a)') ' start time for data editing. Default is the start time in rinex file'
    write (*, '(2x,a)') '-len length_in_second'
    write (*, '(2x,a)') ' length of data to be edited. Default is all data in rinex file'
    write (*, '(2x,a)') '-int interval_in_second'
    write (*, '(2x,a)') ' sampling interval for data editing. Default is 30 seconds'
    write (*, '(2x,a)') '-short length_in_second'
    write (*, '(2x,a)') ' data piece shorter than this value will be removed. Default is 600s'
    write (*, '(2x,a)') '-trunc_dbd yes/no/cont'
    write (*, '(2x,a)') ' truncated all ambiguities at day-boundary' 
    write (*, '(2x,a)') '    yes = truncate'
    write (*, '(2x,a)') '     no = no special action'
    write (*, '(2x,a)') '   cont = keep continuous'
    write (*, '(2x,a)') '-elev cutoff_elevation'
    write (*, '(2x,a)') ' cutoff_elevation in degree. Default is to use all data.'
    write (*, '(2x,a)') '-freq GNSS with frequency number'
    write (*, '(2x,a)') ' dual-frequency combination for G E C J. Default is G12 E15 C26 J12'
    write (*, '(2x,a)') '-lc_check yes/no/only'
    write (*, '(2x,a)') '    yes = check LC and edit WL and IONO and try to connect WL and ION'
    write (*, '(2x,a)') '     no = edit WL and IONO and try to connect WL and ION observation'
    write (*, '(2x,a)') '   only = check LC only'
    if (message .ne. ' ') write (*, '(/a120)') message
    call exit(1)
  end if
  return
end subroutine
