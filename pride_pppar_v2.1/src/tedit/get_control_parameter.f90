!
!! get_control_parameter.f90
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
!
subroutine get_control_parameter(flnrnx, flneph, flnrhd, &
                                 check_lc, turbo_edit, use_brdeph, check_pc, &
                                 keep_end, tstart, sstart, session_length, length_gap, &
                                 length_short, cutoff_elevation, max_mean_namb, min_percent, &
                                 min_mean_nprn, interval, lclimit, pclimit, lglimit, lgrmslimit, &
                                 stanam)
  implicit none
  include '../header/const.h'

  logical*1 check_pc, check_lc, use_brdeph, turbo_edit, keep_end
  character*256 flnrnx, flneph, flnrhd, flnman, flnsta
  character*4   stanam
  integer*4 tstart(5), length_gap, length_short
  real*8 sstart, cutoff_elevation, interval, pclimit, lclimit, lglimit, lgrmslimit, &
    limit_man, max_mean_namb, min_percent, min_mean_nprn, session_length

! local
  character*256 arg(50, 20), arg0, message, chr*1
  integer*4 i0, i, j, ipar, nval(20), npar, narg, ioerr
  logical*1 lexist, display_help
  character*1 mode

! function called
  character*1 lower_string

  flneph = ' '
  flnrhd = ' '
  flnman = ' '
  flnsta = ' '


! set default parameters
  do i = 1, 5
    tstart(i) = 0
  enddo
  sstart = 0.d0
  session_length = 0.d0 
  interval = 0.d0
  cutoff_elevation = 0.d0

  length_gap = 600
  length_short = 1800

  check_lc = .true.
  check_pc = .false.
  turbo_edit = .true.
  use_brdeph = .true.
  keep_end = .true.
  lclimit = 1.d0
  pclimit = 250.d0

  lglimit = -10.d0
  lgrmslimit = -10.d0

  max_mean_namb = 3.d0
  min_percent = 80.d0
  min_mean_nprn = 5.d0

! arguement
  mode = ' '
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
  endif
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
      endif
      arg(npar, nval(npar)) = arg0
    endif
    i = i + 1
  enddo

!! read arguments
  ipar = 1
  do while (ipar .le. npar .and. .not. display_help)
    message = 'Error: '//arg(ipar, 1)
    if (arg(ipar, 1) (1:5) .eq. '-rnxn') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        flneph = arg(ipar, 2)
        inquire (file=flneph, exist=lexist)
        display_help = .not. lexist
      endif
    else if (arg(ipar, 1) (1:4) .eq. '-rhd') then
      display_help = nval(ipar) .lt. 1
      if (.not. display_help) flnrhd = arg(ipar, nval(ipar))
      if (.not. display_help) stanam = flnrhd(13:16)
    else if (arg(ipar, 1) (1:5) .eq. '-time') then
      display_help = nval(ipar) .ne. 7
      do j = 1, nval(ipar) - 2
        read (arg(ipar, j + 1), *, iostat=ioerr) tstart(j)
        if (ioerr .ne. 0) display_help = .true.
      enddo
      read (arg(ipar, nval(ipar)), *, iostat=ioerr) sstart
      if (ioerr .ne. 0) display_help = .true.
    else if (arg(ipar, 1) (1:4) .eq. '-int') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) interval
        if (ioerr .ne. 0) display_help = .true.
      endif
    else if (arg(ipar, 1) (1:4) .eq. '-len') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) session_length
        if (ioerr .ne. 0) display_help = .true.
      endif
    else if (arg(ipar, 1) (1:6) .eq. '-short') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) length_short
        if (ioerr .ne. 0) display_help = .true.
      endif
    else if (arg(ipar, 1) (1:9) .eq. '-lc_check') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        if (arg(ipar, 2) (1:1) .eq. 'y') then
          check_lc = .true.
        else if (arg(ipar, 2) (1:1) .eq. 'n') then
          check_lc = .false.
        else if (arg(ipar, 2) (1:1) .eq. 'o') then
          check_lc = .true.
          turbo_edit = .false.
        else
          display_help = .true.
        endif
      endif
    else if (arg(ipar, 1) (1:9) .eq. '-pc_check') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        check_pc = .true.
        read (arg(ipar, 2), *, iostat=ioerr) pclimit
        display_help = ioerr .ne. 0
      endif
    else if (arg(ipar, 1) (1:4) .eq. '-xyz') then
      display_help = nval(ipar) .ne. 4 .and. nval(ipar) .ne. 2
    else if (arg(ipar, 1) (1:5) .eq. '-elev') then
      display_help = nval(ipar) .ne. 2
      if (.not. display_help) then
        read (arg(ipar, 2), *, iostat=ioerr) cutoff_elevation
        display_help = ioerr .ne. 0
      endif
    else
      display_help = .true.
    endif
    ipar = ipar + 1
  enddo

! check consistence
  if (.not. display_help) then
    if (flneph(1:1) .eq. ' ') then
      if (check_lc) then
        message = '***ERROR(get_control_parameter): LC check need broadcast ephem. file. &
                  Use -rnxn or -no_lc_check'
        display_help = .true.
      else if (cutoff_elevation .ne. 0) then
        message = '***ERROR(get_control_parameter): elevation edit need broadcast ephem. &
                  file. Use -rnxn or set elev. to zero'
        display_help = .true.
      else
        use_brdeph = .false.
      endif
    endif
    i = len_trim(flnrnx)
    if (flnrhd(1:1) .eq. '-') flnrhd = flnrnx(1:i)//'.rhd'
    if (lglimit .eq. -10.d0) then
      if (check_lc) then
        lglimit = 1.d0
        lgrmslimit = 1.8d0
      else
        lglimit = 1.d0
        lgrmslimit = 0.3d0
      endif
    endif
  endif

100 continue
  if (display_help) then
    call getarg(0,arg0)
    if(arg0(1:1).eq.' ') call getarg(1,arg0)
    write(*,'(/)')
    write(*,'(a12,a)') arg0, ' Version 1.0, Wuhan University, Aug. 2007 '
    write(*,'(a,a12,a,/)') ' Usage: ',arg0,' rinex_obs_file [option] / '
    write(*,'(a/a)') ' -rnxn filename ','  broadcast ephem. file. If -check_lc &
                      is active or -elev is on, this file is required. '
    write(*,'(a/a)') ' -rhd [filename]','  output rhd file. Default is no. &
                     Default is no. Default name is {rinx_obs_file}.rhd '
    write(*,'(a/a)') ' -time year month day hour minut second ', &
                     '  start time for data editing. Default is the start time in rinex file '
    write(*,'(a/a)') ' -len length_in_second','  length of data to be edited. &
                     Default is all data in rinex file'
    write(*,'(a/a)') ' -int interval_in_second','  sampling interval for data editing.&
                     Default is 30 seconds'
    write(*,'(a/a)') ' -short length_in_second','  data piece shorter than this value &
                     will be removed. Default is 600s'
    write(*,'(a/a)') ' -elev cutoff_elevation','  cutoff_elevation in degree. &
                     Default is to use all data.'
    write(*,'(a/a/a/a)') ' -lc_check yes/no/only','  yes = check LC and edit WL and &
                         IONO and try to connect WL and ION.','  no = edit WL and IONO &
                         and try to connect WL and ION observation.','  only= check LC only'
    write(*,'(/a120)')  message
    call exit(1)
  endif

  return
end
