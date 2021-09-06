!
!! tedit.f90
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
!! Turboedit
!
program tedit
  implicit none
  include '../header/const.h'
  include '../header/brdeph.h'
  include 'data_flag.h'

  integer*4 nsat, jd0, nobs(MAXSAT)
  integer*4, pointer :: flagall(:, :)
  real*8, pointer :: ti(:), ts(:)
  real*8 bias(MAXSAT, 44)
  real*8, pointer :: obs(:, :, :)
  integer*4 nepo, ierr

! broadcast ephemeris
  integer*4 neph
  type(brdeph) ephem(MAXEPH)

  integer*4 lfnsd
  logical*1 again, use_brdeph, check_pc, check_lc, turbo_edit, &
    debug_sd, debug_tb, keep_end
! control parameters
  character*256 flnrnx, flneph, flnrhd, string*20, stanam*4
  integer*4 tstart(5), length_gap, length_short
  real*8 sstart, cutoff_elevation, interval, pclimit, lclimit, lglimit, lgrmslimit, &
    max_mean_namb, min_percent, min_mean_nprn, session_length
! local
  integer*4 isat, iepo, nused
  real*8, pointer :: x(:), y(:), z(:)
  real*8 t_first_in_rinex, t_last_in_rinex
  real*4, pointer :: v(:, :)
! for one or single difference
  real*8, pointer :: pg(:), sigpg(:), respg(:)
  integer*4 ieph, glschn(30), lfnchn
  character*3 glssvn(30)
!! function used
  logical*1 istrue
  integer*4 get_valid_unit, modified_julday
  
!
!! get input from command line
  call get_control_parameter(flnrnx, flneph, flnrhd, check_lc, &
                             turbo_edit, use_brdeph, check_pc, keep_end, &
                             tstart, sstart, session_length, length_gap, length_short, cutoff_elevation, &
                             max_mean_namb, min_percent, min_mean_nprn, interval, lclimit, pclimit, &
                             lglimit, lgrmslimit, stanam)
  debug_tb = .true.
  debug_sd = .false.

  nepo = nint(session_length/interval)+1
  allocate (x(1:nepo),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): x allocation ', nepo
    call exit(1)
  endif
  allocate (y(1:nepo),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): y allocation ', nepo
    call exit(1)
  endif
  allocate (z(1:nepo),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): z allocation ', nepo
    call exit(1)
  endif

  call get_xyz(use_brdeph, tstart, sstart, session_length, nepo, interval, x, y, z)
  allocate (ti(1:nepo),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): ti allocation ', nepo
    call exit(1)
  endif
  allocate (ts(1:nepo),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): ts allocation ', nepo
    call exit(1)
  endif
  allocate (flagall(1:nepo,1:MAXSAT),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): flagall allocation ', nepo
    call exit(1)
  endif
  allocate (obs(1:nepo,1:MAXSAT,1:6),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): obs allocation ', nepo
    call exit(1)
  endif
  allocate (v(1:nepo,1:MAXSAT),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): v allocation ', nepo
    call exit(1)
  endif
  allocate (pg(1:nepo),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): pg allocation ', nepo
    call exit(1)
  endif
  allocate (sigpg(1:nepo),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): sigpg allocation ', nepo
    call exit(1)
  endif
  allocate (respg(1:nepo),stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(tedit): respg allocation ', nepo
    call exit(1)
  endif
! read broadcast ephemeris
  if (use_brdeph) call rdrnxn(flneph, 0.d0, 0.d0, neph, ephem)
  glssvn = ''
  glschn = 999
  do ieph = 1, neph
    if (ephem(ieph)%svn(1:1) .eq. 'R') then
      read(ephem(ieph)%svn(2:3),'(i2)') isat
      glssvn(isat) = ephem(ieph)%svn
      glschn(isat) = ephem(ieph)%fn
    endif
  enddo
  lfnchn = get_valid_unit(10)
  open(lfnchn, file="glonass_chn")
  do isat = 1, 30
    if (glschn(isat) .eq. 999) cycle
    if (glschn(isat) .gt. 50 .or. glschn(isat) .lt. -50) glschn(isat) = 999
    write(lfnchn,'(a3,1x,i10)') glssvn(isat),glschn(isat)
  enddo
  close(lfnchn)

! read rinex observation file
  bias = 0.d0
  call read_rinex_file(flnrnx, tstart, sstart, interval, &
                       check_pc, pclimit, cutoff_elevation, use_brdeph, neph, &
                       ephem, stanam, x, y, z, t_first_in_rinex, t_last_in_rinex, v, &
                       nepo, nsat, jd0, nobs, flagall, ti, ts, obs, bias)

! remove short piece and flag gap
  do isat = 1, nsat
    again = .true.
    do while (again)
! ********************************************************************** !
!               remove short piece and mark large gap                    !
! ********************************************************************** !
      call remove_short(keep_end, nepo, ti, flagall(1, isat), length_short, length_gap, interval, flag_shrt, again)
    enddo
  enddo

! SD LC checking when necessary
  lfnsd = -1
  if (debug_sd) lfnsd = 0
  if (check_lc) then
! ********************************************************************** !
!                  check single-difference LC                            !
! ********************************************************************** !
    call check_sd(lfnsd, neph, ephem, 2, 3, 20, 8, x, y, z, interval, lclimit, &
                  nepo, nsat, jd0, nobs, flagall, ti, ts, obs)
  endif
!
! widelane checking etc.
  do isat = 1, nsat
    nused = 0
    do iepo = 1, nepo
      obs(iepo, isat, 4) = 0.d0
      obs(iepo, isat, 5) = 0.d0
      if (istrue(flagall(iepo, isat), 'ok')) then
        pg(iepo) = obs(iepo, isat, 1)       ! geometry-free
        nused = nused + 1
      endif
    enddo
    if (nused .ne. 0 .and. turbo_edit) then
!! ****************************************************************** !!
!!                       Used MW obs to check epochs
!!           Nov. 1, 2007. "limit" should not be too small
!! ****************************************************************** !!
      call edit_widelane(nepo, ti, obs(1, isat, 2), flagall(1, isat), 3.0d0)
      if (check_lc) then
!! ****************************************************************** !!
!!                Check ionosphere observations: LG
!!           Nov. 1, 2007. "limit" should not be too small
!! ****************************************************************** !!
        call check_ionosphere(nepo, flagall(1, isat), ti, pg, respg, sigpg, interval, lglimit, lgrmslimit)
        call lc_help(nepo, flagall(1, isat))
      endif
    endif
    if (nused .ne. 0) then
      again = .true.
      do while (again)
        call remove_short(keep_end, nepo, ti, flagall(1, isat), length_short, length_gap, interval, flag_shrt, again)
      enddo
      if (debug_tb) then
        do iepo = 1, nepo
          if (.not. istrue(flagall(iepo, isat), 'nodata')) then
            string = 'ok'
            if (istrue(flagall(iepo, isat), 'no4')) then
              string = 'No 4'
            else if (istrue(flagall(iepo, isat), 'lowele')) then
              string = 'Low Elevation'
            else if (istrue(flagall(iepo, isat), 'shrt')) then
              string = 'Short Piece'
            else if (istrue(flagall(iepo, isat), 'lwbad')) then
              string = 'Bad Widelane'
            else if (istrue(flagall(iepo, isat), 'lgbad')) then
              string = 'Bad Ionosphere'
            else if (istrue(flagall(iepo, isat), 'lccheck')) then
              string = 'Can not check LC'
            else if (istrue(flagall(iepo, isat), 'pcbad')) then
              string = 'Bad Range   '
            else if (istrue(flagall(iepo, isat), 'pc1ms')) then
              string = 'Bad Range 1ms'
            else if (istrue(flagall(iepo, isat), 'amb')) then
              string = 'AMB'
            endif
          endif
        enddo
      endif
    endif
  enddo

! write rhd file for further processing
  if (flnrhd(1:1) .ne. ' ') then
    call write_diag_rpt(flnrhd, nepo, nsat, stanam, jd0, ts, flagall, interval)
  endif

  deallocate (ti)
  deallocate (ts)
  deallocate (flagall)
  deallocate (obs)
  deallocate (v)
  deallocate (x)
  deallocate (y)
  deallocate (z)
  deallocate (pg)
  deallocate (sigpg)
  deallocate (respg)

end
