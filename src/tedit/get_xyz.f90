!
!! get_xyz.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang, Jihang Lin
!! 
!!
!
subroutine get_xyz(use_brdeph, tstart, sstart, session_length, nepo, interval, x1, x2, x3, dwnd)
  implicit none
  include '../header/const.h'
  include '../header/station.h'

! parameter
  logical*1     use_brdeph
  integer*4     tstart(5), nepo
  real*8        sstart, interval, session_length
  real*8        x1(nepo), x2(nepo), x3(nepo), dwnd
! local
  character*256 arg(50, 20), arg0, chr*1
  integer*4     i0, i, j, ipar, nval(20), npar, narg, ioerr
  integer*4     jd, jd0, jd1, iepo
  real*8        sod, sod0, sod1
  real*8        x1_tmp, x2_tmp, x3_tmp
  logical*1     display_help
  character*1   mode
  type(station) SITE
! function called
  character*1   lower_string
  integer*4     get_valid_unit
  integer*4     modified_julday
  real*8        timdif

  x1_tmp = 0.d0
  x2_tmp = 0.d0
  x3_tmp = 0.d0

  x1 = 0.d0
  x2 = 0.d0
  x3 = 0.d0

! arguement
  mode = ' '
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

!! read arguments
  ipar = 1
  do while (ipar .le. npar .and. .not. display_help)
    if (arg(ipar, 1) (1:4) .eq. '-xyz') then
      display_help = nval(ipar) .ne. 4 .and. nval(ipar) .ne. 2
      if (.not. display_help) then
        if (nval(ipar) .eq. 4) then
          mode = 'S'
          do i = 2, 4
            i0 = index(arg(ipar, i), ',')
            if (i0 .gt. 0) arg(ipar, i)(i0:i0) = '.'
          end do
          read (arg(ipar, 2), *, iostat=ioerr) x1_tmp
          read (arg(ipar, 3), *, iostat=ioerr) x2_tmp
          read (arg(ipar, 4), *, iostat=ioerr) x3_tmp
          x1 = x1_tmp
          x2 = x2_tmp
          x3 = x3_tmp
          if (use_brdeph) then
            if (x1_tmp*x2_tmp*x3_tmp .eq. 0.d0) then
              write(*,'(a)') '***ERROR(get_control_parameter): no position when using brdeph '
              call exit(1)
            end if
          end if
          display_help = ioerr .ne. 0
        else if (nval(ipar) .eq. 2) then
          mode = 'K'
          jd0 = modified_julday(tstart(3), tstart(2), tstart(1))
          sod0 = tstart(4)*3600.d0 + tstart(5)*60.d0 + sstart
          call timinc(jd0, sod0, session_length, jd1, sod1)
          jd = jd0
          sod = sod0
          iepo = 0
          SITE%kinfil = arg(ipar, 2)
          do while (timdif(jd, sod, jd1, sod1) .lt. MAXWND)
            iepo = iepo + 1
            call read_kinpos(SITE, jd+sod/864.d2, jd+sod/864.d2, x1(iepo), x2(iepo), x3(iepo), dwnd)
            call timinc(jd, sod, interval, jd, sod)
          end do
        end if
      end if
    end if
    ipar = ipar + 1
  end do

100 continue
  return
end subroutine
