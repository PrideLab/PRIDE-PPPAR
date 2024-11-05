!
!! rdrnxoh.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Jihang Lin
!!
!!
!!
!! purpose  : read header of a RINEX o-file
!!
!! parameter: lfn -- file unit
!!            HD  -- rinex head structure, see `rinex_observation.h`
!!            ierr -- error code
!!
!!
subroutine rdrnxoh(lfn, HD, ierr)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'

  integer*4     lfn, ierr
  type(rnxhdr)  HD
! local
  integer*4     i, ioerr, ii, i0
  real*8        ver
  character*80  msg, line, keyword*20
  character*1   type_sys

!
!! initialize
  HD%t0  = 0.d0
  HD%t0s = 0.d0
  HD%t1  = 0.d0
  HD%t1s = 0.d0
  HD%intv = 0.d0
  HD%x = 0.d0
  HD%y = 0.d0
  HD%z = 0.d0
  HD%h = 0.d0
  HD%e = 0.d0
  HD%n = 0.d0

  ierr = 0
  do while (.true.)
    msg = '   '
    read (lfn, '(a)', end=200) line
    keyword = line(61:80)
!
!! end of header
    if ((index(keyword, 'END OF HEADER') .ne. 0 .and. HD%ver .ge. 200 .and.  HD%ver .lt. 300) &
        .or. (len_trim(keyword) .eq. 0 .and. HD%ver .lt. 200) &
        .or. (index(keyword, 'END OF HEADER') .ne. 0 .and. HD%ver .ge. 300 .and.  HD%ver .lt. 400) &
        .or. (index(keyword, 'END OF HEADER') .ne. 0 .and. HD%ver .ge. 400 .and.  HD%ver .lt. 500)) return
!
!! RINEX version
    if (index(keyword, 'RINEX VERSION') .ne. 0) then
      read (line(1:10), *, iostat=ioerr) ver
      HD%ver = int(ver*100)
      if (ioerr .ne. 0) msg = 'read RINEX VERSION error.'
      if (HD%ver .lt. 100 .or. HD%ver .ge. 420) &
        write (msg, '(a,f4.2)') 'invalid RINEX VERSION ', HD%ver/100
!
!! site name
    else if (index(keyword, 'MARKER NAME') .ne. 0) then
      read (line, '(a4)') HD%mark
!
!! receiver number and type
    else if (index(keyword, 'REC #') .ne. 0) then
      HD%recnum = line(1:20)
      HD%rectyp = line(21:40)
!
!! antenna number and type
    else if (index(keyword, 'ANT #') .ne. 0) then
      HD%antnum = line(1:20)
      HD%anttyp = line(21:40)
      if (HD%anttyp(17:20) .eq. '    ') HD%anttyp(17:20) = 'NONE'
!
!! station coordinate
    else if (index(keyword, 'APPROX POSITION') .ne. 0) then
      read (line, '(3f14.4)', iostat=ioerr) HD%x, HD%y, HD%z
      if (ioerr .ne. 0) msg = 'read APPROX POSITION error'
!
!! antenna offset
    else if (index(keyword, 'ANTENNA: DELTA') .ne. 0) then
      read (line, '(3f14.4)', iostat=ioerr) HD%h, HD%e, HD%n
      if (ioerr .ne. 0) msg = 'read ANTENNA: DELTA error.'
!
!! wavelength fact
    else if (index(keyword, 'WAVELENGTH FACT') .ne. 0) then
      read (line, '(2i6)', iostat=ioerr) HD%fact1, HD%fact2
      if (ioerr .ne. 0) msg = 'read WAVELENGTH FACT error'
!
!! type of obwrvations (2.10)
    else if (index(keyword, 'TYPES OF OBSERV') .ne. 0) then
      read (line, '(i6,9(4x,a2))', iostat=ioerr) HD%nobstyp, (HD%obstyp(i), i=1, min(HD%nobstyp, 9))
      if (ioerr .ne. 0) msg = 'read TYPES OF OBSERV error'
      ii = 1
      do while (HD%nobstyp .gt. 9*ii)
        read (lfn, '(6x,9(4x,a2))', iostat=ioerr) (HD%obstyp(9*ii + i), i=1, min(HD%nobstyp - 9*ii, 9))
        if (ioerr .ne. 0) msg = 'read TYPES OF OBSERV error'
        ii = ii + 1
      end do
!
!! type of obwrvations (3.03)
    else if (index(keyword, 'OBS TYPES') .ne. 0) then
      i0 = index(GNSS_PRIO, line(1:1))
      if (i0 .eq. 0) cycle
      read (line, '(3x,i3,13(1x,a3))', iostat=ioerr) HD%nobstyp3_sys(i0), &
        (HD%cobstyp3_sys(i, i0), i=1, min(HD%nobstyp3_sys(i0), 13))
      ii = 1
      do while (HD%nobstyp3_sys(i0) .gt. 13*ii)
        read (lfn, '(6x,13(1x,a3))', iostat=ioerr, end=200) &
          (HD%cobstyp3_sys(13*ii + i, i0), i=1, min(HD%nobstyp3_sys(i0) - 13*ii, 13))
        if (ioerr .ne. 0) msg = 'read TYPES OF OBSERV error'
        ii = ii + 1
      end do
!
!! whether cc2nc has already been applied
    else if (index(keyword, 'COMMENT') .ne. 0) then
      if (index(line, 'C1 forced to P1') .ne. 0) then
        HD%lc1p1 = .true.
      else if (index(line, 'C2 forced to P2') .ne. 0) then
        HD%lc2p2 = .true.
      end if
!
!! interval
    else if (index(keyword, 'INTERVAL') .ne. 0) then
      read (line, *, iostat=ioerr) HD%intv
      if (ioerr .ne. 0) msg = 'read INTERVAL error'
!
!! start time
    else if (index(keyword, 'TIME OF FIRST OBS') .ne. 0) then
      read (line, '(5i6,f12.6)', iostat=ioerr) (HD%t0(i), i=1, 5), HD%t0s
      if (ioerr .ne. 0) msg = ' read TIME OF FIRST OBS error'
!
!! stop  time
    else if (index(keyword, 'TIME OF LAST OBS') .ne. 0) then
      read (line, '(5i6,f12.6)', iostat=ioerr) (HD%t1(i), i=1, 5), HD%t1s
      if (ioerr .ne. 0) msg = ' read TIME OF LAST OBS error'
    end if
    if (len_trim(msg) .ne. 0) goto 100
  end do
!
!! error
100 continue
  write (*, '(a)') '***ERROR: '//msg
  write (*, '(a)') '   line : '//line
  ierr = 1
  return
!
!! end of file
200 continue
  ierr = 2
  return
end subroutine
