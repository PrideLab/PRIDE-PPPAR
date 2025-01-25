!
!! read_kinpos.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
!!
!!
!!
!! purpose  : read kinematic position file
!! parameter:
!!    input : SITE%  -- station struct
!!            fjd0,fjd1 -- requested time period
!!    output: x,y,z  -- kinematic position in m
!
subroutine read_kinpos(SITE, fjd0, fjd1, x, y, z, dwnd)
  implicit none
  include '../header/const.h'
  include '../header/station.h'

  real*8 fjd0, fjd1, x, y, z, dwnd
  type(station) SITE
!
!! local
  logical*1 lopened, lexist
  integer*4 i, jdx, iepo, ierr
  real*8 fjdx, sodx, xt(3)
  character*256 line
!
!! function called
  integer*4 get_valid_unit
  real*8 timdif

!
!! check and open kin file
  inquire (file=SITE%kinfil, opened=lopened)
  if (.not. lopened) then
    SITE%ikin = get_valid_unit(10)
    open (SITE%ikin, file=SITE%kinfil, status='old', iostat=ierr)
    if (ierr .eq. 0) then
      write (*, '(2a)') '%%%MESSAGE(read_kinpos): kinematic read ', trim(SITE%kinfil)
    else
      write (*, '(2a)') '***ERROR(read_kinpos): open file ', trim(SITE%kinfil)
      call exit(1)
    end if
!
!! read header
    do while (.true.)
      read (SITE%ikin, '(a)', end=100) line
      if (index(line, 'END OF HEADER') .ne. 0) exit
    end do
  end if
!
!! initialize
  x = 1.d0
  y = 1.d0
  z = 1.d0
  if (SITE%ikin .eq. 0) return
!
!! read kinematic position
  line = ' '
  iepo = 0
  do while (.true.)
    read (SITE%ikin, '(a)', end=100) line
!
!! skip the header line but not the bad records
    if (line(1:1) .eq. '*') cycle
    read (line(:15), *, err=200) jdx, sodx
    read (line(18:), *, err=200) (xt(i), i = 1, 3)
    fjdx = jdx + sodx/864.d2
    if ((fjdx - fjd0) * 864.d2 .gt. - dwnd) then
      if ((fjdx - fjd1) * 864.d2 .lt. - dwnd) then
        x = x + xt(1)
        y = y + xt(2)
        z = z + xt(3)
        iepo = iepo + 1
      else
        if ((iepo .eq. 0) .and. &
            (fjdx - fjd1) * 864.d2 .lt. + dwnd) then
          x = xt(1)
          y = xt(2)
          z = xt(3)
          return
        end if
        backspace SITE%ikin
        goto 110
      end if
    end if
  end do

100 write (*, '(2a)') '###WARNING(read_kinpos): end of file ', trim(SITE%kinfil)
  if (SITE%ikin .ne. 0) close (SITE%ikin)
  SITE%ikin = 0

110 continue
  if (iepo .eq. 0) iepo = 1
  x = x / iepo
  y = y / iepo
  z = z / iepo
  return

200 write (*, '(2a)') '***ERROR(read_kinpos): read file ', trim(SITE%kinfil)
  call exit(1)
end subroutine
