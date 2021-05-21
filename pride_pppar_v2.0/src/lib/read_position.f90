!
!! read_position.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License (version 3) as
!!    published by the Free Software Foundation.
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
!! purpose  : read position for static stations
!! parameter:
!!    input : flnpos -- position file
!!            name   -- station name
!!            jd,sod -- time tag
!!            seslen -- session length
!!    output: x      -- requested positions in km
!!            flag   -- read from sit.xyz (0) or pos-file (1)
!
subroutine read_position(flnpos, name, jd, sod, seslen, x, flag)
  implicit none
  include '../header/const.h'

  integer*4 jd, flag
  real*8 sod, seslen, x(3), x1(3)
  character*4 name, name1
  character*(*) flnpos
!
!! local
  logical*1 lfirst
  integer*4 i, k, lfnsit, lfnpos, iyear, idoy, ierr
  real*8 v(3), tref
  character*256 line
!
!! function called
  integer*4 get_valid_unit
  character*4 upper_string

  data lfirst/.true./
  save lfirst, lfnsit, lfnpos

  if (lfirst) then
    lfirst = .false.
    lfnpos = get_valid_unit(10)
    open (lfnpos, file=flnpos, status='old', iostat=ierr)
    if (ierr .ne. 0) lfnpos = 0
    lfnsit = get_valid_unit(10)
    open (lfnsit, file='sit.xyz', status='old', iostat=ierr)
    if (ierr .ne. 0) lfnsit = 0
  endif
  flag = -1
  x = 0.d0
!
!! first search pos file
  if (lfnpos .ne. 0) then
    rewind lfnpos
    do while (.true.)
      read (lfnpos, '(a)', end=100) line
      name1 = upper_string(line(1:4))
      name = upper_string(name)
      if (name1 .eq. name) then
        read (line(17:61), *) (x(i), i=1, 3)
        x(1:3) = x(1:3)*1.d-3
        exit
      endif
    enddo
  endif
100 continue
  if (x(1) .ne. 0.d0 .or. x(2) .ne. 0.d0 .or. x(3) .ne. 0.d0) then
    flag = 1
    return
  endif
!
!! find request position in sit.xyz
  if (lfnsit .ne. 0) then
    rewind lfnsit
    do while (.true.)
      read (lfnsit, '(a)', end=200) line
      if (line(1:1) .ne. ' ' .or. len_trim(line) .eq. 0) cycle
      read (line, *) name1, (x1(i), i=1, 3)
      name1 = upper_string(name1)
      name = upper_string(name)
      if (name1 .eq. name) then
        x(1:3) = x1(1:3)*1.d-3
        if (jd .ne. 0) then
          call mjd2doy(jd, iyear, idoy)
        endif
      endif
    enddo
  endif
200 continue
  if (x(1) .ne. 0.d0 .or. x(2) .ne. 0.d0 .or. x(3) .ne. 0.d0) then
    flag = 0
    return
  endif
!
!! position not found
  write (*, '(2a)') '###WARNING(read_position): position not found ', name
  do i = 1, 3
    x(i) = 1.d0
  enddo

  return
end
