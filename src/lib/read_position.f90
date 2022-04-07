!
!! read_position.f90
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
!! purpose  : read position for static stations
!! parameter:
!!    input : flnpos -- position file
!!            name   -- station name
!!    output: x      -- requested positions in km
!
subroutine read_position(flnpos, name, skd, x, dx0)
  implicit none
  include '../header/const.h'

  real*8 x(3), dx0(3)
  character*4 name, name1
  character*(*) flnpos, skd
!
!! local
  logical*1 lfirst
  integer*4 i, k, lfnsit, lfnpos, ierr
  real*8 v(3), tref, x1(3), dx1(3)
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
  x = 0.d0
!
!! first search pos file
  if (lfnpos .ne. 0 .and. skd(1:1) .ne. 'F') then
    rewind lfnpos
!
!! read header
    do while (.true.)
      read (lfnpos, '(a)', end=100) line
      if (index(line, 'END OF HEADER') .ne. 0) exit
    enddo
    do while (.true.)
      read (lfnpos, '(a)', end=100) line
      if (line(1:1) .eq. '*') cycle
      name1 = upper_string(line(2:5))
      if (name1 .eq. upper_string(name)) then
        read (line(18:), *) (x(i), i=1, 3)
        x(1:3) = x(1:3)*1.d-3
        exit
      endif
    enddo
  endif
100 continue
  if (x(1) .ne. 0.d0 .or. x(2) .ne. 0.d0 .or. x(3) .ne. 0.d0) return
!
!! find request position in sit.xyz
  if (lfnsit .ne. 0) then
    rewind lfnsit
    do while (.true.)
      read (lfnsit, '(a)', end=200) line
      if (line(1:1) .ne. ' ' .or. len_trim(line) .eq. 0) cycle
      if (skd(1:1) .eq. 'F') then
        read (line, *, err=300) name1, (x1(i), i=1, 3), (dx1(i), i=1, 3)
      else
        read (line, *) name1, (x1(i), i=1, 3)
      endif
      name1 = upper_string(name1)
      if (name1 .eq. upper_string(name)) then
        x(1:3) = x1(1:3)*1.d-3
        if (skd(1:1).eq.'F') dx0(1:3)=dx1(1:3)
      endif
    enddo
  endif
200 continue
  if (x(1) .ne. 0.d0 .or. x(2) .ne. 0.d0 .or. x(3) .ne. 0.d0) return
!
!! position not found
  write (*, '(2a)') '###WARNING(read_position): position not found ', name
  do i = 1, 3
    x(i) = 1.d0
  enddo

  return
300 write(*,'(a)') '***ERROR(read_position): read sit.xyz'
  call exit(1)
end
