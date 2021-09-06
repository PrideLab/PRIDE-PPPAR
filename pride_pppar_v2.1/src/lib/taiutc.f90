!
!! taiutc.f90
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
!! Contributor: Maorong Ge
!! 
!!
!!
!! purpose : return the leap seconds at mjd
!!
real*8 function taiutc(mjd)
  implicit none
  include '../header/const.h'

  integer*4 mjd
!
!! local
  integer*4 lun, ios, i
  logical*1 first
  integer*4 jdt(50), ls(50), nls
  character*80 line
!
!! funciton called
  integer*4 get_valid_unit

  data first/.true./
!
!! save
  save nls, jdt, ls, first

  if (first) then
    lun = get_valid_unit(10)
    open (lun, file='leap.sec', status='old', iostat=ios)
    if (ios .ne. 0) then
      write (*, *) '***ERROR(taiutc): open leap.sec, ', ios
      call exit(1)
    endif
    first = .false.
!
!! first line
    do while (index(line, '+leap sec') .eq. 0)
      read (lun, '(a)', end=100) line
    enddo
!
!! read file
    nls = 0
    read (lun, '(a)', end=100) line
    do while (index(line, '-leap sec') .eq. 0)
      nls = nls + 1
      if (nls .gt. 50) then
        write (*, *) '***ERROR(taiutc): more than 50 LS records, change dimmension'
        call exit(1)
      endif
      read (line, *, err=200) jdt(nls), ls(nls)
      read (lun, '(a)', end=100) line
    enddo
    close (lun)
  endif

  if (mjd .le. jdt(1)) then
    write (*, *) '***ERROR(taiutc): epoch before table start,', mjd, jdt(1)
    call exit(1)
  endif
  if (mjd .gt. jdt(nls)) then
    write (*, *) '***ERROR(taiutc): epoch after table end,', mjd, jdt(nls)
    call exit(1)
  endif

  do i = 1, nls - 1
    if (mjd .gt. jdt(i) .and. mjd .le. jdt(i + 1)) then
      taiutc = ls(i)
      return
    endif
  enddo

  return
100 continue
  write (*, *) '***ERROR(taiutc): sinex endline `-leap sec` not found'
  call exit(1)

200 continue
  write (*, '(a/a)') '***ERROR(taiutc): read leap.sec error,', line
  call exit(1)
end
