!
!! mergesp3.f90
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
!! Contributor: Maorong Ge, Jianghui Geng
!! 
!!
!!
!
!  FUNCTIONS:
!  mergesp3 - mrege oribit.

program mergesp3
  implicit none

  integer*4 nargs, lfn1, lfn2, lfn3, merlfn
  integer*4 ierr1, ierr2, ierr3
  character*256 sp3fil1      ! last day
  character*256 sp3fil2      ! processing day
  character*256 sp3fil3      ! next day
  character*256 mersp3fil
  character*256 line1, line2, line3

!! functions called
  integer*4 get_valid_unit
!! read command arguments
  nargs = iargc()
  if (nargs .eq. 0) then
    write (*, '(a)') 'mergesp3 command arguments error!'
    call exit(4)
  endif

  sp3fil1 = ' '
  sp3fil2 = ' '
  sp3fil3 = ' '
  mersp3fil = ' '

!! get arguements
  call get_sp3orb_args(sp3fil1, sp3fil2, sp3fil3, mersp3fil)

  lfn1 = get_valid_unit(10)
  lfn2 = get_valid_unit(20)
  lfn3 = get_valid_unit(30)
  merlfn = get_valid_unit(40)
  open (lfn1, file=sp3fil1, status='old', iostat=ierr1)
  open (lfn2, file=sp3fil2, status='old', iostat=ierr2)
  open (lfn3, file=sp3fil3, status='old', iostat=ierr3)
  open (merlfn, file=mersp3fil)

  if (ierr1 .ne. 0) then
    write (*, '(2a)') '***ERROR(mergesp3): open file ', trim(sp3fil1)
    call exit(1)
  endif

  if (ierr2 .ne. 0) then
    write (*, '(2a)') '***ERROR(mergesp3): open file ', trim(sp3fil2)
    call exit(1)
  endif

  if (ierr3 .ne. 0) then
    write (*, '(2a)') '***ERROR(mergesp3): open file ', trim(sp3fil3)
    call exit(1)
  endif

  line2 = ' '
  do while (index(line2, 'EOF') .eq. 0)
    read (lfn2, '(a)') line2
    if (line2(1:1) == '*') then
      line1 = ' '
      line3 = ' '
      !! add last day orbit data
      do while (index(line1, 'EOF') .eq. 0)
        read (lfn1, '(a)') line1
        if (line1(1:1) == '*') then
          write (merlfn, '(a80)') line1
          read (lfn1, '(a)') line1
          do while (index(line1, 'EOF') .eq. 0)
            write (merlfn, '(a80)') line1
            read (lfn1, '(a)') line1
          enddo
        endif
      enddo

      !! add processing day orbit data
      do while (index(line2, 'EOF') .eq. 0)
        write (merlfn, '(a80)') line2
        read (lfn2, '(a)') line2
      enddo

      !! add next day orbit data
      do while (index(line3, 'EOF') .eq. 0)
        read (lfn3, '(a)') line3
        if (line3(1:1) == '*') then
          write (merlfn, '(a80)') line3
          read (lfn3, '(a)') line3
          do while (index(line3, 'EOF') .eq. 0)
            write (merlfn, '(a80)') line3
            read (lfn3, '(a)') line3
          enddo
        endif
      enddo
    endif
    write (merlfn, '(a80)') line2
  enddo
  close (merlfn)
  close (lfn1)
  close (lfn2)
  close (lfn3)
end program mergesp3
