!
!! rdorbh.f90
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
!! Contributor: Maorong Ge
!! 
!!
!!
!! purpose   : read header of PANDA orbit file
!! parameters: orbfil -- orbit file name
!!             iunit  -- unit of file
!!             OH     -- header part 1 (oi control data)
!
subroutine rdorbh(orbfil, iunit, OH)
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'

  integer*4 iunit
  character*(*) orbfil
  type(orbhdr) OH
!
!! local
  integer*4 ierr
!
!! function called
  integer*4 get_valid_unit

  iunit = get_valid_unit(10)
  open (iunit, file=orbfil, status='old', form='unformatted', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(rdorbh): open orbit ', trim(orbfil)
    call exit(1)
  endif
!
!! read header of the orbit file
  read (iunit, iostat=ierr) OH
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(rdorbh): read orbit ', trim(orbfil)
    call exit(1)
  endif

  return
end
