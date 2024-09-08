!
!! get_otl_args.f90
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
!! Contributor: Yingda Deng
!
subroutine get_otl_args(OTL)
  implicit none
  include 'otl.h'
  type(otls),intent(out) :: OTL ! ocean thermal loading settings
!
!! local
  integer*4 nargs, i, n, ymdhm(5), stat
  real*8 sec, mjd0, mjd
  character*256 :: line
  character*100 :: buffer(6)
!
!! functions called
  integer*4 get_valid_unit
!
!! read command arguments
  nargs = iargc()
  if (nargs .lt. 14) then
    write (*, '(a)') 'Usage: oceanload -b lat -l lon -h hgt -s ts -e te -i ti -o outfile'
    call exit(4)
  endif
  i = 1
  do while (i .le. nargs)
    call getarg(i, line)
    i = i + 1
    if (line(1:2) .eq. '-b') then
      call getarg(i, line)
      read (line, *) OTL%plat
    elseif (line(1:2) .eq. '-l') then
      call getarg(i, line)
      read (line, *) OTL%plon
    elseif (line(1:2) .eq. '-h') then
      call getarg(i, line)
      read (line, *) OTL%phgt
    elseif (line(1:2) .eq. '-s') then
      call getarg(i, line)
      call split_string(.false., line, ' ', ' ', '/', n, buffer)
      read(buffer(1:3),*) ymdhm(1), ymdhm(2),ymdhm(3)
      i=i+1
      call getarg(i, line)
      call split_string(.false., line, ' ', ' ', ':', n, buffer)
      read(buffer(1:3),*) ymdhm(4), ymdhm(5),sec
      call CAL2JD(ymdhm(1), ymdhm(2), ymdhm(3), mjd0, mjd, stat)
      OTL%mjd0=INT(mjd)
      OTL%sod0=ymdhm(4)*3600+ymdhm(5)*60+sec
    elseif (line(1:2) .eq. '-e') then
      call getarg(i, line)
      call split_string(.false., line, ' ', ' ', '/', n, buffer)
      read(buffer(1:3),*) ymdhm(1), ymdhm(2),ymdhm(3)
      i=i+1
      call getarg(i, line)
      call split_string(.false., line, ' ', ' ', ':', n, buffer)
      read(buffer(1:3),*) ymdhm(4), ymdhm(5),sec
      call CAL2JD(ymdhm(1), ymdhm(2), ymdhm(3), mjd0, mjd, stat)
      OTL%mjd1=INT(mjd)
      OTL%sod1=ymdhm(4)*3600+ymdhm(5)*60+sec
    elseif (line(1:2) .eq. '-i') then
      call getarg(i, line)
      read (line, *) OTL%interval    
    elseif (line(1:2) .eq. '-o') then
      call getarg(i, line)
      OTL%otlname=line
    endif
    i = i + 1
  enddo

end

